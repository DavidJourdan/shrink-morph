#include "LocalGlobalSolver.h"

#include <geometrycentral/utilities/vector2.h>
#include <igl/cat.h>
#include <igl/columnize.h>
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/project_isometrically_to_plane.h>
#include <igl/repdiag.h>

LocalGlobalSolver::LocalGlobalSolver(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
  init(V, F);
}

void LocalGlobalSolver::init(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
  using namespace Eigen;

  // number of vertices
  nV = V.rows();
  nF = F.rows();

  _F = F;
  s1.resize(nF);
  s2.resize(nF);
  stressX.resize(nF, 2);
  stressY.resize(nF, 2);

  MatrixX2d plane_V;
  MatrixX3i plane_F;
  SparseMatrix<double> ref_map, ref_map_dim;

  igl::project_isometrically_to_plane(V, F, plane_V, plane_F, ref_map);
  igl::repdiag(ref_map, 2, ref_map_dim);

  int nF = F.rows();
  _Ainv.resize(2, 2 * nF);
  for(int i = 0; i < nF; ++i)
  {
    Matrix2d A;
    A.col(0) = plane_V.row(nF + i) - plane_V.row(i);
    A.col(1) = plane_V.row(2 * nF + i) - plane_V.row(i);
    _Ainv.block<2, 2>(0, 2 * i) = A.inverse();
  }

  buildRHS(plane_V, plane_F);
  _K = (ref_map_dim * _K).eval();
  _K.prune([&](const int row, const int col, const double val) { return row != 0 && row != nV; });

  assert(_K.rows() == nV * 2);

  SparseMatrix<double> L; // cotan Laplacian
  igl::cotmatrix(V, F, L);
  SparseMatrix<double> Q = (-L).eval();
  Q.prune([](const int row, const int col, const double val) { return row > 0 && col > 0; });
  Q.reserve(1);
  Q.insert(0, 0) = 1;

  _solver = std::make_unique<geometrycentral::PositiveDefiniteSolver<double>>(Q);
}

void LocalGlobalSolver::solveOneStep(Eigen::Ref<Eigen::MatrixX2d> U, double sMin, double sMax)
{
  using namespace Eigen;

  Matrix<double, 2, -1> R(2, 2 * nF);
#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads() - 1)
  for(int i = 0; i < nF; ++i)
  {
    Matrix2d B;
    B.col(0) = U.row(_F(i, 1)) - U.row(_F(i, 0));
    B.col(1) = U.row(_F(i, 2)) - U.row(_F(i, 0));

    R.block<2, 2>(0, 2 * i) = project(B * _Ainv.block<2, 2>(0, 2 * i), sMin, sMax, i);
  }

  VectorXd Rcol;
  igl::columnize(R, nF, 2, Rcol);
  VectorXd Bcol = -_K * Rcol;

  U.col(0) = _solver->solve(-Bcol.segment(0, nV));
  U.col(1) = _solver->solve(-Bcol.segment(nV, nV));
}

void LocalGlobalSolver::solve(Eigen::Ref<Eigen::MatrixX2d> U, double sMin, double sMax, int nbIter)
{
  using namespace Eigen;

  if(nbIter > 0)
    for(int i = 0; i < nbIter; ++i)
      solveOneStep(U, sMin, sMax);
  else
  {
    VectorXd s1_prev = VectorXd::Constant(s1.size(), sMin);
    VectorXd s2_prev = VectorXd::Constant(s2.size(), sMax);
    solveOneStep(U, sMin, sMax);
    while(std::max((s1 - s1_prev).norm() / s1.size(), (s2 - s2_prev).norm() / s2.size()) > 1e-8)
    {
      s1_prev = s1;
      s2_prev = s2;
      solveOneStep(U, sMin, sMax);
    }
  }
}

Eigen::Matrix2d LocalGlobalSolver::project(const Eigen::Matrix2d& M, double sMin, double sMax, int i)
{
  using namespace Eigen;
  JacobiSVD<Matrix2d> svd;
  svd.compute(M, ComputeFullU | ComputeFullV);
  Matrix2d S = svd.singularValues().asDiagonal();

  // save data
  s1(i) = S(0, 0);
  s2(i) = S(1, 1);
  stressX.row(i) = svd.matrixU().col(0);
  stressY.row(i) = svd.matrixU().col(1);

  S(1, 1) = sMin;
  S(0, 0) = sMax;
  return svd.matrixU() * S * svd.matrixV().transpose();
}

void LocalGlobalSolver::buildRHS(const Eigen::Ref<Eigen::MatrixX2d> V, const Eigen::Ref<Eigen::MatrixX3i> F)
{
  using namespace Eigen;

  SparseMatrix<double> KX, KY;
  buildLinearBlock(V, F, 0, KX);
  buildLinearBlock(V, F, 1, KY);
  _K = igl::cat(2, igl::repdiag(KX, 2), igl::repdiag(KY, 2));
}

void LocalGlobalSolver::buildLinearBlock(const Eigen::Ref<Eigen::MatrixXd> V,
                                         const Eigen::Ref<Eigen::MatrixXi> F,
                                         const int d,
                                         Eigen::SparseMatrix<double>& Kd)
{
  using namespace Eigen;

  // Temporary output
  std::vector<Triplet<double>> Kd_IJV;
  Kd_IJV.reserve(3 * 2 * F.rows());

  // gather cotangent weights
  MatrixXd C;
  igl::cotmatrix_entries(V, F, C);
  // should have weights for each edge
  assert(C.cols() == 3);
  // loop over elements
  for(int i = 0; i < F.rows(); ++i)
  {
    // loop over edges of element
    for(int j = 0; j < 3; ++j)
    {
      int source = F(i, (j + 1) % 3);
      int dest = F(i, (j + 2) % 3);
      double v = C(i, j) * (V(source, d) - V(dest, d));
      Kd_IJV.emplace_back(source, i, v);
      Kd_IJV.emplace_back(dest, i, -v);
    }
  }
  Kd.resize(V.rows(), F.rows());
  Kd.setFromTriplets(Kd_IJV.begin(), Kd_IJV.end());
}

Eigen::VectorXd LocalGlobalSolver::stretchAngles()
{
  using namespace geometrycentral;

  Eigen::VectorXd angles(nF);
  for(int i = 0; i < nF; ++i)
  {
    Vector2 dir{stressX(i, 0), stressX(i, 1)};
    angles(i) = dir.pow(2).arg() / 2;
  }
  return angles;
}
