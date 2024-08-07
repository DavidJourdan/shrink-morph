#include "simulation_utils.h"

#include <Eigen/Core>
#include <igl/colon.h>
#include <igl/slice.h>
#include <igl/slice_into.h>

Eigen::SparseMatrix<double> projectionMatrix(const std::vector<int>& fixedIdx, int size)
{
  using namespace Eigen;

  VectorXi indices(size - fixedIdx.size());
  int k = 0;
  for(int i = 0; i < size; ++i)
  {
    if(!std::binary_search(fixedIdx.begin(), fixedIdx.end(), i))
    {
      indices(k) = i;
      ++k;
    }
  }

  SparseMatrix<double> P, Id(size, size);
  Id.setIdentity();

  igl::slice(Id, indices, 1, P);

  return P;
}

Eigen::SparseMatrix<double> buildHGN(const Eigen::VectorXd& masses,
                                     const Eigen::SparseMatrix<double>& P,
                                     const Eigen::SparseMatrix<double>& M_theta,
                                     const Eigen::SparseMatrix<double>& H)
{
  using namespace Eigen;

  int n = P.cols();
  int m = P.rows();
  int nTheta = M_theta.rows();

  // Mass matrix x (P * M * P')
  SparseMatrix<double> D(n, n);
  D.reserve(n);
  for(int i = 0; i < n; ++i)
    D.insert(i, i) = masses(i);
  D = (P * D * P.transpose()).eval();

  std::vector<Triplet<double>> v;
  for(int i = 0; i < D.outerSize(); i++)
    for(typename SparseMatrix<double>::InnerIterator it(D, i); it; ++it)
      v.emplace_back(it.row(), it.col(), it.value());

  // Mass matrix theta
  for(int i = 0; i < M_theta.outerSize(); i++)
    for(typename SparseMatrix<double>::InnerIterator it(M_theta, i); it; ++it)
      v.emplace_back(it.row() + m, it.col() + m, it.value());

  // P * (df / dx)' * P'
  SparseMatrix<double> A = (P * H.block(0, 0, n, n) * P.transpose()).eval();
  for(int i = 0; i < A.outerSize(); i++)
    for(typename SparseMatrix<double>::InnerIterator it(A, i); it; ++it)
    {
      v.emplace_back(it.row(), it.col() + m + nTheta, it.value());
      v.emplace_back(it.col() + m + nTheta, it.row(), it.value());
    }

  // (df / dθ)' * P'
  A = (H.block(n, 0, nTheta, n) * P.transpose()).eval();
  for(int i = 0; i < A.outerSize(); i++)
    for(typename SparseMatrix<double>::InnerIterator it(A, i); it; ++it)
    {
      v.emplace_back(it.row() + m, it.col() + m + nTheta, it.value());
      v.emplace_back(it.col() + m + nTheta, it.row() + m, it.value());
    }

  SparseMatrix<double> HGN(2 * m + nTheta, 2 * m + nTheta);
  HGN.setFromTriplets(v.begin(), v.end());

  return HGN;
}

void updateHGN(Eigen::SparseMatrix<double>& HGN,
               const Eigen::SparseMatrix<double>& P,
               const Eigen::SparseMatrix<double>& H)
{
  using namespace Eigen;

  int n = P.cols();
  int m = P.rows();
  int nTheta = H.rows() - n;

  // P * (df / dx)' * P'
  Eigen::SparseMatrix<double> A = (P * H.block(0, 0, n, n) * P.transpose()).eval();
  igl::slice_into(A, igl::colon<int>(0, m - 1), igl::colon<int>(m + nTheta, 2 * m + nTheta - 1), HGN);

  // (df / dθ)' * P'
  A = (H.block(n, 0, nTheta, n) * P.transpose()).eval();
  igl::slice_into(A, igl::colon<int>(m, m + nTheta - 1), igl::colon<int>(m + nTheta, 2 * m + nTheta - 1), HGN);

  HGN = SparseMatrix<double>(HGN.selfadjointView<Upper>());
}

std::vector<int> findCenterFaceIndices(const Eigen::MatrixXd& P, const Eigen::MatrixXi& F)
{
  int centerIdx = 0;
  double dist = (P.row(F(0, 0)) + P.row(F(0, 1)) + P.row(F(0, 2))).norm();
  for(int i = 0; i < F.rows(); ++i)
  {
    if((P.row(F(i, 0)) + P.row(F(i, 1)) + P.row(F(i, 2))).norm() < dist)
    {
      dist = (P.row(F(i, 0)) + P.row(F(i, 1)) + P.row(F(i, 2))).norm();
      centerIdx = i;
    }
  }

  // Fixed Indices
  std::vector<int> fixedIdx = {3 * F(centerIdx, 0), 3 * F(centerIdx, 0) + 1, 3 * F(centerIdx, 0) + 2,
                               3 * F(centerIdx, 1), 3 * F(centerIdx, 1) + 1, 3 * F(centerIdx, 1) + 2,
                               3 * F(centerIdx, 2), 3 * F(centerIdx, 2) + 1, 3 * F(centerIdx, 2) + 2};
  std::sort(fixedIdx.begin(), fixedIdx.end());

  return fixedIdx;
}