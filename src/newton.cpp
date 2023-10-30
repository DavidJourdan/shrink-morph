#include "newton.h"

#include "functions.h"
#include "parameterization.h"
#include "simulation_utils.h"
#include "timer.h"

#include <TinyAD/Utils/NewtonDecrement.hh>

#ifdef USE_PARDISO
#include "Eigen/PardisoSupport"
using LLTSolver = Eigen::PardisoLLT<Eigen::SparseMatrix<double>>;
using LUSolver = Eigen::PardisoLU<Eigen::SparseMatrix<double>>;
#else
using LLTSolver = Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>;
using LUSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>>;
#endif // USE_PARDISO

using namespace geometrycentral::surface;

template <class Func, class Solver>
void newton(Eigen::VectorXd& x,
            Func& func,
            Solver& solver,
            int max_iters,
            double lim,
            bool verbose,
            std::vector<int> fixedIdx,
            const std::function<void(const Eigen::VectorXd&)>& callback)
{
  Timer timer("Newton", !verbose);

  if(verbose)
    std::cout << "Initial energy: " << func.eval(x) << std::endl;

  Eigen::SparseMatrix<double> P = projectionMatrix(fixedIdx, x.size());

  for(int i = 0; i < max_iters; ++i)
  {
    auto [f, g, H] = func.eval_with_derivatives(x);

    for(int j = 0; j < H.cols(); ++j)
      H.coeffRef(j, j) += 1e-10;

    // restrict H and g to free variables
    H = (P * H * P.transpose()).eval();
    g = P * g;

    // Newton direction
    if(i == 0)
      solver.compute(H);
    else
      solver.factorize(H);

    bool exact = true;
    if(solver.info() != Eigen::Success)
    {
      exact = false;
      auto [f, g, H_proj] = func.eval_with_hessian_proj(x);
      H_proj = (P * H_proj * P.transpose()).eval();

      H = 0.9 * H + 0.1 * H_proj;
      solver.factorize(H);
      if(solver.info() != Eigen::Success)
        solver.factorize(H_proj);
    }

    Eigen::VectorXd d = -solver.solve(g);

    if(verbose)
    {
      if(exact)
        std::cout << "Decrement in iteration " << i << ": " << TinyAD::newton_decrement(d, g)
                  << "\tFactorization = Exact\n";
      else
        std::cout << "Decrement in iteration " << i << ": " << TinyAD::newton_decrement(d, g)
                  << "\tFactorization = Project\n";
    }

    d = P.transpose() * d;
    g = P.transpose() * g;

    double s = lineSearch(x, d, f, g, func);
    if(s < 0)
      break;
    x += s * d;

    if(TinyAD::newton_decrement(d, g) < lim && exact)
      break;

    callback(x);
  }
  if(verbose)
    std::cout << "Final energy: " << func.eval(x) << "\n";
}

template <class Func>
void newton(IntrinsicGeometryInterface& geometry,
            Eigen::MatrixXd& V,
            Func& func,
            int max_iters,
            double lim,
            bool verbose,
            std::vector<int> fixedIdx,
            const std::function<void(const Eigen::VectorXd&)>& callback)
{
  // Assemble inital x vector
  geometry.requireVertexIndices();
  Eigen::VectorXd x = func.x_from_data([&](Vertex v) { return V.row(geometry.vertexIndices[v]); });

  LLTSolver solver;

  // Newton algorithm
  newton(x, func, solver, max_iters, lim, verbose, fixedIdx, callback);

  func.x_to_data(x, [&](Vertex v, const auto& row) { V.row(geometry.vertexIndices[v]) = row; });
}

void sparse_gauss_newton(IntrinsicGeometryInterface& geometry,
                         Eigen::MatrixXd& V,
                         const Eigen::VectorXd& xTarget,
                         const FaceData<Eigen::Matrix2d>& MrInv,
                         FaceData<double>& theta1,
                         VertexData<double>& theta2,
                         TinyAD::ScalarFunction<1, double, Eigen::Index>& adjointFunc,
                         std::vector<int>& fixedIdx,
                         int max_iters,
                         double lim,
                         double wM,
                         double wL,
                         double E1,
                         double lambda1,
                         double lambda2,
                         double deltaLambda,
                         double thickness,
                         const std::function<void(const Eigen::VectorXd&)>& callback)
{
  geometry.requireCotanLaplacian();
  geometry.requireVertexLumpedMassMatrix();
  geometry.requireFaceAreas();
  geometry.requireVertexIndices();

  SurfaceMesh& mesh = geometry.mesh;

  // build vector of Voronoi areas for the distance metric
  Eigen::VectorXd masses = Eigen::VectorXd::Zero(V.size());

  // divide M by the total mesh area
  double totalArea = 0;
  for(Face f: mesh.faces())
  {
    for(Vertex v: f.adjacentVertices())
    {
      masses(3 * geometry.vertexIndices[v]) += geometry.faceAreas[f] / 3.;
      masses(3 * geometry.vertexIndices[v] + 1) += geometry.faceAreas[f] / 3.;
      masses(3 * geometry.vertexIndices[v] + 2) += geometry.faceAreas[f] / 3.;
    }
    totalArea += geometry.faceAreas[f];
  }
  masses /= totalArea;

  Eigen::SparseMatrix<double> L = geometry.cotanLaplacian;

  // Mass matrix theta
  Eigen::SparseMatrix<double> M_theta(V.rows(), V.rows());
  M_theta.reserve(V.rows());
  for(int i = 0; i < V.rows(); ++i)
    M_theta.insert(i, i) = totalArea * masses(3 * i);

  Eigen::VectorXd theta = theta2.toVector();
  Eigen::VectorXd x = xTarget;

  LLTSolver adjointSolver;

  auto distance = [&](const Eigen::VectorXd& th) {
    theta2.fromVector(th);
    auto simFunc = simulationFunction(geometry, MrInv, theta1, theta2, E1, lambda1, lambda2, deltaLambda, thickness);
    newton(x, simFunc, adjointSolver, 1000, lim, false, fixedIdx);

    return (x - xTarget).dot(masses.cwiseProduct(x - xTarget)) + wM * th.dot(M_theta * th) + wL * th.dot(L * th);
  };

  // Build matrix P
  Eigen::SparseMatrix<double> P = projectionMatrix(fixedIdx, x.size());

  // Hessian matrix H
  Eigen::VectorXd X(V.size() + theta.size());
  X.head(V.size()) = x;
  X.tail(theta.size()) = theta;
  Eigen::SparseMatrix<double> H = adjointFunc.eval_hessian(X);

  // Build HGN matrix
  Eigen::SparseMatrix<double> HGN = buildHGN(2 * masses, P, 2 * wM * M_theta + 2 * wL * L, H);

  auto distanceGrad = [&](const Eigen::VectorXd& th) -> Eigen::VectorXd {
    Eigen::VectorXd X(V.size() + th.size());
    X.head(V.size()) = x;
    X.tail(th.size()) = th;
    H = adjointFunc.eval_hessian(X);

    for(int j = 0; j < V.size(); ++j)
      H.coeffRef(j, j) += 1e-10;

    Eigen::SparseMatrix<double> A = (P * H.block(0, 0, V.size(), V.size()) * P.transpose()).eval();

    adjointSolver.factorize(A);
    if(adjointSolver.info() != Eigen::Success)
    {
      auto [f, g, A_proj] = adjointFunc.eval_with_hessian_proj(X);
      A_proj = (P * A_proj.block(0, 0, V.size(), V.size()) * P.transpose()).eval();

      A = 0.9 * A + 0.1 * A_proj;
      adjointSolver.factorize(A);
      if(adjointSolver.info() != Eigen::Success)
        adjointSolver.factorize(A_proj);
    }

    Eigen::VectorXd b = P * masses.cwiseProduct(x - xTarget);
    Eigen::VectorXd dir = adjointSolver.solve(b);
    if(adjointSolver.info() != Eigen::Success)
      std::cout << "Solver error\n";

    dir = P.transpose() * dir;

    return -2 * H.block(V.size(), 0, th.size(), V.size()) * dir + 2 * wM * M_theta * th + 2 * wL * L * th;
  };

  std::cout << "Initial energy: " << distance(theta) << std::endl;

  LUSolver solver;

  for(int i = 0; i < max_iters; ++i)
  {
    double f = distance(theta);
    Eigen::VectorXd g = distanceGrad(theta);

    Eigen::VectorXd b(2 * x.size() - 2 * fixedIdx.size() + theta.size());
    b.setZero();
    b.segment(x.size() - fixedIdx.size(), theta.size()) = -g;

    // Update HGN
    updateHGN(HGN, P, H);

    if(i == 0)
      solver.compute(HGN);
    else
      solver.factorize(HGN);

    if(solver.info() != Eigen::Success)
    {
      std::cout << "Solver error\n";
      return;
    }

    Eigen::VectorXd d = solver.solve(b);
    Eigen::VectorXd deltaTheta = d.segment(x.size() - fixedIdx.size(), theta.size());
    Eigen::VectorXd deltaX = d.segment(0, x.size() - fixedIdx.size());
    deltaX = P.transpose() * deltaX;

    // LINE SEARCH
    Eigen::VectorXd x_old = x;
    double s = lineSearch(theta, deltaTheta, f, g, distance, [&](double s) { x = x_old + s * deltaX; });
    if(s < 0)
    {
      std::cout << "Line search failed\n";
      break;
    }
    theta += s * deltaTheta;

    std::cout << "Decrement in iteration " << i << ": " << TinyAD::newton_decrement(deltaTheta, g)
              << "\tDistance: " << (x - xTarget).dot(masses.cwiseProduct(x - xTarget)) << "\tStep size: " << s
              << std::endl;
    if(TinyAD::newton_decrement(deltaTheta, g) < lim || solver.info() != Eigen::Success)
      break;

    callback(x);
  }

  std::cout << "Final energy: " << distance(theta) << "\n";

  for(int i = 0; i < V.rows(); ++i)
    for(int j = 0; j < 3; ++j)
      V(i, j) = x(3 * i + j);
  theta2.fromVector(theta);
}

template void newton<>(IntrinsicGeometryInterface&,
                       Eigen::MatrixXd&,
                       TinyAD::ScalarFunction<2, double, Vertex>&,
                       int,
                       double,
                       bool,
                       std::vector<int>,
                       const std::function<void(const Eigen::VectorXd&)>&);

template void newton<>(IntrinsicGeometryInterface&,
                       Eigen::MatrixXd&,
                       TinyAD::ScalarFunction<3, double, Vertex>&,
                       int,
                       double,
                       bool,
                       std::vector<int>,
                       const std::function<void(const Eigen::VectorXd&)>&);
