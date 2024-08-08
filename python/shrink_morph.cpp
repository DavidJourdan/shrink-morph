#include "functions.h"
#include "newton.h"
#include "parameterization.h"
#include "path_extraction.h"
#include "simulation_utils.h"
#include "stretch_angles.h"
#include "timer.h"

#include <Eigen/Dense>
#include <igl/loop.h>
#include <igl/readOBJ.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <string>
#include <tuple>

namespace nb = nanobind;
using namespace nb::literals;

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> readFromOBJ(std::string fileName)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(!igl::readOBJ(fileName, V, F))
  {
    std::cout << "File " << fileName << " not found\n";
  }
  V.col(1) *= -1;
  return std::make_tuple(V, F);
}

// std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXi>
// subdivide(const nb::DRef<Eigen::MatrixXd>& V, const nb::DRef<Eigen::MatrixXd>& P, const nb::DRef<Eigen::MatrixXi>& F)
// {
//   Eigen::VectorXd x;
//   const auto& [NV, NP, NF, nx] = subdivide(V, P, F, x);

//   return std::make_tuple(NV, NP, NF);
// }

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXd>
subdivide(const nb::DRef<Eigen::MatrixXd>& V,
          const nb::DRef<Eigen::MatrixXd>& P,
          const nb::DRef<Eigen::MatrixXi>& F,
          const nb::DRef<Eigen::VectorXd>& x,
          double targetLength)
{
  Eigen::MatrixXi NF;
  Eigen::SparseMatrix<double> S;

  igl::loop(V.rows(), F, S, NF);

  Eigen::MatrixXd NV = S * V;
  Eigen::MatrixXd NP = S * P;
  Eigen::MatrixXd nx;
  if(x.size() == S.cols())
    nx = S * x;

  if(targetLength > 0)
  {
    // create geometry-central objects
    using namespace geometrycentral::surface;
    ManifoldSurfaceMesh mesh(NF);
    VertexPositionGeometry geometry(mesh, NV);

    double avgEdgeLength = 0;
    geometry.requireVertexPositions();
    for(Edge e: mesh.edges())
      avgEdgeLength += norm(geometry.vertexPositions[e.firstVertex()] - geometry.vertexPositions[e.secondVertex()]);
    avgEdgeLength /= mesh.nEdges();

    while(avgEdgeLength > targetLength)
    {
      Eigen::MatrixXi tempF = NF;
      igl::loop(NV.rows(), tempF, S, NF);
      NV = S * NV;
      NP = S * NP;
      if(nx.size() == S.cols())
        nx = S * nx;

      avgEdgeLength /= 2;
    }
  }
  return std::make_tuple(NV, NP, NF, nx);
}

void simulation(nb::DRef<Eigen::MatrixXd> V,
                nb::DRef<Eigen::MatrixXd> P,
                const nb::DRef<Eigen::MatrixXi>& F,
                const nb::DRef<Eigen::VectorXd>& theta2,
                double E1,
                double lambda1,
                double lambda2,
                double deltaLambda,
                double thickness,
                double width,
                int n_iter,
                double lim)
{
  Timer timer("Simulation");

  using namespace geometrycentral::surface;

  const double scaleFactor = width / (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
  V *= scaleFactor;
  P *= scaleFactor;

  // create geometry-central objects
  ManifoldSurfaceMesh mesh(F);
  VertexPositionGeometry geometry(mesh, V);

  // Find face closest to mesh center and fix its vertices
  std::vector<int> fixedIdx = findCenterFaceIndices(P, F);

  FaceData<Eigen::Matrix2d> MrInv = precomputeSimData(mesh, P, F);
  FaceData<double> theta1 = computeStretchAngles(mesh, V, F, MrInv);
  VertexData<double> _theta2(mesh, theta2);

  // Define simulation function
  auto func = simulationFunction(geometry, MrInv, theta1, _theta2, E1, lambda1, lambda2, deltaLambda, thickness);

  // ---------------- (Projected) Newton optimization ----------------
  // Assemble inital x vector
  geometry.requireVertexIndices();
  Eigen::VectorXd x = func.x_from_data([&, V](Vertex v) { return V.row(geometry.vertexIndices[v]); });

  LLTSolver solver;

  // Newton algorithm
  newton(x, func, solver, n_iter, lim, true, fixedIdx);

  V = x.reshaped<Eigen::RowMajor>(V.rows(), 3);
}

Eigen::VectorXd directionsOptimization(nb::DRef<Eigen::MatrixXd> V,
                                       nb::DRef<Eigen::MatrixXd> targetV,
                                       nb::DRef<Eigen::MatrixXd> P,
                                       const Eigen::MatrixXi& F,
                                       double E1,
                                       double lambda1,
                                       double lambda2,
                                       double deltaLambda,
                                       double thickness,
                                       double width,
                                       int n_iter,
                                       double lim,
                                       double wM,
                                       double wL)
{

  Timer timer("Directions optimization");

  using namespace geometrycentral::surface;

  const double scaleFactor = width / (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
  V *= scaleFactor;
  targetV *= scaleFactor;
  P *= scaleFactor;

  // create geometry-central objects
  ManifoldSurfaceMesh mesh(F);
  VertexPositionGeometry geometry(mesh, V);

  // Find face closest to mesh center and fix its vertices
  std::vector<int> fixedIdx = findCenterFaceIndices(P, F);

  FaceData<Eigen::Matrix2d> MrInv = precomputeSimData(mesh, P, F);
  FaceData<double> theta1 = computeStretchAngles(mesh, V, F, MrInv);
  VertexData<double> theta2(mesh, 0);

  // Define optimization function
  auto adjointFunc = adjointFunction(geometry, F, MrInv, theta1, E1, lambda1, lambda2, deltaLambda, thickness);

  // Optimize this energy function using SGN [Zehnder et al. 2021]
  V = sparse_gauss_newton(geometry, targetV, MrInv, theta1, theta2, adjointFunc, fixedIdx, n_iter, lim, wM, wL, E1,
                          lambda1, lambda2, deltaLambda, thickness);

  return theta2.toVector();
}

std::vector<Eigen::MatrixXd> generateTrajectories(const nb::DRef<Eigen::MatrixXd>& V,
                                                  const nb::DRef<Eigen::MatrixXd>& P,
                                                  const nb::DRef<Eigen::MatrixXi>& F,
                                                  const nb::DRef<Eigen::VectorXd>& theta2,
                                                  double layerHeight,
                                                  double spacing,
                                                  int nLayers)
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  // create geometry-central objects
  ManifoldSurfaceMesh mesh(F);
  VertexPositionGeometry geometry(mesh, V);

  Eigen::MatrixXd P_3D(P.rows(), 3);
  P_3D.leftCols(2) = P;
  P_3D.col(2).setZero();
  VertexPositionGeometry geometryUV(mesh, P_3D);

  FaceData<double> theta1 = computeStretchAngles(mesh, V, P, F);

  // convert face angles into vertex angles
  Eigen::VectorXd th1(V.rows());
  for(int i = 0; i < V.rows(); ++i)
  {
    double sumAngles = 0;
    int nFaces = 0;
    Vertex v = mesh.vertex(i);
    for(Face f: v.adjacentFaces())
    {
      // add face orientations in global coordinates
      if(!f.isBoundaryLoop())
      {
        sumAngles += theta1[f];
        nFaces += 1;
      }
    }
    th1(i) = sumAngles / nFaces;
  }

  std::vector<std::vector<std::vector<Vector3>>> paths =
      generatePaths(geometryUV, th1, theta2, layerHeight, nLayers, spacing);

  // convert data to the right format
  std::vector<Eigen::MatrixXd> dataArray;
  for(auto layer: paths)
  {
    for(auto path: layer)
    {
      Eigen::MatrixXd traj(path.size(), 3);
      for(int i = 0; i < path.size(); ++i)
        for(int j = 0; j < 3; ++j)
          traj(i, j) = path[i][j];
      dataArray.push_back(traj);
    }
  }
  return dataArray;
}

NB_MODULE(shrink_morph_py, m)
{
  m.def("read_from_OBJ", &readFromOBJ);
  m.def("subdivide", &subdivide, "V"_a, "P"_a, "F"_a, "x"_a, "target_edge_length"_a = 0);
  m.def("compute_SVD_data", [](const nb::DRef<Eigen::MatrixXd>& V, const nb::DRef<Eigen::MatrixXd>& P,
                               const nb::DRef<Eigen::MatrixXi>& F) { return computeSVDdata(V, P, F); });
  m.def("parameterization",
        [](const nb::DRef<Eigen::MatrixXd>& V, Eigen::MatrixXi& F, double lambda1, double lambda2, double wD,
           int n_iter, double lim) { return parameterization(V, F, lambda1, lambda2, wD, n_iter, lim); });
  m.def("reparameterization",
        [](const nb::DRef<Eigen::MatrixXd>& V, Eigen::MatrixXd& P, const nb::DRef<Eigen::MatrixXi>& F, double lambda1,
           double lambda2, double wD, int n_iter, double lim) {
          parameterization(V, P, F, lambda1, lambda2, wD, n_iter, lim);
          return P;
        });
  m.def("simulation", &simulation);
  m.def("directions_optimization", &directionsOptimization);
  m.def("generate_trajectories", &generateTrajectories);
}
