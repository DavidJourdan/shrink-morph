/*
 * This file is part of TinyAD and released under the MIT license.
 * Author: Patrick Schmidt
 */
#include "LocalGlobalSolver.h"
#include "constants.h"
#include "functions.h"
#include "newton.h"
#include "parameterization.h"
#include "path_extraction.h"
#include "simulation_utils.h"
#include "stretch_angles.h"
#include "timer.h"

#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <igl/loop.h>
#include <igl/readOBJ.h>
#include <thread>

int main(int argc, char* argv[])
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  std::string filename = "beetle";
  double wD = 0;
  double width = 200;
  double lim = 1e-6;
  int n_iter = 1000;
  double timeLimit = 5;
  double wM = 0.01;
  double wL = 0.01;
  double E1 = 10;
  int nFmin = 1000;

  // Material data
  double lambda1 = 0.58;
  double lambda2 = 1.08;
  double thickness = 1.218;
  double deltaLambda = 0.0226764665509417;

  if(argc >= 2)
    filename = argv[1];
  if(argc >= 3)
    wD = std::stod(argv[2]);
  if(argc >= 4)
    width = std::stod(argv[3]);
  if(argc >= 5)
    nFmin = std::stoi(argv[4]);
  if(argc >= 6)
    timeLimit = std::stod(argv[5]);

  // Load a mesh in OBJ format
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  if(!igl::readOBJ(std::string(DATA_PATH_STR) + filename + ".obj", V, F))
  {
    std::cout << "File " << DATA_PATH_STR << filename << ".obj not found\n";
    return 0;
  }
  // resize mesh to a 100mm width
  V *= width / (V.colwise().maxCoeff() - V.colwise().minCoeff()).maxCoeff();
  V.col(1) *= -1;

  // subdivide until number of mesh faces is at least nFmin
  while(F.rows() < nFmin)
  {
    Eigen::MatrixXd tempV = V;
    Eigen::MatrixXi tempF = F;
    igl::loop(tempV, tempF, V, F);
  }

  // create geometry-central objects
  ManifoldSurfaceMesh mesh(F);
  VertexPositionGeometry geometry(mesh, V);

  // Run local-global parameterization algorithm
  Timer paramTimer("Parameterization");
  LocalGlobalSolver LGsolver(V, F);
  Eigen::MatrixXd P = localGlobal(V, F, lambda1, lambda2, LGsolver);
  paramTimer.stop();

  if(wD > 0) // smoothing
  {
    std::cout << "*********\nSMOOTHING\n*********\n";
    paramTimer.start();
    auto func = parameterizationFunction(geometry, wD, lambda1, lambda2);

    // optimize the parameterization function with Newton's method
    newton(geometry, P, func, n_iter, lim);
    paramTimer.stop();
  }

  // Resize to required width
  const double scaleFactor = width / (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
  V *= scaleFactor;
  P *= scaleFactor;
  geometry.inputVertexPositions *= scaleFactor;
  geometry.refreshQuantities();

  // precompute data
  FaceData<Eigen::Matrix2d> MrInv = precomputeSimData(mesh, P, F);
  FaceData<double> theta1 = computeStretchAngles(mesh, V, F, MrInv);

  // initialize theta2 (stored on vertices for convenience and speed)
  VertexData<double> theta2(mesh, 0);

  // Find face closest to mesh center and fix its vertices
  std::vector<int> fixedIdx = findCenterFaceIndices(P, F);

  // Save the input mesh DOFs
  Eigen::VectorXd xTarget(V.size());
  for(int i = 0; i < V.rows(); ++i)
    for(int j = 0; j < 3; ++j)
      xTarget(3 * i + j) = V(i, j);


  // Simulation
  std::cout << "**********\nSIMULATION\n**********\n";

  Timer optimTimer("Directions optimization");
  // Define simulation function
  auto func = simulationFunction(geometry, MrInv, theta1, theta2, E1, lambda1, lambda2, deltaLambda, thickness);

  // (Projected) Newton optimization
  newton(geometry, V, func, n_iter, lim, true, fixedIdx);
  optimTimer.stop();

  // Display distance with target mesh
  Eigen::MatrixXd VTarget = xTarget.reshaped<Eigen::RowMajor>(V.rows(), 3);
  Eigen::VectorXd d = (V - VTarget).cwiseProduct((V - VTarget)).rowwise().sum();
  d = d.array().sqrt();
  std::cout << "Avg distance = "
            << 100 * d.sum() / d.size() / (VTarget.colwise().maxCoeff() - VTarget.colwise().minCoeff()).norm()
            << "\n";
  std::cout << "Max distance = "
            << 100 * d.lpNorm<Eigen::Infinity>() /
                  (VTarget.colwise().maxCoeff() - VTarget.colwise().minCoeff()).norm()
            << "\n";

  // Directions optimizatiom
  std::cout << "***********************\nDIRECTIONS OPTIMIZATION\n***********************\n";

  optimTimer.start();
  // Define optimization function
  auto adjointFunc = adjointFunction(geometry, F, MrInv, theta1, E1, lambda1, lambda2, deltaLambda, thickness);
  // Optimize this energy function using SGN [Zehnder et al. 2021]
  sparse_gauss_newton(geometry, V, xTarget, MrInv, theta1, theta2, adjointFunc, fixedIdx, n_iter, lim, wM, wL, E1,
                      lambda1, lambda2, deltaLambda, thickness);
  optimTimer.stop();

  // Display distance with target mesh
  d = (V - VTarget).cwiseProduct((V - VTarget)).rowwise().sum();
  d = d.array().sqrt();
  std::cout << "Avg distance = "
            << 100 * d.sum() / d.size() / (VTarget.colwise().maxCoeff() - VTarget.colwise().minCoeff()).norm()
            << "\n";
  std::cout << "Max distance = "
            << 100 * d.lpNorm<Eigen::Infinity>() /
                  (VTarget.colwise().maxCoeff() - VTarget.colwise().minCoeff()).norm()
            << "\n";

  // Generate trajectories
  std::cout << "************\nTRAJECTORIES\n************\n";

  // Center flat mesh at 0
  Eigen::RowVector2d center = P.colwise().sum() / P.rows();
  for(int i = 0; i < P.rows(); ++i)
    P.row(i) = P.row(i) - center;

  Timer trajTimer("Trajectory optimization");

  double layerHeight = 0.08;
  double spacing = 0.4;
  int nLayers = 10;

  std::vector<Eigen::SparseMatrix<double>> subdivMat;
  subdivideMesh(geometry, V, P, F, subdivMat, spacing);

  ManifoldSurfaceMesh subdividedMesh(F);
  Eigen::MatrixXd P_3D(P.rows(), 3);
  P_3D.leftCols(2) = P;
  P_3D.col(2).setZero();
  VertexPositionGeometry geometryUV(subdividedMesh, P_3D);

  theta1 = computeStretchAngles(subdividedMesh, V, P, F);

  // convert face angles into vertex angles
  Eigen::VectorXd th1(V.rows());
  for(int i = 0; i < V.rows(); ++i)
  {
    double sumAngles = 0;
    int nFaces = 0;
    Vertex v = subdividedMesh.vertex(i);
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

  Eigen::VectorXd th2 = theta2.toVector();
  for(auto& mat: subdivMat)
    th2 = mat * th2;

  std::vector<std::vector<std::vector<Vector3>>> paths =
      generatePaths(geometryUV, th1, th2, layerHeight, nLayers, spacing, timeLimit);

  std::ofstream s(std::string(DATA_PATH_STR) + filename + ".path");
  for(int i = 0; i < nLayers; ++i)
  {
    writePaths(std::string(DATA_PATH_STR) + filename + ".path", paths[i], (i + 1) * layerHeight);
  }
}
