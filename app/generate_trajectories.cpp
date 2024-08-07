#include "generate_trajectories.h"
#include "parameterization.h"
#include "path_extraction.h"
#include "stretch_angles.h"


#include <polyscope/curve_network.h>
#include <polyscope/polyscope.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

void stripePattern(VertexPositionGeometry& geometry,
                   const Eigen::MatrixXd& _V,
                   const Eigen::MatrixXd& _P,
                   const Eigen::MatrixXi& _F,
                   VertexData<double>& vTheta2,
                   std::string filename,
                   double timeLimit,
                   double layerHeight,
                   double spacing,
                   int nLayers)
{
  Eigen::MatrixXd V = _V;
  Eigen::MatrixXd P = _P;
  Eigen::MatrixXi F = _F;
  std::vector<Eigen::SparseMatrix<double>> subdivMat;

  subdivideMesh(geometry, V, P, F, subdivMat, spacing);

  ManifoldSurfaceMesh subdividedMesh(F);
  Eigen::MatrixXd P_3D(P.rows(), 3);
  P_3D.leftCols(2) = P;
  P_3D.col(2).setZero();
  VertexPositionGeometry geometryUV(subdividedMesh, P_3D);

  FaceData<double> theta1 = computeStretchAngles(subdividedMesh, V, P, F);

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

  Eigen::VectorXd th2 = vTheta2.toVector();
  for(auto& mat: subdivMat)
    th2 = mat * th2;

  std::vector<std::vector<std::vector<Vector3>>> paths =
      generatePaths(geometryUV, th1, th2, layerHeight, nLayers, spacing, timeLimit);

  // geometryUV.requireDECOperators();
  // SparseMatrix<double> A = geometryUV.d0.transpose() * geometryUV.hodge1 * geometryUV.d0;
  // LDLTSolver solver;
  // solver.compute(A);

  // FaceData<double> theta2(subdividedMesh);
  // geometryUV.requireVertexIndices();
  // for(Face f: subdividedMesh.faces())
  // {
  //   double t = 0;
  //   for(Vertex v: f.adjacentVertices())
  //     t += th2(geometryUV.vertexIndices[v]);
  //   t /= 3;
  //   theta2[f] = t;
  // }

  // // #pragma omp parallel for schedule(static) num_threads(omp_get_max_threads() - 1)
  // for(int i = 0; i < nLayers; ++i)
  // {
  //   FaceData<Vector3> vectorField(subdividedMesh);
  //   for(Face f: subdividedMesh.faces())
  //   {
  //     double theta = theta1[f] + (i / (nLayers - 1.f) - 0.5) * theta2[f];
  //     vectorField[f] = {cos(theta), sin(theta), 0};
  //   }
  //   Eigen::VectorXd vCoordinate = curlFreeParameterization(geometryUV, vectorField, solver);
  //   paths[i] = orderPolylinesNew(paths[i], P, vCoordinate);
  // }

  for(int i = 0; i < nLayers; ++i)
  {
    paths[i] = orderPolylines(paths[i], timeLimit);
  }

  // Update height every layers
  static double height = 0;
  for(int i = 0; i < nLayers; ++i)
  {
    height += layerHeight + static_cast<float>(i) / (nLayers - 1) * 2 * (0.8 / nLayers - layerHeight);
    writePaths(filename + ".path", paths[i], height);
    drawPathsAndTravels(paths[i], spacing, i + 1);
  }
}


void drawPathsAndTravels(const std::vector<std::vector<Vector3>>& polylines, double spacing, int id)
{
  std::vector<Vector3> nodes;
  std::vector<std::array<size_t, 2>> edges;
  std::vector<double> colors;

  double c = 0;
  size_t k = 0;
  for(auto& polyline: polylines)
  {
    nodes.insert(nodes.end(), polyline.begin(), polyline.end());
    colors.insert(colors.end(), polyline.size() - 1, c);
    c += 1;

    for(size_t j = 0; j < polyline.size() - 1; ++j)
    {
      edges.push_back({k + j, k + j + 1});
    }
    k += polyline.size();
  }

  std::vector<Vector3> nodes2;
  std::vector<std::array<size_t, 2>> edges2;
  double distance = 0;
  for(size_t i = 0; i < polylines.size() - 1; ++i)
  {
    nodes2.push_back(polylines[i][polylines[i].size() - 1]);
    nodes2.push_back(polylines[i + 1][0]);
    edges2.push_back({2 * i, 2 * i + 1});
    distance += norm(polylines[i][polylines[i].size() - 1] - polylines[i + 1][0]);
  }
  std::cout << "Layer " << id << ": total travel distance = " << distance << "mm.\n";

  polyscope::registerCurveNetwork("Layer " + std::to_string(id), nodes, edges)->setRadius(spacing / 2, false);
  polyscope::getCurveNetwork("Layer " + std::to_string(id))
      ->addEdgeScalarQuantity("order", colors)
      ->setColorMap("blues")
      ->setEnabled(true);

  polyscope::registerCurveNetwork("Travel " + std::to_string(id), nodes2, edges2)
      ->setRadius(spacing / 2, false)
      ->setEnabled(false);
}
