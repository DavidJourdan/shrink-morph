#include "path_extraction.h"

#include "parameterization.h"
#include "stretch_angles.h"
#include "stripe_patterns.h"
#include "timer.h"

#include <geometrycentral/surface/direction_fields.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/ramer_douglas_peucker.h>
#include <polyscope/curve_network.h>
#include <polyscope/polyscope.h>

#ifdef USE_ORTOOLS
#include <ortools/constraint_solver/routing.h>
#include <ortools/constraint_solver/routing_enums.pb.h>
#include <ortools/constraint_solver/routing_index_manager.h>
#include <ortools/constraint_solver/routing_parameters.h>
#endif

#ifdef USE_PARDISO
#include <Eigen/PardisoSupport>
using LDLTSolver = Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>;
using LLTSolver = Eigen::PardisoLLT<Eigen::SparseMatrix<double>>;
using LUSolver = Eigen::PardisoLU<Eigen::SparseMatrix<double>>;
#else
using LDLTSolver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
using LLTSolver = Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>;
using LUSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>>;
#endif // USE_PARDISO

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace
{
#ifdef USE_ORTOOLS
std::vector<int64_t> computeReordering(const std::vector<Vector3>& endPoints, double timeLimit)
{
  using namespace operations_research;

  // Create the VRP routing model. The 1 means we are only looking for a single path.
  int nPaths = endPoints.size() / 4;
  RoutingIndexManager manager(2 * nPaths + 1, 1, RoutingIndexManager::NodeIndex{2 * nPaths});
  RoutingModel routing(manager);

  // For every path node, add a disjunction so that we do not also trace its reverse
  for(int64_t i = 0; i < nPaths; ++i)
    routing.AddDisjunction({2 * i, 2 * i + 1});

  // Wrap the distance function so that it converts to an integer, as or-tools requires.
  // Distances are rounded to the nearest nanometer value, should be precise enough
  const double COST_MULTIPLIER = 1e6;
  int64_t transitCallbackIndex = routing.RegisterTransitCallback([&](int64_t i, int64_t j) -> int64_t {
    int64_t fromNode = manager.IndexToNode(i).value();
    int64_t toNode = manager.IndexToNode(j).value();
    if(fromNode == 2 * nPaths || toNode == 2 * nPaths)
      return 0;
    double dist = (endPoints[2 * fromNode + 1] - endPoints[2 * toNode]).norm();
    return static_cast<int64_t>(COST_MULTIPLIER * dist);
  });
  routing.SetArcCostEvaluatorOfAllVehicles(transitCallbackIndex);

  RoutingSearchParameters searchParameters = DefaultRoutingSearchParameters();
  searchParameters.mutable_time_limit()->set_seconds(timeLimit);
  searchParameters.set_solution_limit(kint64max);
  searchParameters.set_first_solution_strategy(FirstSolutionStrategy::SAVINGS);
  searchParameters.set_local_search_metaheuristic(LocalSearchMetaheuristic::GUIDED_LOCAL_SEARCH);

  // solve routing problem
  auto finalAssignment = std::make_unique<Assignment>(routing.SolveWithParameters(searchParameters));

  // Iterate over the result
  std::vector<int64_t> result;
  result.reserve(nPaths);
  int64_t index = routing.Start(0);
  for(int i = 0; i < nPaths; ++i)
  {
    index = manager.IndexToNode(finalAssignment->Value(routing.NextVar(index))).value();
    assert(index <= 2 * nPaths);
    result.push_back(index);
  }

  // build polyline list from the list of endpoint indices
  return result;
}
#else
std::vector<int> greedyIndexing(const std::vector<Vector3>& endPoints, EmbeddedGeometryInterface& geometry)
{
  // find closest end-point to home
  Vector3 homePos = {-80, -80};
  double minDist = std::numeric_limits<double>::max();
  int minIdx;
  for(int j = 0; j < endPoints.size(); ++j)
  {
    if((endPoints[j] - homePos).norm() < minDist)
    {
      minDist = (endPoints[j] - homePos).norm();
      minIdx = j;
    }
  }

  // indices into endPoints paired with the direction of isoline (true = standard, false = reversed)
  std::vector<int> endPointIndices(endPoints.size() / 2);

  // tells whether isolines[i] has been visited or not
  std::vector<bool> visited(endPoints.size() / 2, false);
  visited[minIdx / 2] = true;

  // order end-point indices by proximity
  endPointIndices[0] = minIdx;
  Vector3 current;
  if(minIdx % 2 == 0)
    current = endPoints[minIdx + 1];
  else
    current = endPoints[minIdx - 1];

  for(int i = 1; i < endPoints.size() / 2; ++i)
  {
    // find closest point to current
    double minDist = std::numeric_limits<double>::max();
    int minIdx;
    for(int j = 0; j < endPoints.size(); ++j)
    {
      if(!visited[j / 2] && (current - endPoints[j]).norm() < minDist)
      {
        minDist = (current - endPoints[j]).norm();
        minIdx = j;
      }
    }
    // list index of isoline + whether it should be reversed or not
    endPointIndices[i] = minIdx;
    // set current to the other end-point
    if(minIdx % 2 == 0)
      current = endPoints[minIdx + 1];
    else
      current = endPoints[minIdx - 1];

    visited[minIdx / 2] = true;
  }

  return endPointIndices;
}
#endif
} // namespace

std::vector<std::vector<Vector3>>
orderPolylines(const std::vector<std::vector<Vector3>>& isolines, EmbeddedGeometryInterface& geometry, double timeLimit)
{
#ifdef USE_ORTOOLS
  // build list of isoline end-points, enumerated twice (one for each direction)
  std::vector<Vector3> endPoints(4 * isolines.size());

  for(size_t i = 0; i < isolines.size(); ++i)
  {
    endPoints[4 * i] = isolines[i][0];
    endPoints[4 * i + 1] = isolines[i][isolines[i].size() - 1];

    endPoints[4 * i + 2] = endPoints[4 * i + 1];
    endPoints[4 * i + 3] = endPoints[4 * i];
  }

  // TSP solving
  std::vector<int64_t> indices = computeReordering(endPoints, timeLimit);
#else
  std::vector<Vector3> endPoints(2 * isolines.size());

  for(int i = 0; i < isolines.size(); ++i)
  {
    endPoints[2 * i] = isolines[i][0];
    endPoints[2 * i + 1] = isolines[i][isolines[i].size() - 1];
  }

  std::vector<int> indices = greedyIndexing(endPoints, geometry);
#endif
  // assemble final polyline list
  std::vector<std::vector<Vector3>> polylines(isolines.size());
  for(size_t i = 0; i < isolines.size(); ++i)
  {
    auto isoline = isolines.at(indices[i] / 2);

    if(indices[i] % 2 == 0)
    {
      polylines[i] = isoline;
    }
    else // enumerate vertices from last to first
    {
      polylines[i].reserve(isoline.size());
      for(int j = isoline.size() - 1; j >= 0; --j)
      {
        polylines[i].push_back(isoline[j]);
      }
    }
  }

  return polylines;
}

std::vector<std::vector<Vector3>> simplifyPolylines(const std::vector<std::vector<Vector3>>& polylines, double z)
{
  std::vector<std::vector<Vector3>> simplified;

  for(auto& polyline: polylines)
  {
    // convert to Eigen
    Eigen::MatrixXd P(polyline.size(), 2);
    for(size_t i = 0; i < polyline.size(); ++i)
      P.row(i) << polyline[i].x, polyline[i].y;

    // simplify polyline (RDP)
    Eigen::MatrixXd S;
    Eigen::VectorXi J;
    igl::ramer_douglas_peucker(P, 1e-2, S, J);

    // convert back to geometrycentral
    std::vector<Vector3> simplifiedPolyline(S.rows());
    for(int i = 0; i < S.rows(); ++i)
      simplifiedPolyline[i] = Vector3{S(i, 0), S(i, 1), z};

    simplified.push_back(simplifiedPolyline);
  }
  return simplified;
}

void writePaths(const std::string& filename, const std::vector<std::vector<Vector3>>& paths, double height)
{
  std::ofstream s(filename, std::ios_base::app);
  for(const auto& path: paths)
  {
    s << path.size() << "\n";
    for(Vector3 p: path)
      s << p.x << " " << p.y << " " << height << "\n";
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

void stripePattern(VertexPositionGeometry& geometry,
                   const Eigen::MatrixXd& _V,
                   Eigen::MatrixXd& _P,
                   const Eigen::MatrixXi& _F,
                   VertexData<double>& vTheta2,
                   std::string filename,
                   double timeLimit,
                   double layerHeight,
                   double spacing,
                   int nLayers)
{
  Eigen::MatrixXd V = _V;
  Eigen::MatrixXi F = _F;
  std::vector<Eigen::SparseMatrix<double>> subdivMat;

  subdivideMesh(geometry, V, _P, F, subdivMat, spacing);

  ManifoldSurfaceMesh subdividedMesh(F);
  Eigen::MatrixXd P_3D(_P.rows(), 3);
  P_3D.leftCols(2) = _P;
  P_3D.col(2).setZero();
  VertexPositionGeometry geometryUV(subdividedMesh, P_3D);

  FaceData<double> theta1 = computeStretchAngles(subdividedMesh, V, _P, F);

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

  std::ofstream s(std::string(DATA_PATH_STR) + filename + ".path");
  for(int i = 0; i < nLayers; ++i)
  {
    writePaths(std::string(DATA_PATH_STR) + filename + ".path", paths[i], (i + 1) * layerHeight);
    drawPathsAndTravels(paths[i], spacing, i + 1);
  }
}

std::vector<std::vector<std::vector<geometrycentral::Vector3>>> generatePaths(EmbeddedGeometryInterface& geometry,
                                                                              const Eigen::VectorXd& theta1,
                                                                              const Eigen::VectorXd& theta2,
                                                                              double layerHeight,
                                                                              int nLayers,
                                                                              double spacing,
                                                                              double timeLimit)
{
  std::vector<std::vector<std::vector<Vector3>>> paths;

  LDLTSolver solver;

  SparseMatrix<double> massMatrix = computeRealVertexMassMatrix(geometry);
  SurfaceMesh& mesh = geometry.mesh;
  VertexData<double> frequencies(mesh, 1.0 / spacing);

  Vector<double> u = Vector<double>::Random(mesh.nVertices() * 2);
  for(int i = 0; i < nLayers / 2; ++i)
  {
    Timer timer("Layers " + std::to_string(2 * i) + " and " + std::to_string(2 * i + 1));
    geometry.requireVertexTangentBasis();

    // compute direction field
    VertexData<Vector2> directionField(mesh);
    for(size_t j = 0; j < mesh.nVertices(); ++j)
    {
      // interpolate orientations w.r.t layer
      directionField[j] = Vector2::fromAngle(2 * (theta1(j) + (i / (nLayers / 2 - 1.) - 1 / 2.) * theta2(j)));
      // express vectors in their respective vertex local bases (the stripes algorithm expects that)
      auto basisVector = geometry.vertexTangentBasis[j][0];
      Vector2 base{basisVector.x, basisVector.y};
      directionField[j] = -directionField[j] / base.pow(2);
    }

    // build matrices and factorize solver
    FaceData<int> fieldIndices = computeFaceIndex(geometry, directionField, 2);
    SparseMatrix<double> energyMatrix =
        buildVertexEnergyMatrix(geometry, directionField, fieldIndices, 2 * PI * frequencies);
    if(i == 0)
      solver.compute(energyMatrix);
    else
      solver.factorize(energyMatrix);

    // solve the eigenvalue problem
    Vector<double> x = u;
    double residual = eigenvectorResidual(energyMatrix, massMatrix, x);
    const double tol = 1e-5;
    while(residual > tol)
    {
      // Solve
      Vector<double> rhs = massMatrix * u;
      x = solver.solve(rhs);

      // Re-normalize
      double scale = std::sqrt(std::abs(x.dot(massMatrix * x)));
      x /= scale;

      // Update
      u = x;
      residual = eigenvectorResidual(energyMatrix, massMatrix, x);
    }

    // Copy the result to a VertexData vector
    VertexData<Vector2> parameterization(mesh);
    for(size_t j = 0; j < mesh.nVertices(); ++j)
    {
      parameterization[j].x = u(2 * j);
      parameterization[j].y = u(2 * j + 1);
    }

    CornerData<double> stripeValues;
    FaceData<int> stripeIndices;
    std::tie(stripeValues, stripeIndices) =
        computeTextureCoordinates(geometry, directionField, 2 * PI * frequencies, parameterization);

    // extract isolines
    for(int j = 0; j < 2; ++j)
    {
      const auto& [points, edges] = extractPolylinesFromStripePattern(geometry, stripeValues + j * PI, stripeIndices,
                                                                      fieldIndices, directionField);
      auto polylines = edgeToPolyline(points, edges);
      polylines = orderPolylines(polylines, geometry, timeLimit);
      paths.push_back(simplifyPolylines(polylines, (2 * i + j) * layerHeight));
    }
  }

  return paths;
}

std::vector<std::vector<Vector3>> edgeToPolyline(const std::vector<Vector3>& points,
                                                 const std::vector<std::array<size_t, 2>>& edges)
{
  // build up connectivity information (neighoring indices on polyline)
  std::vector<std::vector<size_t>> connectivity(points.size());

  for(auto& edge: edges)
  {
    connectivity[edge[0]].push_back(edge[1]);
    connectivity[edge[1]].push_back(edge[0]);
  }

  std::vector<std::vector<Vector3>> polylines;
  std::vector<bool> visited(points.size(), false);
  for(size_t i = 0; i < points.size(); ++i)
    if(!visited[i])
    {
      std::vector<Vector3> isoline;

      size_t nbOfPieces = connectivity[i].size();
      if(nbOfPieces == 0)
        continue; // ignore isolated points

      size_t currIdx = connectivity[i][0];
      size_t prevIdx = i;
      bool open = true;
      if(nbOfPieces == 2) // walk to the end of the line (if open)
      {
        while(connectivity[currIdx].size() == 2 && currIdx != i)
        {
          if(connectivity[currIdx][0] == prevIdx)
          {
            prevIdx = currIdx;
            currIdx = connectivity[currIdx][1];
          }
          else
          {
            prevIdx = currIdx;
            currIdx = connectivity[currIdx][0];
          }
        }
        open = currIdx != i;
        std::swap(prevIdx, currIdx); // revert the order to start from this end
      }

      isoline.push_back(points[prevIdx]);
      visited[prevIdx] = true;

      while((connectivity[currIdx].size() == 2 && open) || (currIdx != i && !open))
      {
        isoline.push_back(points[currIdx]);
        visited[currIdx] = true;

        if(connectivity[currIdx][0] == prevIdx)
        {
          prevIdx = currIdx;
          currIdx = connectivity[currIdx][1];
        }
        else
        {
          prevIdx = currIdx;
          currIdx = connectivity[currIdx][0];
        }
      }
      assert(connectivity[currIdx].size() == 1 || currIdx == i && !open);

      isoline.push_back(points[currIdx]);
      visited[currIdx] = true;

      polylines.push_back(isoline);
    }
  return polylines;
}
