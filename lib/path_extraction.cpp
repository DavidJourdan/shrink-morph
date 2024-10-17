#include "path_extraction.h"

#include "stripe_patterns.h"
#include "timer.h"

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/direction_fields.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/ramer_douglas_peucker.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace
{
std::vector<int> greedyIndexing(const std::vector<Vector3>& endPoints)
{
  // find closest end-point to home
  Vector3 homePos = {-150, -150};
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
} // namespace

std::vector<std::vector<Vector3>> orderPolylines(const std::vector<std::vector<Vector3>>& isolines)
{
  std::vector<Vector3> endPoints(2 * isolines.size());

  for(int i = 0; i < isolines.size(); ++i)
  {
    endPoints[2 * i] = isolines[i][0];
    endPoints[2 * i + 1] = isolines[i][isolines[i].size() - 1];
  }

  std::vector<int> indices = greedyIndexing(endPoints);
  // assemble final polyline list
  std::vector<std::vector<Vector3>> polylines(isolines.size());
  for(size_t i = 0; i < isolines.size(); ++i)
  {
    auto isoline = isolines[indices[i] / 2];

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
  // std::vector<std::vector<Vector3>> polylines;
  // for(size_t i = 0; i < isolines.size(); ++i)
  // {
  //   int n = isolines[indices[i] / 2].size();
  //   auto previous = [&](int i) {
  //     if(indices[i - 1] % 2 == 0)
  //       return indices[i - 1] + 1;
  //     else
  //       return indices[i - 1] - 1;
  //   };
  //   if(i > 0 && norm(endPoints[previous(i)] - endPoints[indices[i]]) < 2)
  //   { // append to previous isoline
  //     if(indices[i] % 2 == 0)
  //       for(int j = 0; j < n; ++j)
  //         polylines[polylines.size() - 1].push_back(isolines[indices[i] / 2][j]);
  //     else
  //       for(int j = n - 1; j >= 0; --j)
  //         polylines[polylines.size() - 1].push_back(isolines[indices[i] / 2][j]);
  //   }
  //   else
  //   { // create new isoline
  //     if(indices[i] % 2 == 0)
  //     {
  //       polylines.push_back(isolines[indices[i] / 2]);
  //     }
  //     else // enumerate vertices from last to first
  //     {
  //       std::vector<Vector3> isoline;
  //       isoline.reserve(n);
  //       for(int j = n - 1; j >= 0; --j)
  //       {
  //         isoline.push_back(isolines[indices[i] / 2][j]);
  //       }
  //       polylines.push_back(isoline);
  //     }
  //   }
  // }
  // std::cout << isolines.size() << " " << polylines.size() << "\n";

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

std::vector<std::vector<geometrycentral::Vector3>> generateOneLayer(EmbeddedGeometryInterface& geometry,
                                                                    const Eigen::VectorXd& theta1,
                                                                    const Eigen::VectorXd& theta2,
                                                                    const Eigen::SparseMatrix<double> &massMatrix,
                                                                    Eigen::VectorXd& u,
                                                                    LDLTSolver &solver,
                                                                    int i,
                                                                    int nLayers,
                                                                    double layerHeight,
                                                                    double spacing)
{
  geometry.requireVertexTangentBasis();
  SurfaceMesh& mesh = geometry.mesh;
  VertexData<double> frequencies(mesh, 1.0 / spacing);

  // compute direction field
  VertexData<Vector2> directionField(mesh);
  for(size_t j = 0; j < mesh.nVertices(); ++j)
  {
    // interpolate orientations w.r.t layer
    directionField[j] = Vector2::fromAngle(2 * (theta1(j) + (i / (nLayers - 1.) - 1 / 2.) * theta2(j)));
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
  Vector<double> x;
  double residual = eigenvectorResidual(energyMatrix, massMatrix, u);
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
  const auto& [points, edges] =
      extractPolylinesFromStripePattern(geometry, stripeValues, stripeIndices, fieldIndices, directionField);
  auto polylines = edgeToPolyline(points, edges);
  polylines = orderPolylines(polylines);
  return simplifyPolylines(polylines, (i + 1) * layerHeight);
}

std::vector<std::vector<std::vector<geometrycentral::Vector3>>> generatePaths(EmbeddedGeometryInterface& geometry,
                                                                              const Eigen::VectorXd& theta1,
                                                                              const Eigen::VectorXd& theta2,
                                                                              double layerHeight,
                                                                              int nLayers,
                                                                              double spacing)
{
  std::vector<std::vector<std::vector<Vector3>>> paths;
  SparseMatrix<double> massMatrix = computeRealVertexMassMatrix(geometry);
  LDLTSolver solver;
  Vector<double> u = Vector<double>::Random(geometry.mesh.nVertices() * 2);
  for(int i = 0; i < nLayers; ++i)
  {
    Timer timer("Layer " + std::to_string(i));
    paths.push_back(generateOneLayer(geometry, theta1, theta2, massMatrix, u, solver, i, nLayers, layerHeight, spacing));
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

std::tuple<VertexData<double>, FaceData<double>, EdgeData<double>> hodgeDecomposition(VertexPositionGeometry& geometry,
                                                                                      const EdgeData<double>& oneForm)
{
  SurfaceMesh& mesh = geometry.mesh;
  geometry.requireDECOperators();

  Vector<double> omega = oneForm.toVector();

  // Solve for the curl-free part alpha
  SparseMatrix<double> A = geometry.d0.transpose() * geometry.hodge1 * geometry.d0;

  Vector<double> rhs = geometry.d0.transpose() * geometry.hodge1 * omega;
  Vector<double> alpha = solvePositiveDefinite(A, rhs);

  // Solve for the divergence free part beta
  SparseMatrix<double> B = geometry.d1 * geometry.hodge1Inverse * geometry.d1.transpose();

  rhs = geometry.d1 * omega;
  Vector<double> beta = solveSquare(B, rhs);

  beta = geometry.hodge2 * beta;

  EdgeData<double> gamma(mesh);
  gamma.fromVector(omega - geometry.d0 * alpha -
                   geometry.hodge1Inverse * geometry.d1.transpose() * geometry.hodge2 * beta);

  VertexData<double> vAlpha(mesh);
  vAlpha.fromVector(alpha);

  FaceData<double> fBeta(mesh);
  fBeta.fromVector(beta);

  return std::make_tuple(vAlpha, fBeta, gamma);
}

Eigen::VectorXd
curlFreeParameterization(const Eigen::MatrixXd& P, const Eigen::MatrixXi& F, const Eigen::VectorXd& theta)
{
  // create geometry-central objects
  ManifoldSurfaceMesh mesh(F);
  Eigen::MatrixXd P_3D(P.rows(), 3);
  P_3D.leftCols(2) = P;
  P_3D.col(2).setZero();
  VertexPositionGeometry geometry(mesh, P_3D);

  FaceData<Vector3> vectorField(mesh);
  for(int i = 0; i < F.rows(); ++i)
  {
    vectorField[i] = {cos(theta(i)), sin(theta(i)), 0};
  }

  return curlFreeParameterization(geometry, vectorField);
}

Eigen::VectorXd curlFreeParameterization(VertexPositionGeometry& geometry, const FaceData<Vector3>& vectorField)
{
  geometry.requireDECOperators();

  SparseMatrix<double> A = geometry.d0.transpose() * geometry.hodge1 * geometry.d0;
  LDLTSolver solver;
  solver.compute(A);

  return curlFreeParameterization(geometry, vectorField, solver);
}

Eigen::VectorXd
curlFreeParameterization(VertexPositionGeometry& geometry, const FaceData<Vector3>& vectorField, LDLTSolver& solver)
{
  SurfaceMesh& mesh = geometry.mesh;
  geometry.requireDECOperators();

  // Construct two one forms: one from the vector field, and the other from its 90 degree rotation
  EdgeData<double> oneFormPerp(mesh);
  geometry.requireVertexNormals();
  for(Edge e: mesh.edges())
  {
    Vector3 v = geometry.vertexPositions[e.secondVertex()] - geometry.vertexPositions[e.firstVertex()];
    Vector3 w;
    if(e.isBoundary())
    {
      w = normalize(vectorField[e.halfedge().face()]);
    }
    else
    {
      w = normalize(vectorField[e.halfedge().face()] + vectorField[e.halfedge().twin().face()]);
    }

    w = cross(w, geometry.vertexNormals[e.firstVertex()]);
    oneFormPerp[e] = dot(v, w);
  }

  // Solve for the curl-free part (Helmholtz-Hodge decomposition)
  Eigen::VectorXd rhs = geometry.d0.transpose() * geometry.hodge1 * oneFormPerp.toVector();

  return solver.solve(rhs);
}

std::vector<std::vector<Vector3>> orderPolylinesNew(const std::vector<std::vector<Vector3>>& isolines,
                                                    const Eigen::MatrixX2d& P,
                                                    const Eigen::MatrixXi& F,
                                                    const Eigen::VectorXd& theta)
{
  Eigen::VectorXd vCoordinate = curlFreeParameterization(P, F, theta);
  return orderPolylinesNew(isolines, P, vCoordinate);
}

std::vector<std::vector<Vector3>> orderPolylinesNew(const std::vector<std::vector<Vector3>>& isolines,
                                                    const Eigen::MatrixX2d& P,
                                                    const Eigen::VectorXd& vCoordinate)
{
  Timer timer("Order polylines");

  // // build grid
  // Eigen::Vector2d dimensions = P.colwise().maxCoeff() - P.colwise().minCoeff();
  // Eigen::MatrixXd values;
  // values.setZero(round(dimensions(0)) + 1, round(dimensions(1)) + 1);

  // for(int i = 0; i < P.rows(); ++i)
  // {
  //   int x = round(P(i, 0) + dimensions(0) / 2);
  //   int y = round(P(i, 1) + dimensions(1) / 2);
  //   values(x, y) = vCoordinate(i, 1);
  // }

  std::vector<std::pair<double, int>> toSort(isolines.size());

  // #pragma omp parallel for schedule(static) num_threads(omp_get_max_threads() - 1)
  for(int k = 0; k < isolines.size(); ++k)
  {
    const auto& isoline = isolines[k];
    double value = 0;
    for(const Vector3& v: isoline)
    {
      // find closest vertex in mesh
      Eigen::RowVector2d p(v.x, v.y);
      double minDist = (p - P.row(0)).norm();
      int minIdx = 0;
      for(int i = 0; i < P.rows(); ++i)
      {
        double dist = (p - P.row(i)).norm();
        if(dist < minDist)
        {
          minDist = dist;
          minIdx = i;
        }
      }
      value += vCoordinate(minIdx, 1);
      // int x = round(v.x) + dimensions(0) / 2;
      // int y = round(v.y) + dimensions(1) / 2;
      // value += values(x, y);
    }
    value /= isoline.size();
    toSort[k].first = value;
    toSort[k].second = k;
  }

  std::sort(toSort.begin(), toSort.end(), [](auto& a, auto& b) { return a.first < b.first; });
  // for(const auto& [value, k]: toSort)
  //   std::cout << value << ", " << k << " " << "\n";
  std::vector<std::vector<Vector3>> sortedIsolines(isolines.size());
  Vector3 prevEndpoint = isolines[toSort[0].second][0];

  for(int i = 0; i < isolines.size(); ++i)
  {
    double dist1 = norm(prevEndpoint - isolines[toSort[i].second][0]);
    double dist2 = norm(prevEndpoint - isolines[toSort[i].second].back());
    if(dist1 < dist2)
    {
      sortedIsolines[i] = isolines[toSort[i].second];
    }
    else // enumerate vertices from last to first
    {
      sortedIsolines[i].reserve(isolines[toSort[i].second].size());
      for(int j = isolines[toSort[i].second].size() - 1; j >= 0; --j)
        sortedIsolines[i].push_back(isolines[toSort[i].second][j]);
    }
    prevEndpoint = sortedIsolines[i].back();
  }
  return sortedIsolines;
}