
#include "parameterization.h"

#include <TinyAD/Utils/Helpers.hh>
#include <igl/bfs_orient.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/loop.h>
#include <igl/map_vertices_to_circle.h>

Eigen::MatrixXd tutte_embedding(const Eigen::MatrixXd &_V, const Eigen::MatrixXi &_F)
{
  Eigen::VectorXi b;  // #constr boundary constraint indices
  Eigen::MatrixXd bc; // #constr-by-2 2D boundary constraint positions
  Eigen::MatrixXd P;  // #V-by-2 2D vertex positions
  Eigen::MatrixXi F = _F;

  // Identify boundary loops
  std::vector<std::vector<int>> loops;
  igl::boundary_loop(F, loops);

  // find longest boundary loop
  size_t idxMax = -1;
  double maxLen = 0;
  for(size_t i = 0; i < loops.size(); ++i)
  {
    double len = 0;
    for(size_t j = 0; j < loops[i].size(); ++j)
      len += (_V.row(loops[i][(j + 1) % loops[i].size()]) - _V.row(loops[i][j])).norm();

    if(len > maxLen)
    {
      maxLen = len;
      idxMax = i;
    }
  }

  // Set boundary vertex positions
  b.resize(loops[idxMax].size(), 1);
  for(size_t i = 0; i < loops[idxMax].size(); ++i)
    b(i) = loops[idxMax][i];
  igl::map_vertices_to_circle(_V, b, bc);

  // Fill-in other holes
  for(size_t i = 0; i < loops.size(); ++i)
    if(i != idxMax)
    {
      int nB = loops[i].size();
      int nF = F.rows();
      F.conservativeResize(nF + nB - 2, 3);
      for(int j = 0; j < nB - 2; ++j)
        F.row(nF + j) << loops[i][0], loops[i][j + 1], loops[i][j + 2];
    }

  // make sure normals are consistent
  Eigen::VectorXi C;
  igl::bfs_orient(F, F, C);

  // Compute interior vertex positions
  igl::harmonic(F, b, bc, 1, P);

  return P;
}

Eigen::MatrixXd localGlobal(const Eigen::MatrixXd &V,
                            const Eigen::MatrixXi &F,
                            double lambda1,
                            double lambda2,
                            LocalGlobalSolver &solver)
{
  using namespace Eigen;

  MatrixXd P = tutte_embedding(V, F);

  // run ARAP
  solver.solve(P, 1., 1.);

  // find vertex closest to the X axis in the input mesh
  RowVector3d center3D = V.colwise().sum() / V.rows();
  int maxIdx;
  double maxX = 0;
  for(int i = 0; i < V.rows(); ++i)
  {
    Vector3d dir = (V.row(i) - center3D).normalized();
    if(dir(0) > maxX)
    {
      maxX = dir(0);
      maxIdx = i;
    }
  }

  // center flat mesh
  RowVector2d center2D = P.colwise().sum() / P.rows();
  P.rowwise() -= center2D;

  // rotate flat mesh (to make the orientation consistent with the input mesh)
  Matrix2d R;
  R.col(0) = P.row(maxIdx).normalized();
  R(1, 1) = R(0, 0);
  R(0, 1) = -R(1, 0);
  P = P * R;

  // change aspect ratio
  P.col(0) /= lambda1;
  P.col(1) /= lambda2;

  // run Local-Global algorithm
  solver.solve(P, 1 / lambda2, 1 / lambda1);

  // Repeat center and rotate operations
  center2D = P.colwise().sum() / P.rows();
  P.rowwise() -= center2D;

  R.col(0) = P.row(maxIdx).normalized();
  R(1, 1) = R(0, 0);
  R(0, 1) = -R(1, 0);
  P = P * R;

  return P;
}

geometrycentral::surface::FaceData<Eigen::Matrix2d>
precomputeParamData(geometrycentral::surface::VertexPositionGeometry &geometry)
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  SurfaceMesh &mesh = geometry.mesh;

  FaceData<Eigen::Matrix2d> rest_shapes(mesh);
  geometry.requireVertexPositions();

  auto toEigen = [](const geometrycentral::Vector3 &_v) { return Eigen::Vector3d(_v.x, _v.y, _v.z); };

  for(Face f: mesh.faces())
  {
    // Get 3D vertex positions
    Eigen::Vector3d ar_3d = toEigen(geometry.vertexPositions[f.halfedge().vertex()]);
    Eigen::Vector3d br_3d = toEigen(geometry.vertexPositions[f.halfedge().next().vertex()]);
    Eigen::Vector3d cr_3d = toEigen(geometry.vertexPositions[f.halfedge().next().next().vertex()]);

    // Set up local 2D coordinate system
    Eigen::Vector3d n = (br_3d - ar_3d).cross(cr_3d - ar_3d);
    Eigen::Vector3d b1 = (br_3d - ar_3d).normalized();
    Eigen::Vector3d b2 = n.cross(b1).normalized();

    // Express a, b, c in local 2D coordiante system
    Eigen::Vector2d ar_2d(0.0, 0.0);
    Eigen::Vector2d br_2d((br_3d - ar_3d).dot(b1), 0.0);
    Eigen::Vector2d cr_2d((cr_3d - ar_3d).dot(b1), (cr_3d - ar_3d).dot(b2));

    // Save 2-by-2 matrix with edge vectors as colums
    rest_shapes[f] = TinyAD::col_mat(br_2d - ar_2d, cr_2d - ar_2d);
  };
  return rest_shapes;
}

geometrycentral::surface::FaceData<Eigen::Matrix2d>
precomputeSimData(geometrycentral::surface::ManifoldSurfaceMesh &mesh,
                  const Eigen::MatrixXd &P,
                  const Eigen::MatrixXi &F)
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  FaceData<Eigen::Matrix2d> rest_shapes(mesh);

  for(int i = 0; i < F.rows(); ++i)
  {
    // Get 2D vertex positions
    Eigen::Vector2d a = P.row(F(i, 0));
    Eigen::Vector2d b = P.row(F(i, 1));
    Eigen::Vector2d c = P.row(F(i, 2));
    Eigen::Matrix2d Mr = TinyAD::col_mat(b - a, c - a);

    // Save 2-by-2 matrix with edge vectors as colums
    rest_shapes[i] = Mr.inverse();
  }
  return rest_shapes;
}

geometrycentral::surface::EdgeData<double>
computeDualCotanWeights(geometrycentral::surface::IntrinsicGeometryInterface &geometry)
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  SurfaceMesh &mesh = geometry.mesh;

  geometry.requireEdgeLengths();
  geometry.requireHalfedgeCotanWeights();
  geometry.requireFaceAreas();

  EdgeData<double> dualCotanWeights(mesh);
  for(Edge e: mesh.edges())
  {
    if(e.isBoundary())
    {
      dualCotanWeights[e] = 0;
      continue;
    }

    Halfedge he = e.halfedge();

    // cot(x/2) = cotx + 1 / sinx
    double sinaInv = 0.5 * geometry.edgeLengths[he.edge()] * geometry.edgeLengths[he.next().next().edge()] /
                     geometry.faceAreas[he.face()];
    double cota = geometry.halfedgeCotanWeights[he.next()] + sinaInv;

    double sinbInv = 0.5 * geometry.edgeLengths[he.twin().edge()] * geometry.edgeLengths[he.twin().next().edge()] /
                     geometry.faceAreas[he.twin().face()];
    double cotb = geometry.halfedgeCotanWeights[he.twin().next().next()] + sinbInv;

    double sincInv =
        0.5 * geometry.edgeLengths[he.edge()] * geometry.edgeLengths[he.next().edge()] / geometry.faceAreas[he.face()];
    double cotc = geometry.halfedgeCotanWeights[he.next().next()] + sincInv;

    double sindInv = 0.5 * geometry.edgeLengths[he.twin().edge()] *
                     geometry.edgeLengths[he.twin().next().next().edge()] / geometry.faceAreas[he.twin().face()];
    double cotd = geometry.halfedgeCotanWeights[he.twin().next()] + sindInv;

    // cot((a+b)/2) = (cot(a/2) cot(b/2) - 1) / (cot(a/2) + cot(b/2))
    double cotab = (cota * cotb - 1) / (cota + cotb);
    double cotcd = (cotc * cotd - 1) / (cotc + cotd);

    dualCotanWeights[e] = (cotab + cotcd) / 2;
  }
  return dualCotanWeights;
}

void subdivideMesh(geometrycentral::surface::VertexPositionGeometry &geometry,
                   Eigen::MatrixXd &V,
                   Eigen::MatrixXd &P,
                   Eigen::MatrixXi &F,
                   std::vector<Eigen::SparseMatrix<double>> &subdivMat,
                   double threshold)
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  SurfaceMesh &mesh = geometry.mesh;

  double avgEdgeLength = 0;
  geometry.requireVertexPositions();
  for(Edge e: mesh.edges())
    //     avgEdgeLength += norm(geometry.vertexPositions[e.firstVertex()] -
    //     geometry.vertexPositions[e.secondVertex()]);
    // avgEdgeLength /= mesh.nEdges();
    if(norm(geometry.vertexPositions[e.firstVertex()] - geometry.vertexPositions[e.secondVertex()]) > avgEdgeLength)
      avgEdgeLength = norm(geometry.vertexPositions[e.firstVertex()] - geometry.vertexPositions[e.secondVertex()]);

  while(avgEdgeLength > threshold)
  {
    Eigen::MatrixX3i tempF = F;
    Eigen::SparseMatrix<double> S;
    igl::loop(V.rows(), tempF, S, F);

    V = S * V;
    P = S * P;

    subdivMat.push_back(S);

    avgEdgeLength /= 2;
  }
}

TinyAD::ScalarFunction<2, double, geometrycentral::surface::Vertex>
parameterizationFunction(geometrycentral::surface::VertexPositionGeometry &geometry,
                         double wPhi,
                         double lambda1,
                         double lambda2)
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  SurfaceMesh &mesh = geometry.mesh;

  // precompute inverse jacobians
  FaceData<Eigen::Matrix2d> restShapes = precomputeParamData(geometry);

  // Set up function with 2D vertex positions as variables.
  TinyAD::ScalarFunction<2, double, Vertex> func = TinyAD::scalar_function<2>(mesh.vertices());

  // Main objective term
  func.add_elements<3>(mesh.faces(), [&, lambda1, lambda2, restShapes](auto &element) -> TINYAD_SCALAR_TYPE(element) {
    // Evaluate element using either double or TinyAD::Double
    using T = TINYAD_SCALAR_TYPE(element);

    // Get variable 2D vertex positions
    Face f = element.handle;
    Eigen::Vector2<T> a = element.variables(f.halfedge().vertex());
    Eigen::Vector2<T> b = element.variables(f.halfedge().next().vertex());
    Eigen::Vector2<T> c = element.variables(f.halfedge().next().next().vertex());

    Eigen::Matrix2<T> M = TinyAD::col_mat(b - a, c - a);

    // Get constant 2D rest shape of f
    Eigen::Matrix2d Mr = restShapes[f];
    double A = 0.5 * Mr.determinant();

    // Compute symmetric Dirichlet energy
    Eigen::Matrix2<T> J = M * Mr.inverse();

    T x = hypot(J(0, 0) + J(1, 1), J(0, 1) - J(1, 0)) / 2;
    T y = hypot(J(0, 1) + J(1, 0), J(0, 0) - J(1, 1)) / 2;

    T s0 = x + y;
    T s1 = x - y;

    return 0.5 * A * (pow(s1 - 1 / lambda2, 2) + pow(s0 - 1 / lambda1, 2));
  });

  geometry.requireFaceIndices();
  // compute dual halfedge cotan weigths
  EdgeData<double> dualCotanWeights = computeDualCotanWeights(geometry);
  // Smoothness regularization (Dirichlet energy)
  func.add_elements<4>(
      mesh.edges(),
      [&, lambda1, lambda2, wPhi, dualCotanWeights, restShapes](auto &element) -> TINYAD_SCALAR_TYPE(element) {
        // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);

        auto edge = element.handle;

        if(edge.isBoundary())
          return (T)0;

        std::vector<T> angles;

        for(Face f: edge.adjacentFaces())
        {
          Eigen::Vector2<T> a = element.variables(f.halfedge().vertex());
          Eigen::Vector2<T> b = element.variables(f.halfedge().next().vertex());
          Eigen::Vector2<T> c = element.variables(f.halfedge().next().next().vertex());

          Eigen::Matrix2<T> M = TinyAD::col_mat(b - a, c - a);

          // Get constant 2D rest shape of f
          Eigen::Matrix2d Mr = restShapes[f];

          // Compute symmetric Dirichlet energy
          Eigen::Matrix2<T> J = M * Mr.inverse();

          T a1 = atan2(J(0, 1) + J(1, 0), J(0, 0) - J(1, 1));
          T a2 = atan2(J(1, 0) - J(0, 1), J(0, 0) + J(1, 1));
          angles.push_back((a1 + a2) / 2);
        }

        return dualCotanWeights[edge] * wPhi * pow(sin(angles[1] - angles[0]), 2);
      });

  return func;
}
