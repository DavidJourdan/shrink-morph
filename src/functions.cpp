#include "functions.h"

#include "constants.h"
#include "simulation_utils.h"

#include <TinyAD/Utils/Helpers.hh>

using namespace geometrycentral::surface;

TinyAD::ScalarFunction<3, double, Vertex> simulationFunction(IntrinsicGeometryInterface& geometry,
                                                             const FaceData<Eigen::Matrix2d>& MrInv,
                                                             FaceData<double>& theta1,
                                                             VertexData<double>& theta2,
                                                             double E1,
                                                             double lambda1,
                                                             double lambda2,
                                                             double deltaLambda,
                                                             double thickness,
                                                             double E2)
{
  SurfaceMesh& mesh = geometry.mesh;

  // Set up function with 3D vertex positions as variables.
  TinyAD::ScalarFunction<3, double, Vertex> func = TinyAD::scalar_function<3>(mesh.vertices());

  // 1st fundamental form
  func.add_elements<3>(mesh.faces(), [&, E1, E2, lambda1, lambda2](auto& element) -> TINYAD_SCALAR_TYPE(element) {
    // Evaluate element using either double or TinyAD::Double
    using T = TINYAD_SCALAR_TYPE(element);

    // Get variable 2D vertex positions
    Face f = element.handle;
    Eigen::Vector3<T> a = element.variables(f.halfedge().vertex());
    Eigen::Vector3<T> b = element.variables(f.halfedge().next().vertex());
    Eigen::Vector3<T> c = element.variables(f.halfedge().next().next().vertex());

    double dA = 0.5 / MrInv[f].determinant();

    Eigen::Matrix<T, 3, 2> M = TinyAD::col_mat(b - a, c - a);

    // Compute deformation gradient
    Eigen::Matrix<T, 3, 2> F = M * (MrInv[f] * Eigen::Rotation2D<double>(theta1[f]));

    Eigen::Matrix2d S;
    S << 1 / lambda1 / lambda1, 0, 0, 1 / lambda2 / lambda2;
    Eigen::Matrix2<T> E = S * F.transpose() * F - Eigen::Matrix2d::Identity();

    Eigen::Vector3<T> EVoigt(E(0, 0), E(1, 1), E(0, 1) + E(1, 0));
    // clang-format off
    Eigen::Matrix3d C;
    C <<  E1, 0.4 * std::sqrt(E1 * E2), 0,
          0.4 * std::sqrt(E1 * E2), E2, 0,
          0, 0, 0.6 / 2 * std::sqrt(E1 * E2);
    // clang-format on
    return EVoigt.dot(C * EVoigt) * dA;
  });

  // 2nd fundamental form
  func.add_elements<6>(mesh.faces(), [&, E1, E2, lambda1, lambda2, deltaLambda](auto& element) -> TINYAD_SCALAR_TYPE(element) {
    // Evaluate element using either double or TinyAD::Double
    using T = TINYAD_SCALAR_TYPE(element);

    // Get variable 2D vertex positions
    Face f = element.handle;
    Eigen::Vector3<T> a = element.variables(f.halfedge().vertex());
    Eigen::Vector3<T> b = element.variables(f.halfedge().next().vertex());
    Eigen::Vector3<T> c = element.variables(f.halfedge().next().next().vertex());

    Eigen::Matrix<T, 3, 2> M = TinyAD::col_mat(b - a, c - a);

    double dA = 0.5 / MrInv[f].determinant();
    double deltaTheta = 0;
    for(Vertex v: f.adjacentVertices())
      deltaTheta += theta2[v] / 3;

    // Compute deformation gradient
    Eigen::Matrix<T, 3, 2> F = M * (MrInv[f] * Eigen::Rotation2D<double>(theta1[f]));

    Eigen::Matrix2d J;
    J << 0, 1, 1, 0;
    J *= 0.5 * deltaTheta * (lambda2 * lambda2 - lambda1 * lambda1);

    J(0, 0) += -deltaLambda * lambda1;

    Eigen::Matrix2d RJRt = Eigen::Rotation2D<double>(theta1[f]) * J * Eigen::Rotation2D<double>(-theta1[f]) / thickness;

    // Compute normal
    Eigen::Vector3<T> n = M.col(0).cross(M.col(1));

    Eigen::Matrix3<T> L = Eigen::Matrix3<T>::Zero();
    for(Halfedge he: f.adjacentHalfedges())
    {
      if(he.edge().isBoundary())
        continue;
      Eigen::Vector3<T> e = element.variables(he.next().vertex()) - element.variables(he.vertex());

      // compute dihedral angle theta
      Eigen::Vector3<T> nf =
          (element.variables(he.twin().next().next().vertex()) - element.variables(he.vertex())).cross(e);
      T theta = atan2(n.cross(nf).dot(e), e.norm() * nf.dot(n));

      Eigen::Vector3<T> t = n.cross(e);

      // add edge contribution
      L += theta * t.normalized() * t.transpose();
    }
    L /= n.squaredNorm();

    Eigen::Matrix2d S;
    S << 1 / lambda1 / lambda1, 0, 0, 1 / lambda2 / lambda2;
    Eigen::Matrix2<T> E = S * (F.transpose() * L * F - RJRt);

    Eigen::Vector3<T> EVoigt(E(0, 0), E(1, 1), E(0, 1) + E(1, 0));
    // clang-format off
    Eigen::Matrix3d C;
    C <<  E1, 0.4 * std::sqrt(E1 * E2), 0,
          0.4 * std::sqrt(E1 * E2), E2, 0,
          0, 0, 0.6 / 2 * std::sqrt(E1 * E2);
    // clang-format on
    return std::pow(thickness, 2) * EVoigt.dot(C * EVoigt) * dA / 3;
  });

  return func;
}

TinyAD::ScalarFunction<1, double, Eigen::Index> adjointFunction(IntrinsicGeometryInterface& geometry,
                                                                const Eigen::MatrixXi& F,
                                                                const FaceData<Eigen::Matrix2d>& MrInv,
                                                                const FaceData<double>& theta1,
                                                                double E1,
                                                                double lambda1,
                                                                double lambda2,
                                                                double deltaLambda,
                                                                double thickness,
                                                                double E2)
{
  SurfaceMesh& mesh = geometry.mesh;

  // Set up function with 3D vertex positions as variables.
  TinyAD::ScalarFunction<1, double, Eigen::Index> func =
      TinyAD::scalar_function<1>(TinyAD::range(3 * mesh.nVertices() + mesh.nVertices()));

  // 1st fundamental form
  func.add_elements<9>(
      TinyAD::range(F.rows()), [&, E1, E2, lambda1, lambda2](auto& element) -> TINYAD_SCALAR_TYPE(element) {
    // Evaluate element using either double or TinyAD::Double
    using T = TINYAD_SCALAR_TYPE(element);
    Eigen::Index f_idx = element.handle;

    // Get 3D vertex positions
    Eigen::Matrix<T, 3, 2> M;
    M << element.variables(3 * F(f_idx, 1) + 0) - element.variables(3 * F(f_idx, 0) + 0),
        element.variables(3 * F(f_idx, 2) + 0) - element.variables(3 * F(f_idx, 0) + 0),
        element.variables(3 * F(f_idx, 1) + 1) - element.variables(3 * F(f_idx, 0) + 1),
        element.variables(3 * F(f_idx, 2) + 1) - element.variables(3 * F(f_idx, 0) + 1),
        element.variables(3 * F(f_idx, 1) + 2) - element.variables(3 * F(f_idx, 0) + 2),
        element.variables(3 * F(f_idx, 2) + 2) - element.variables(3 * F(f_idx, 0) + 2);

    double dA = 0.5 / MrInv[f_idx].determinant();

    // Compute deformation gradient
    Eigen::Matrix<T, 3, 2> F = M * (MrInv[f_idx] * Eigen::Rotation2D<double>(theta1[f_idx]));

    Eigen::Matrix2d S;
    S << 1 / lambda1 / lambda1, 0, 0, 1 / lambda2 / lambda2;
    Eigen::Matrix2<T> E = S * F.transpose() * F - Eigen::Matrix2d::Identity();

    Eigen::Vector3<T> EVoigt(E(0, 0), E(1, 1), E(0, 1) + E(1, 0));
    // clang-format off
    Eigen::Matrix3d C;
    C <<  E1, 0.4 * std::sqrt(E1 * E2), 0,
          0.4 * std::sqrt(E1 * E2), E2, 0,
          0, 0, 0.6 / 2 * std::sqrt(E1 * E2);
    // clang-format on
    return EVoigt.dot(C * EVoigt) * dA;
  });

  // 2nd fundamental form
  geometry.requireVertexIndices();
  func.add_elements<
      3 * 6 + 3>(TinyAD::range(F.rows()), [&, E1, E2, lambda1, lambda2, deltaLambda](auto& element) -> TINYAD_SCALAR_TYPE(element) {
    // Evaluate element using either double or TinyAD::Double
    using T = TINYAD_SCALAR_TYPE(element);
    Eigen::Index f_idx = element.handle;

    // Get 3D vertex positions
    Eigen::Matrix<T, 3, 2> M;
    M << element.variables(3 * F(f_idx, 1) + 0) - element.variables(3 * F(f_idx, 0) + 0),
        element.variables(3 * F(f_idx, 2) + 0) - element.variables(3 * F(f_idx, 0) + 0),
        element.variables(3 * F(f_idx, 1) + 1) - element.variables(3 * F(f_idx, 0) + 1),
        element.variables(3 * F(f_idx, 2) + 1) - element.variables(3 * F(f_idx, 0) + 1),
        element.variables(3 * F(f_idx, 1) + 2) - element.variables(3 * F(f_idx, 0) + 2),
        element.variables(3 * F(f_idx, 2) + 2) - element.variables(3 * F(f_idx, 0) + 2);

    Eigen::Vector3<T> theta2;
    theta2 << element.variables(3 * mesh.nVertices() + F(f_idx, 0)),
        element.variables(3 * mesh.nVertices() + F(f_idx, 1)), element.variables(3 * mesh.nVertices() + F(f_idx, 2));

    double dA = 0.5 / MrInv[f_idx].determinant();

    // Compute deformation gradient
    Eigen::Matrix<T, 3, 2> F = M * (MrInv[f_idx] * Eigen::Rotation2D<T>(theta1[f_idx]));

    Eigen::Matrix2<T> J;
    J << 0, 1, 1, 0;
    J *= 0.5 * theta2.sum() / 3 * (lambda2 * lambda2 - lambda1 * lambda1);

    J(0, 0) += -deltaLambda * lambda1;

    Eigen::Matrix2<T> RJRt = Eigen::Rotation2D<T>(theta1[f_idx]) * J * Eigen::Rotation2D<T>(-theta1[f_idx]) / thickness;

    // Compute normal
    Eigen::Vector3<T> n = M.col(0).cross(M.col(1));

    Eigen::Matrix3<T> L = Eigen::Matrix3<T>::Zero();
    Face f = mesh.face(f_idx);
    for(Halfedge he: f.adjacentHalfedges())
    {
      if(he.edge().isBoundary())
        continue;
      
      // rotate edge e around n
      Eigen::Vector3<T> e;
      e << element.variables(3 * geometry.vertexIndices[he.next().vertex()]) -
               element.variables(3 * geometry.vertexIndices[he.vertex()]),
          element.variables(3 * geometry.vertexIndices[he.next().vertex()] + 1) -
              element.variables(3 * geometry.vertexIndices[he.vertex()] + 1),
          element.variables(3 * geometry.vertexIndices[he.next().vertex()] + 2) -
              element.variables(3 * geometry.vertexIndices[he.vertex()] + 2);

      // compute dihedral angle
      Eigen::Vector3<T> nf;
      nf << element.variables(3 * geometry.vertexIndices[he.twin().next().next().vertex()]) -
                element.variables(3 * geometry.vertexIndices[he.vertex()]),
          element.variables(3 * geometry.vertexIndices[he.twin().next().next().vertex()] + 1) -
              element.variables(3 * geometry.vertexIndices[he.vertex()] + 1),
          element.variables(3 * geometry.vertexIndices[he.twin().next().next().vertex()] + 2) -
              element.variables(3 * geometry.vertexIndices[he.vertex()] + 2);
      nf = nf.cross(e);
      T theta = atan2(n.cross(nf).dot(e), e.norm() * nf.dot(n));

      Eigen::Vector3<T> t = n.cross(e);

      // add edge contribution
      L += theta * t.normalized() * t.transpose();
    }
    L /= n.squaredNorm();

    Eigen::Matrix2d S;
    S << 1 / lambda1 / lambda1, 0, 0, 1 / lambda2 / lambda2;
    Eigen::Matrix2<T> E = S * (F.transpose() * L * F - RJRt);

    Eigen::Vector3<T> EVoigt(E(0, 0), E(1, 1), E(0, 1) + E(1, 0));
    // clang-format off
    Eigen::Matrix3d C;
    C <<  E1, 0.4 * std::sqrt(E1 * E2), 0,
          0.4 * std::sqrt(E1 * E2), E2, 0,
          0, 0, 0.6 / 2 * std::sqrt(E1 * E2);
    // clang-format on
    return std::pow(thickness, 2) * EVoigt.dot(C * EVoigt) * dA / 3;
  });
   
  return func;
}
