#include "stretch_angles.h"

#include <TinyAD/Utils/Helpers.hh>

geometrycentral::surface::FaceData<double>
computeStretchAngles(geometrycentral::surface::ManifoldSurfaceMesh& mesh,
                     const Eigen::MatrixXd& V,
                     const Eigen::MatrixXi& F,
                     const geometrycentral::surface::FaceData<Eigen::Matrix2d>& Minv)
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  FaceData<double> angles(mesh);

#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads() - 1)
  for(int i = 0; i < F.rows(); ++i)
  {
    // Get 3D vertex positions
    Eigen::Vector3d ar = V.row(F(i, 0));
    Eigen::Vector3d br = V.row(F(i, 1));
    Eigen::Vector3d cr = V.row(F(i, 2));
    Eigen::Matrix<double, 3, 2> Mr = TinyAD::col_mat(br - ar, cr - ar);

    // Compute deformation gradient
    Eigen::Matrix<double, 3, 2> F = Mr * Minv[i];

    // Compute first fundamental form
    Eigen::Matrix2d I = F.transpose() * F;

    // when considering the metric tensor I (instead of the differential J), the angles should be multiplied by two
    angles[i] = atan2(I(0, 1) + I(1, 0), I(0, 0) - I(1, 1));
    if(angles[i] < 0)
      angles[i] += PI;
    else
      angles[i] -= PI;
    angles[i] /= 2;
  }

  return angles;
}

geometrycentral::surface::FaceData<double> computeStretchAngles(geometrycentral::surface::ManifoldSurfaceMesh& mesh,
                                                                const Eigen::MatrixXd& V,
                                                                const Eigen::MatrixXd& P,
                                                                const Eigen::MatrixXi& F)
{
  using namespace geometrycentral;
  using namespace geometrycentral::surface;

  FaceData<double> angles(mesh);

#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads() - 1)
  for(int i = 0; i < F.rows(); ++i)
  {
    // Get 3D vertex positions
    Eigen::Vector3d ar = V.row(F(i, 0));
    Eigen::Vector3d br = V.row(F(i, 1));
    Eigen::Vector3d cr = V.row(F(i, 2));
    Eigen::Matrix<double, 3, 2> Mr = TinyAD::col_mat(br - ar, cr - ar);

    // Get 2D vertex positions
    Eigen::Vector2d a = P.row(F(i, 0));
    Eigen::Vector2d b = P.row(F(i, 1));
    Eigen::Vector2d c = P.row(F(i, 2));
    Eigen::Matrix2d M = TinyAD::col_mat(b - a, c - a);

    // Compute deformation gradient
    Eigen::Matrix<double, 3, 2> F = Mr * M.inverse();

    // Compute first fundamental form
    Eigen::Matrix2d I = F.transpose() * F;

    // when considering the metric tensor I (instead of the differential J), the angles should be multiplied by two
    angles[i] = atan2(I(0, 1) + I(1, 0), I(0, 0) - I(1, 1));
    if(angles[i] < 0)
      angles[i] += PI;
    else
      angles[i] -= PI;
    angles[i] /= 2;
  }

  return angles;
}
