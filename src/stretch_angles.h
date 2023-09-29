#pragma once

#include <geometrycentral/surface/manifold_surface_mesh.h>

geometrycentral::surface::FaceData<double> computeStretchAngles(geometrycentral::surface::ManifoldSurfaceMesh& mesh,
                                                                const Eigen::MatrixXd& V,
                                                                const Eigen::MatrixXd& P,
                                                                const Eigen::MatrixXi& F);

geometrycentral::surface::FaceData<double>
computeStretchAngles(geometrycentral::surface::ManifoldSurfaceMesh& mesh,
                     const Eigen::MatrixXd& V,
                     const Eigen::MatrixXi& F,
                     const geometrycentral::surface::FaceData<Eigen::Matrix2d>& Minv);
