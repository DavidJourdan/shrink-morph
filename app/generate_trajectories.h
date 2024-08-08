#pragma once

#include <geometrycentral/surface/vertex_position_geometry.h>

void generateTrajectories(geometrycentral::surface::VertexPositionGeometry& geometry,
                          const Eigen::MatrixXd& _V,
                          const Eigen::MatrixXd& _P,
                          const Eigen::MatrixXi& _F,
                          geometrycentral::surface::VertexData<double>& vTheta2,
                          std::string filename,
                          double layerHeight = 0.08,
                          double spacing = 0.4,
                          int nLayers = 10);

void drawPathsAndTravels(const std::vector<std::vector<geometrycentral::Vector3>>& paths, double spacing, int id);
