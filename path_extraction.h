#pragma once

#include <geometrycentral/surface/vertex_position_geometry.h>

std::vector<std::vector<geometrycentral::Vector3>>
orderPolylines(const std::vector<std::vector<geometrycentral::Vector3>>& isolines,
               geometrycentral::surface::EmbeddedGeometryInterface& geometry,
               double timeLimit = 5);

std::vector<std::vector<geometrycentral::Vector3>>
simplifyPolylines(const std::vector<std::vector<geometrycentral::Vector3>>& polylines, double z = 0);

std::vector<std::vector<std::vector<geometrycentral::Vector3>>>
generatePaths(geometrycentral::surface::EmbeddedGeometryInterface& geometry,
              const Eigen::VectorXd& theta1,
              const Eigen::VectorXd& theta2,
              double layerHeight,
              int nLayers,
              double spacing,
              double timeLimit = 5);

void writePaths(const std::string& filename,
                const std::vector<std::vector<geometrycentral::Vector3>>& paths,
                double height);

void drawPathsAndTravels(const std::vector<std::vector<geometrycentral::Vector3>>& paths, double spacing, int id);

void stripePattern(geometrycentral::surface::VertexPositionGeometry& geometry,
                   const Eigen::MatrixXd& _V,
                   Eigen::MatrixXd& _P,
                   const Eigen::MatrixXi& _F,
                   geometrycentral::surface::VertexData<double>& vTheta2,
                   std::string filename,
                   double timeLimit = 5,
                   double layerHeight = 0.08,
                   double spacing = 0.4,
                   int nLayers = 10);

std::vector<std::vector<geometrycentral::Vector3>> edgeToPolyline(const std::vector<geometrycentral::Vector3>& points,
                                                                  const std::vector<std::array<size_t, 2>>& edges);
