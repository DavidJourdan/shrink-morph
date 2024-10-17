#pragma once

#include "solvers.h"

#include <geometrycentral/surface/vertex_position_geometry.h>

std::vector<std::vector<geometrycentral::Vector3>>
orderPolylines(const std::vector<std::vector<geometrycentral::Vector3>>& isolines);

std::vector<std::vector<geometrycentral::Vector3>>
simplifyPolylines(const std::vector<std::vector<geometrycentral::Vector3>>& polylines, double z = 0);

std::vector<std::vector<std::vector<geometrycentral::Vector3>>>
generatePaths(geometrycentral::surface::EmbeddedGeometryInterface& geometry,
              const Eigen::VectorXd& theta1,
              const Eigen::VectorXd& theta2,
              double layerHeight,
              int nLayers,
              double spacing);

std::vector<std::vector<geometrycentral::Vector3>> generateOneLayer(geometrycentral::surface::EmbeddedGeometryInterface& geometry,
                                                                    const Eigen::VectorXd& theta1,
                                                                    const Eigen::VectorXd& theta2,
                                                                    const Eigen::SparseMatrix<double>& massMatrix,
                                                                    Eigen::VectorXd& u,
                                                                    LDLTSolver& solver,
                                                                    int i,
                                                                    int nLayers,
                                                                    double layerHeight,
                                                                    double spacing);

void writePaths(const std::string& filename,
                const std::vector<std::vector<geometrycentral::Vector3>>& paths,
                double height);

std::vector<std::vector<geometrycentral::Vector3>> edgeToPolyline(const std::vector<geometrycentral::Vector3>& points,
                                                                  const std::vector<std::array<size_t, 2>>& edges);

std::tuple<geometrycentral::surface::VertexData<double>,
           geometrycentral::surface::FaceData<double>,
           geometrycentral::surface::EdgeData<double>>
hodgeDecomposition(geometrycentral::surface::VertexPositionGeometry& geometry,
                   const geometrycentral::surface::EdgeData<double>& oneForm);

Eigen::VectorXd
curlFreeParameterization(geometrycentral::surface::VertexPositionGeometry& geometry,
                         const geometrycentral::surface::FaceData<geometrycentral::Vector3>& vectorField);

Eigen::VectorXd
curlFreeParameterization(geometrycentral::surface::VertexPositionGeometry& geometry,
                         const geometrycentral::surface::FaceData<geometrycentral::Vector3>& vectorField,
                         LDLTSolver& solver);

Eigen::VectorXd
curlFreeParameterization(const Eigen::MatrixXd& P, const Eigen::MatrixXi& F, const Eigen::VectorXd& theta);

std::vector<std::vector<geometrycentral::Vector3>>
orderPolylinesNew(const std::vector<std::vector<geometrycentral::Vector3>>& isolines,
                  const Eigen::MatrixX2d& P,
                  const Eigen::MatrixXi& F,
                  const Eigen::VectorXd& theta);

std::vector<std::vector<geometrycentral::Vector3>>
orderPolylinesNew(const std::vector<std::vector<geometrycentral::Vector3>>& isolines,
                  const Eigen::MatrixX2d& P,
                  const Eigen::VectorXd& vCoordinate);