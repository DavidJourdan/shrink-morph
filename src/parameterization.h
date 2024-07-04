/*
 * This file is part of TinyAD and released under the MIT license.
 * Author: Patrick Schmidt
 */
#pragma once

#include "LocalGlobalSolver.h"

#include <TinyAD/Support/GeometryCentral.hh>
#include <TinyAD/ScalarFunction.hh>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

Eigen::MatrixXd localGlobal(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double lambda1, double lambda2);

geometrycentral::surface::FaceData<Eigen::Matrix2d>
precomputeParamData(geometrycentral::surface::VertexPositionGeometry& geometry);

geometrycentral::surface::FaceData<Eigen::Matrix2d>
precomputeSimData(geometrycentral::surface::ManifoldSurfaceMesh& mesh,
                  const Eigen::MatrixXd& P,
                  const Eigen::MatrixXi& F);

geometrycentral::surface::EdgeData<double>
computeDualCotanWeights(geometrycentral::surface::IntrinsicGeometryInterface& geometry);

void subdivideMesh(geometrycentral::surface::VertexPositionGeometry& geometry,
                   Eigen::MatrixXd& V,
                   Eigen::MatrixXd& P,
                   Eigen::MatrixXi& F,
                   std::vector<Eigen::SparseMatrix<double>>& subdivMat,
                   double threshold);

/**
 * Compute tutte embedding with boundary on circle.
 * Per-vertex 2D coordinates returned as geometrycentral VertexData.
 */
Eigen::MatrixXd tutte_embedding(const Eigen::MatrixXd& _V, const Eigen::MatrixXi& _F);

TinyAD::ScalarFunction<2, double, geometrycentral::surface::VertexRangeF::Etype>
parameterizationFunction(geometrycentral::surface::VertexPositionGeometry& geometry,
                         double wPhi,
                         double lambda1,
                         double lambda2);

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>
computeSVDdata(const Eigen::MatrixXd& V, const Eigen::MatrixXd& P, const Eigen::MatrixXi& F);