#pragma once

#include <geometrycentral/numerical/linear_solvers.h>
#include <geometrycentral/surface/embedded_geometry_interface.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>

// Implementation of "Stripe Patterns on Surfaces" [Knoppel et al. 2015]

// Takes as input a geometry along with vertex-based frequencies and a line field (2-RoSy) and outputs a 2\pi-periodic
// function defined on triangle corners such that the 0 (mod 2\pi) isolines of this function are stripes following the
// direction field spaced according to the target frequencies

std::tuple<std::vector<geometrycentral::Vector3>, std::vector<std::array<size_t, 2>>>
extractPolylinesFromStripePattern(geometrycentral::surface::EmbeddedGeometryInterface& geometry,
                                  const geometrycentral::surface::CornerData<double>& stripesValues,
                                  const geometrycentral::surface::FaceData<int>& stripesIndices,
                                  const geometrycentral::surface::FaceData<int>& fieldIndices,
                                  const geometrycentral::surface::VertexData<geometrycentral::Vector2>& directionField,
                                  bool connectOnSingularities = false);

// Build a Laplace-like matrix with double entries (necessary to represent complex conjugation)
geometrycentral::SparseMatrix<double>
buildVertexEnergyMatrix(geometrycentral::surface::IntrinsicGeometryInterface& geometry,
                        const geometrycentral::surface::VertexData<geometrycentral::Vector2>& directionField,
                        const geometrycentral::surface::FaceData<int>& branchIndices,
                        const geometrycentral::surface::VertexData<double>& frequencies);

// Build a lumped mass matrix with double entries
geometrycentral::SparseMatrix<double>
computeRealVertexMassMatrix(geometrycentral::surface::IntrinsicGeometryInterface& geometry);

// extract the final texture coordinates from the parameterization
std::tuple<geometrycentral::surface::CornerData<double>, geometrycentral::surface::FaceData<int>>
computeTextureCoordinates(geometrycentral::surface::IntrinsicGeometryInterface& geometry,
                          const geometrycentral::surface::VertexData<geometrycentral::Vector2>& directionField,
                          const geometrycentral::surface::VertexData<double>& frequencies,
                          const geometrycentral::surface::VertexData<geometrycentral::Vector2>& parameterization);
