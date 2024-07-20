#pragma once

#include <Eigen/Core>
#include <TinyAD/Support/GeometryCentral.hh>
#include <TinyAD/ScalarFunction.hh>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>

#include <functional>

template <class Func, class Solver>
void newton(
    Eigen::VectorXd& x,
    Func& func,
    Solver& solver,
    int max_iters = 1000,
    double lim = 1e-6,
    bool verbose = true,
    const std::vector<int>& fixedIdx = {},
    const std::function<void(const Eigen::VectorXd&)>& callBack = [](const auto&) {});

template <class Func>
void newton(
    geometrycentral::surface::IntrinsicGeometryInterface& geometry,
    Eigen::MatrixXd& V,
    Func& func,
    int max_iters,
    double lim,
    bool verbose = true,
    const std::vector<int>& fixedIdx = {},
    const std::function<void(const Eigen::VectorXd&)>& callBack = [](const auto&) {});

Eigen::MatrixXd sparse_gauss_newton(
    geometrycentral::surface::IntrinsicGeometryInterface& geometry,
    const Eigen::MatrixXd& targetV,
    const geometrycentral::surface::FaceData<Eigen::Matrix2d>& MrInv,
    const geometrycentral::surface::FaceData<double>& theta1,
    geometrycentral::surface::VertexData<double>& theta2,
    const TinyAD::ScalarFunction<1, double, Eigen::Index>& adjointFunc,
    const std::vector<int>& fixedIdx,
    int max_iters,
    double lim,
    double wM,
    double wL,
    double E1,
    double lambda1,
    double lambda2,
    double deltaLambda,
    double thickness,
    const std::function<void(const Eigen::VectorXd&)>& callback = [](const auto&) {});