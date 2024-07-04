#pragma once

#ifdef USE_PARDISO
#include "Eigen/PardisoSupport"
using LDLTSolver = Eigen::PardisoLDLT<Eigen::SparseMatrix<double>>;
using LLTSolver = Eigen::PardisoLLT<Eigen::SparseMatrix<double>>;
using LUSolver = Eigen::PardisoLU<Eigen::SparseMatrix<double>>;
#else
#ifdef USE_SUITESPARSE
#include <Eigen/CholmodSupport>
#include <Eigen/UmfPackSupport>
using LDLTSolver = Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double>>;
using LLTSolver = Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>>;
using LUSolver = Eigen::UmfPackLU<Eigen::SparseMatrix<double>>;
#else
#include<Eigen/SparseCholesky>
#include<Eigen/SparseLU> 
using LDLTSolver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
using LLTSolver = Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>;
using LUSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>>;
#endif // USE_SUITESPARSE
#endif // USE_PARDISO
