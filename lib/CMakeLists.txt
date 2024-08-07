cmake_minimum_required(VERSION 3.18)

# Dependencies
include(FetchContent)

# eigen
FetchContent_Declare(Eigen EXCLUDE_FROM_ALL
  URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
  SOURCE_DIR "${SHRINK_MORPH_ROOT}/deps/eigen-src"
)
# libigl
FetchContent_Declare(libigl
  URL https://github.com/libigl/libigl/archive/refs/tags/v2.5.0.zip
  SOURCE_DIR "${SHRINK_MORPH_ROOT}/deps/libigl-src"
)
# TinyAD
FetchContent_Declare(TinyAD
  GIT_REPOSITORY https://github.com/patr-schm/TinyAD.git
  GIT_TAG 81fab13c3884b787c0f03bbbbb95b0f794e54434
  SOURCE_DIR "${SHRINK_MORPH_ROOT}/deps/TinyAD-src"
)
# geometry-central
FetchContent_Declare(geometry-central EXCLUDE_FROM_ALL
  URL https://github.com/DavidJourdan/geometry-central/releases/download/v0.1.1/src.zip
  SOURCE_DIR "${SHRINK_MORPH_ROOT}/deps/geometry-central-src"
  BINARY_DIR "${SHRINK_MORPH_ROOT}/deps/geometry-central-build"
)
FetchContent_MakeAvailable(eigen libigl TinyAD geometry-central)

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  set(EXTRA_LIBS ${EXTRA_LIBS} OpenMP::OpenMP_CXX)
endif(OpenMP_CXX_FOUND)

# MKL PARDISO 
set(MKL_LINK dynamic)
set(MKL_INTERFACE_FULL intel_lp64)
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
	set(MKL_THREADING intel_thread)
else()
	set(MKL_THREADING gnu_thread)
endif()	
find_package(MKL CONFIG)

if(MKL_FOUND)
  set(EXTRA_LIBS ${EXTRA_LIBS} MKL::MKL)
else()
  # SUITESPARSE & UMFPACK
  find_package(SuiteSparse QUIET)
  if (SUITESPARSE_FOUND)
    string(REGEX REPLACE "cholmod" "umfpack" UMFPACK_LIBRARY ${CHOLMOD_LIBRARY})
    if(EXISTS ${UMFPACK_LIBRARY})
      list(APPEND SUITESPARSE_LIBRARIES ${UMFPACK_LIBRARY})
    else()
      set(SUITESPARSE_FOUND FALSE)
    endif()
  endif (SUITESPARSE_FOUND)

  if (SUITESPARSE_FOUND)
    set(EXTRA_LIBS ${EXTRA_LIBS} ${SUITESPARSE_LIBRARIES})
  endif (SUITESPARSE_FOUND)
endif(MKL_FOUND)

add_library(shrink_morph 
  ${SHRINK_MORPH_LIBRARY}/functions.cpp
  ${SHRINK_MORPH_LIBRARY}/LocalGlobalSolver.cpp
  ${SHRINK_MORPH_LIBRARY}/newton.cpp
  ${SHRINK_MORPH_LIBRARY}/parameterization.cpp
  ${SHRINK_MORPH_LIBRARY}/path_extraction.cpp
  ${SHRINK_MORPH_LIBRARY}/stripe_patterns.cpp
  ${SHRINK_MORPH_LIBRARY}/stretch_angles.cpp
  ${SHRINK_MORPH_LIBRARY}/simulation_utils.cpp)

target_link_libraries(shrink_morph PUBLIC Eigen3::Eigen geometry-central TinyAD igl::core)
target_include_directories(shrink_morph PUBLIC ${SHRINK_MORPH_LIBRARY})

# Compiler definitions and options
if(MKL_FOUND)
	target_compile_definitions(shrink_morph PUBLIC EIGEN_USE_MKL_ALL)
	target_compile_definitions(shrink_morph PUBLIC USE_PARDISO)
	target_compile_options(shrink_morph PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
	target_include_directories(shrink_morph PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
endif(MKL_FOUND)

if(SUITESPARSE_FOUND)
	target_compile_definitions(shrink_morph PUBLIC USE_SUITESPARSE)
	target_include_directories(shrink_morph PUBLIC ${SUITESPARSE_INCLUDE_DIRS})
endif(SUITESPARSE_FOUND)
