cmake_minimum_required(VERSION 3.29)

project(shrink_morph)

# c++ flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# FindXXX.cmake are in the /cmake directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake") 

### Required dependencies ###

include(FetchContent)

# Look for a system install of Eigen (otherwise dependencies will fetch it)
find_package(Eigen3)

# libigl
FetchContent_Declare(libigl
  URL https://github.com/libigl/libigl/archive/refs/tags/v2.5.0.zip
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl-src"
  BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl-build"
)
# TinyAD
FetchContent_Declare(TinyAD
  GIT_REPOSITORY https://github.com/patr-schm/TinyAD.git
  GIT_TAG 81fab13c3884b787c0f03bbbbb95b0f794e54434
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/TinyAD-src"
  BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/TinyAD-build"
)
# polyscope
FetchContent_Declare(polyscope
  URL https://github.com/DavidJourdan/polyscope/releases/download/v2.2.1/src.zip
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope-src"
  BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope-build"
)
# geometry-central
FetchContent_Declare(geometry-central
  URL https://github.com/DavidJourdan/geometry-central/releases/download/v0.1.1/src.zip
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/geometry-central-src"
  BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/geometry-central-build"
)
FetchContent_MakeAvailable(libigl TinyAD polyscope geometry-central)

### Optional dependencies ###

# OR-Tools
set(ORTOOLS_PATH "/opt/or-tools/") # change this variable to whatever path the binaries have been extracted to
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} "${ORTOOLS_PATH}") 
find_package(ortools CONFIG)
if(ortools_FOUND)
  message(STATUS "Found OR-Tools")
  add_definitions(-DUSE_ORTOOLS)
  set(EXTRA_LIBS ${EXTRA_LIBS} ortools::ortools)
endif(ortools_FOUND)

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

# cpp files minus the mains
file(GLOB SRCFILES 
  src/functions.cpp
  src/LocalGlobalSolver.cpp
  src/newton.cpp
  src/parameterization.cpp
  src/path_extraction.cpp
  src/save.cpp
  src/stripe_patterns.cpp
  src/stretch_angles.cpp
  src/simulation_utils.cpp)
  
# Create executables
add_executable(shrink_morph  src/main.cpp ${SRCFILES})
add_executable(shrink_morph_cli  src/cli.cpp ${SRCFILES})
set(EXECUTABLES 
  shrink_morph 
  shrink_morph_cli
)
  
foreach(exec ${EXECUTABLES})
  target_link_libraries(${exec} Eigen3::Eigen geometry-central polyscope TinyAD igl::core nlohmann_json::nlohmann_json ${EXTRA_LIBS})
    
  # Compiler definitions and options
  if(MKL_FOUND)
    target_compile_definitions(${exec} PUBLIC EIGEN_USE_MKL_ALL)
    target_compile_definitions(${exec} PUBLIC USE_PARDISO)
    target_compile_options(${exec} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(${exec} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  endif(MKL_FOUND)

  if(SUITESPARSE_FOUND)
    target_compile_definitions(${exec} PUBLIC USE_SUITESPARSE)
    target_include_directories(${exec} PUBLIC ${SUITESPARSE_INCLUDE_DIRS})
  endif(SUITESPARSE_FOUND)

  target_compile_definitions(${exec} PUBLIC DATA_PATH_STR="${CMAKE_CURRENT_SOURCE_DIR}/data/")
  # target_compile_options(${exec} PRIVATE $<$<CONFIG:RELEASE>:
  #   $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Ofast>
  #   $<$<CXX_COMPILER_ID:MSVC>:/Ox>>)

  # remove annoying Polyscope-related compiler warning
  target_compile_options(${exec} PRIVATE 
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Wno-nonnull>)
endforeach(exec ${EXECUTABLES})