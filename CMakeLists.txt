cmake_minimum_required(VERSION 3.16)

project(shrink_morph)

# c++ flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

### Required dependencies ###

include(FetchContent)

# libigl
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG fd16e24391f3a85628999fe67b84a9b2201f2bc2
)

# TinyAD
FetchContent_Declare(
    TinyAD
    GIT_REPOSITORY https://github.com/patr-schm/TinyAD.git
    GIT_TAG 81fab13c3884b787c0f03bbbbb95b0f794e54434
)

# polyscope
FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG c772fd9bf05f4d914fbc74fd5f1a69af58687a87
)

# geometry-central
FetchContent_Declare(
    geometry-central
    GIT_REPOSITORY https://github.com/nmwsharp/geometry-central.git
    GIT_TAG b467ce80c4be13b94160f0f579cc71500eae5663
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
  message(STATUS "Found OpenMP")
  set(EXTRA_LIBS ${EXTRA_LIBS}  ${OpenMP_CXX_LIBRARIES})
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
endif(MKL_FOUND)

# cpp files minus the mains
file(GLOB SRCFILES 
  src/functions.cpp
  src/LocalGlobalSolver.cpp
  src/newton.cpp
  src/parameterization.cpp
  src/path_extraction.cpp
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
  target_link_libraries(${exec} Eigen3::Eigen geometry-central polyscope TinyAD igl::core ${EXTRA_LIBS})
    
  # Compiler definitions and options
  if(MKL_FOUND)
    target_compile_definitions(${exec} PUBLIC EIGEN_USE_MKL_ALL)
    target_compile_definitions(${exec} PUBLIC USE_PARDISO)
    target_compile_options(${exec} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(${exec} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  endif(MKL_FOUND)

  target_compile_definitions(${exec} PUBLIC DATA_PATH_STR="${CMAKE_CURRENT_SOURCE_DIR}/data/")
  target_compile_options(${exec} PRIVATE $<$<CONFIG:RELEASE>:
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Ofast>
    $<$<CXX_COMPILER_ID:MSVC>:/Ox>>)

  # reomove annoying Polyscope-related compiler warning
  target_compile_options(${exec} PRIVATE 
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Wno-nonnull>)
endforeach(exec ${EXECUTABLES})