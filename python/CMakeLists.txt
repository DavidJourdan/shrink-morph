cmake_minimum_required(VERSION 3.18)

find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

# nanobind
FetchContent_Declare(nanobind
	GIT_REPOSITORY https://github.com/wjakob/nanobind.git
	GIT_TAG 8d7f1ee0621c17fa370b704b2100ffa0243d5bfb
	SOURCE_DIR "${SHRINK_MORPH_ROOT}/deps/nanobind-src"
)

FetchContent_MakeAvailable(nanobind)

nanobind_add_module(shrink_morph_py
	shrink_morph.cpp)
target_link_libraries(shrink_morph_py PRIVATE Eigen3::Eigen geometry-central TinyAD igl::core shrink_morph)

set(CMAKE_INSTALL_PREFIX ${SHRINK_MORPH_LIBRARY})
install(TARGETS shrink_morph_py
  LIBRARY 
  DESTINATION "${SHRINK_MORPH_ROOT}/python")