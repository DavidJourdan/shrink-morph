cmake_minimum_required(VERSION 3.29)

# c++ flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

### Configuration
set(SHRINK_MORPH_ROOT "${CMAKE_CURRENT_LIST_DIR}")
set(SHRINK_MORPH_LIBRARY "${SHRINK_MORPH_ROOT}/lib")

# FindXXX.cmake are in the /cmake directory
list(APPEND CMAKE_MODULE_PATH "${SHRINK_MORPH_ROOT}/cmake") 

add_subdirectory(lib)

option(PYTHON_APP ON)

if(${PYTHON_BINDINGS})
  add_subdirectory(python)
else(${PYTHON_BINDINGS})
  add_subdirectory(app)
endif(${PYTHON_BINDINGS})