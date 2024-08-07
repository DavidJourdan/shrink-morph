cmake_minimum_required(VERSION 3.29)

option(BUILD_GUI ON)

if(${BUILD_GUI})

  # polyscope
  FetchContent_Declare(polyscope
    URL https://github.com/DavidJourdan/polyscope/releases/download/v2.2.1/src.zip
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope-src"
    BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope-build"
  )
  FetchContent_MakeAvailable(polyscope)


  add_executable(shrink_morph_app  main.cpp save.cpp generate_trajectories.cpp)
  target_link_libraries(shrink_morph_app shrink_morph polyscope nlohmann_json::nlohmann_json)
  target_compile_definitions(shrink_morph_app PUBLIC DATA_PATH_STR="${SHRINK_MORPH_ROOT}/data/")
  
  # remove annoying Polyscope-related compiler warning
  target_compile_options(shrink_morph_app PRIVATE 
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Wno-nonnull>)

endif(${BUILD_GUI})

add_executable(shrink_morph_cli  cli.cpp)
target_link_libraries(shrink_morph_cli shrink_morph)
target_compile_definitions(shrink_morph_cli PUBLIC DATA_PATH_STR="${SHRINK_MORPH_ROOT}/data/")
