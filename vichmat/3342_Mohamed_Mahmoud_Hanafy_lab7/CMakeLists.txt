cmake_minimum_required(VERSION 3.10)
project(gauss)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(Eigen3_DIR "/opt/homebrew/opt/eigen/share/eigen3/cmake")
find_package(Eigen3 REQUIRED)

add_executable(gauss
    main.cpp
)

target_link_libraries(gauss Eigen3::Eigen)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(gauss PRIVATE -O3 -march=native)
endif()