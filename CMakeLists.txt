cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(Fluid2D LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()


add_executable(Fluid2D)

target_sources(Fluid2D PRIVATE src/LevelSetSim.cpp src/main.cpp)

find_package(Eigen3 REQUIRED)

target_include_directories(Fluid2D PUBLIC include PRIVATE ${EIGEN3_INCLUDE_DIR})