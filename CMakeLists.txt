# Require a recent version of cmake
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

SET(HPX_IGNORE_COMPILER_COMPATIBILITY 1)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# This project is C++ based.
project(your_app CXX)

# Instruct cmake to find the HPX settings
find_package(HPX REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/generator")

# build your application using HPX
#add_hpx_component(distributed_control_component
#    ESSENTIAL
#    SOURCES dc-sssp.cpp
#    HEADERS distributed_control.hpp
#    COMPONENT_DEPENDENCIES iostreams)

# build your application using HPX
add_hpx_executable(dc_sssp
    ESSENTIAL
    SOURCES dc_sssp.cpp
    HEADERS distributed_control.hpp
    COMPONENT_DEPENDENCIES iostreams)
