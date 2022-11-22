cmake_minimum_required(VERSION 3.19)

project(phase_separation
    VERSION 0.0.1
    LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
set(LLPS_DEPENDENCIES_SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/")
set(LLPS_DEPENDENCIES_BINARY_DIR "${CMAKE_BINARY_DIR}/dependencies/")

set(LLPS_HEADERS
    "include/grid.hpp"
    "include/calculus/finite_difference.hpp"
    "include/calculus/differentiate.hpp"
    "include/utilities/io.hpp"
    "include/utilities/data_analytics.hpp"
    "include/utilities/meta.hpp")
source_group(TREE ${CMAKE_SOURCE_DIR} FILES ${LLPS_HEADERS})

add_library(LLPS_BASIC INTERFACE ${LLPS_HEADERS})
target_include_directories(LLPS_BASIC INTERFACE "${CMAKE_SOURCE_DIR}/include/")
target_compile_features(LLPS_BASIC INTERFACE cxx_std_20)

#Optional features
option(LLPS_BUILD_TESTS "Builds and runs tests.")
option(LLPS_USE_EIGEN "Uses eigen arrays.")
option(LLPS_USE_MKL "Use MKL FFT.")

if(LLPS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

#Adding source files
add_subdirectory(src)
