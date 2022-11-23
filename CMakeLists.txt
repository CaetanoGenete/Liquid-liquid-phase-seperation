cmake_minimum_required(VERSION 3.19)

include(FetchContent)

project(phase_separation
    VERSION 0.0.1
    LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
set(LLPS_DEPENDENCIES_SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/")
set(LLPS_DEPENDENCIES_BINARY_DIR "${CMAKE_BINARY_DIR}/dependencies/")

set(LLPS_HEADERS
    "include/grid.hpp"
    "include/aligned_allocator.hpp"
    "include/calculus/finite_difference.hpp"
    "include/calculus/differentiate.hpp"
    "include/calculus/fourier_spectral.hpp"
    "include/utilities/io.hpp"
    "include/utilities/data_analytics.hpp"
    "include/utilities/meta.hpp")
source_group(TREE ${CMAKE_SOURCE_DIR} FILES ${LLPS_HEADERS})

add_library(LLPS_BASIC INTERFACE ${LLPS_HEADERS})
target_include_directories(LLPS_BASIC INTERFACE "${CMAKE_SOURCE_DIR}/include/")
target_compile_features(LLPS_BASIC INTERFACE cxx_std_20)

#Look for Boost
find_package(Boost 1.80.0 QUIET)
if(Boost_FOUND)
    target_include_directories(modelB PRIVATE ${Boost_INCLUDE_DIRS})
#If Boost is not found, download only odeint and integrate into the project
else()
    FetchContent_Declare(odeint
        GIT_REPOSITORY https://github.com/boostorg/odeint.git
        GIT_TAG        boost-1.80.0
        SOURCE_DIR     "${LLPS_DEPENDENCIES_SOURCE_DIR}boost/odeint"
        BINARY_DIR     "${LLPS_DEPENDENCIES_BINARY_DIR}boost/odeint/build/"
        SUBBUILD_DIR   "${LLPS_DEPENDENCIES_BINARY_DIR}boost/odeint/sub-build/")

    FetchContent_MakeAvailable(odeint)
    target_include_directories(LLPS_BASIC INTERFACE "${odeint_SOURCE_DIR}/include/")
endif()

#Optional features
option(LLPS_BUILD_TESTS "Builds and runs tests.")
option(LLPS_USE_EIGEN "Uses eigen arrays.")
option(LLPS_USE_MKL "Use MKL FFT.")

if(LLPS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(LLPS_USE_MKL)
    find_package(MKL CONFIG)
    if(MKL_FOUND)
        add_library(LLPS_MKL INTERFACE)

        target_compile_options(LLPS_MKL INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
        target_include_directories(LLPS_MKL INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
        target_compile_definitions(LLPS_MKL INTERFACE "LLPS_USE_MKL")
        target_link_libraries(LLPS_MKL INTERFACE LLPS_BASIC INTERFACE $<LINK_ONLY:MKL::MKL>)
    else()
        message(WARNING "intelMKL was not found! Disabling usage.")
    endif()
endif()

#Adding source files
add_subdirectory(src)
