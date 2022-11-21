include(FetchContent)
include(GoogleTest)
include(add_gtest)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        release-1.12.1
    SOURCE_DIR     "${LLPS_DEPENDENCIES_SOURCE_DIR}/googletest/"
    BINARY_DIR     "${LLPS_DEPENDENCIES_BINARY_DIR}/googletest/build/"
    SUBBUILD_DIR   "${LLPS_DEPENDENCIES_BINARY_DIR}/googletest/sub-build/")

FetchContent_MakeAvailable(googletest)

set_target_properties(gtest      PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)
set_target_properties(gmock      PROPERTIES FOLDER extern)
set_target_properties(gmock_main PROPERTIES FOLDER extern)

add_gtest(test_finite_difference "test_finite_difference.cpp" LLPS_BASIC) 
