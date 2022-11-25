include(FetchContent)

function(include_boost_ncmake target_name link_type lib_name)
    FetchContent_Declare(${lib_name}
        GIT_REPOSITORY "https://github.com/boostorg/${lib_name}.git"
        GIT_TAG        "boost-1.80.0"
        SOURCE_DIR     "${LLPS_DEPENDENCIES_SOURCE_DIR}boost/${lib_name}"
        BINARY_DIR     "${LLPS_DEPENDENCIES_BINARY_DIR}boost/${lib_name}/build/"
        SUBBUILD_DIR   "${LLPS_DEPENDENCIES_BINARY_DIR}boost/${lib_name}/sub-build/")

    FetchContent_MakeAvailable(${lib_name})
    target_include_directories(${target_name} ${link_type} "${${lib_name}_SOURCE_DIR}/include/") 
endfunction()

function(include_boost_cmake target_name link_type lib_name)
    FetchContent_Declare(${lib_name}
        GIT_REPOSITORY "https://github.com/boostorg/${lib_name}.git"
        GIT_TAG        "boost-1.80.0"
        SOURCE_DIR     "${LLPS_DEPENDENCIES_SOURCE_DIR}boost/${lib_name}"
        BINARY_DIR     "${LLPS_DEPENDENCIES_BINARY_DIR}boost/${lib_name}/build/"
        SUBBUILD_DIR   "${LLPS_DEPENDENCIES_BINARY_DIR}boost/${lib_name}/sub-build/")

    FetchContent_MakeAvailable(${lib_name})
    target_link_libraries(${target_name} ${link_type} "Boost::${lib_name}") 
endfunction()