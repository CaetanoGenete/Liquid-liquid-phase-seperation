
add_executable(gen_fd_error_data "generate_fd_error_data.cpp")
target_link_libraries(gen_fd_error_data PRIVATE LLPS_BASIC)
target_compile_definitions(gen_fd_error_data PRIVATE "LLPS_OUTPUT_DIR=\"${CMAKE_SOURCE_DIR}/out/\"")

if(LLPS_USE_MKL)
    add_executable(fft_test "fft_test.cpp")
    target_link_libraries(fft_test PRIVATE LLPS_MKL)
    target_compile_definitions(fft_test PRIVATE "LLPS_OUTPUT_DIR=\"${CMAKE_SOURCE_DIR}/out/\"")
endif()