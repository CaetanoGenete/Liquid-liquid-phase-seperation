
add_executable(gen_fd_error_data "generate_fd_error_data.cpp")
target_link_libraries(gen_fd_error_data PRIVATE LLPS_BASIC)
target_compile_definitions(gen_fd_error_data PRIVATE "LLPS_OUTPUT_DIR=\"${CMAKE_SOURCE_DIR}/out/\"")