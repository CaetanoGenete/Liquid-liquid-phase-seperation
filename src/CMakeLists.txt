set(LLPS_OUTFILE_DEFINE "LLPS_OUTPUT_DIR=\"${CMAKE_SOURCE_DIR}/out/\"")

add_executable(gen_fd_error_data "generate_fd_error_data.cpp")
target_link_libraries(gen_fd_error_data PRIVATE LLPS_BASIC)
target_compile_definitions(gen_fd_error_data PRIVATE ${LLPS_OUTFILE_DEFINE})

add_executable(simulate_modelb_fd "modelb.cpp")
target_link_libraries(simulate_modelb_fd PRIVATE LLPS_BASIC)
target_compile_definitions(simulate_modelb_fd PRIVATE ${LLPS_OUTFILE_DEFINE})

if(LLPS_USE_MKL)
    add_executable(gen_spectral_error_data "gen_spectral_error_data.cpp")
    target_link_libraries(gen_spectral_error_data PRIVATE LLPS_MKL)
    target_compile_definitions(gen_spectral_error_data PRIVATE ${LLPS_OUTFILE_DEFINE})

    add_executable(simulate_modelb_spectral "modelb.cpp")
    target_link_libraries(simulate_modelb_spectral PRIVATE LLPS_MKL)
    target_compile_definitions(simulate_modelb_spectral PRIVATE ${LLPS_OUTFILE_DEFINE})
endif()