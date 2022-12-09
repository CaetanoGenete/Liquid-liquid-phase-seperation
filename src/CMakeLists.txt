set(LLPS_OUTFILE_DEFINE "LLPS_OUTPUT_DIR=\"${CMAKE_SOURCE_DIR}/out/\"")

add_executable(gen_fd_error_data "generate_fd_error_data.cpp")
target_link_libraries(gen_fd_error_data PRIVATE LLPS_BASIC)
target_compile_definitions(gen_fd_error_data PRIVATE ${LLPS_OUTFILE_DEFINE})

add_executable(simulate_modelb_fd "modelb.cpp" "_modelb_common.hpp")
target_link_libraries(simulate_modelb_fd PRIVATE LLPS_BASIC)
target_compile_definitions(simulate_modelb_fd PRIVATE ${LLPS_OUTFILE_DEFINE})
#Temporary remove at some point
target_compile_options(simulate_modelb_fd PRIVATE "/arch:AVX512")

add_executable(coupled_modelb_fd "coupled_modelb.cpp" "_modelb_common.hpp")
target_link_libraries(coupled_modelb_fd PRIVATE LLPS_BASIC)
target_compile_definitions(coupled_modelb_fd PRIVATE ${LLPS_OUTFILE_DEFINE})
#Temporary remove at some point
target_compile_options(coupled_modelb_fd PRIVATE "/arch:AVX512")

add_executable(test_view "test_view.cpp" "_modelb_common.hpp")
target_link_libraries(test_view PRIVATE LLPS_BASIC)
target_compile_definitions(test_view PRIVATE ${LLPS_OUTFILE_DEFINE})

if(LLPS_USE_MKL)
    add_executable(gen_spectral_error_data "gen_spectral_error_data.cpp")
    target_link_libraries(gen_spectral_error_data PRIVATE LLPS_MKL)
    target_compile_definitions(gen_spectral_error_data PRIVATE ${LLPS_OUTFILE_DEFINE})

    add_executable(simulate_modelb_spectral "modelb_spectral.cpp" "_modelb_common.hpp")
    target_link_libraries(simulate_modelb_spectral PRIVATE LLPS_MKL)
    target_compile_definitions(simulate_modelb_spectral PRIVATE ${LLPS_OUTFILE_DEFINE})
    #Temporary remove at some point
target_compile_options(simulate_modelb_spectral PRIVATE "/arch:AVX512")
endif()