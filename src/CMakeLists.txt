set(LLPS_OUTFILE_DEFINE "LLPS_OUTPUT_DIR=\"${CMAKE_SOURCE_DIR}/out/\"")

function(llps_add_executable target_name link_target)
    add_executable(${target_name} ${ARGN})
    target_link_libraries(${target_name} PRIVATE ${link_target})
    target_compile_definitions(${target_name} PRIVATE ${LLPS_OUTFILE_DEFINE})
endfunction()

llps_add_executable(gen_fd_error_data  LLPS_BASIC "generate_fd_error_data.cpp")
llps_add_executable(simulate_modelb_fd LLPS_BASIC "modelb.cpp" "_modelb_common.hpp")
llps_add_executable(coupled_modelb_fd  LLPS_BASIC "coupled_modelb.cpp" "_modelb_common.hpp")
llps_add_executable(coupled_modelb_diffusion  LLPS_BASIC "coupled_modelb_diffusion.cpp" "_modelb_common.hpp" "multi_range_algebra.hpp") 

llps_add_executable(test_view  LLPS_BASIC "test_view.cpp" "_modelb_common.hpp" "multi_range_algebra.hpp")

if(LLPS_USE_MKL)
    llps_add_executable(gen_spectral_error_data  LLPS_MKL "gen_spectral_error_data.cpp")
    llps_add_executable(simulate_modelb_spectral LLPS_MKL "modelb_spectral.cpp" "_modelb_common.hpp")
endif()