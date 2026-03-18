
# to build with cmake
# create a build directory and move into it
# $ mkdir build
# $ cd build
# generate the makefile (to do only ones, if we don't add files or change makefiles)
# don't forget the two points at the end of the command '..'.
# It runs cmake in the 'build' directory
# but with the data from the '..' directory.
# the script should autodetect the CUDA architecture, when run alone
# $ cmake ..
# You can specify the architecture with the SM variable:
# $ cmake -DSM=30 ..

cmake_minimum_required (VERSION 3.15 FATAL_ERROR)

OPTION(USE_PYTHON "Compile BeAR's Python Wrapper" OFF)
OPTION(DEBUG_CUDA "Enable full CUDA synchronization for debugging" OFF)

# check command line args, user can specify SM to override manual setting
set(SM "0" CACHE STRING "GPU SM value")


project(bear C CXX CUDA Fortran)

# if (USE_PYBEAR MATCHES ON)
#   set(USE_PYTHON ON)
# endif()


#some C++/C flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -march=native -pedantic -MMD")
set(CMAKE_C_FLAGS "-Wall -march=native -pedantic -MMD")


set(CMAKE_BUILD_TYPE Release)



if(SM MATCHES "0")
  set(CUDA_ARCH_FLAGS "-arch=sm_75")
else()
  message(STATUS "CUDA Architecture manually set to: -arch=sm_${SM}")
  set(CUDA_ARCH_FLAGS "-arch=sm_${SM}")
endif()
  
  
#output directories
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})



#set CUDA flags
string(APPEND CMAKE_CUDA_FLAGS ${CUDA_ARCH_FLAGS})
string(APPEND CMAKE_CUDA_FLAGS " -std=c++17 -march=native -lineinfo -Xptxas -v")

#string (APPEND CMAKE_CUDA_FLAGS " -cudart shared" )



#define the source files 
set(SRC
  src/additional/aux_functions.cpp
  src/additional/piecewise_poly.cpp
  src/config/global_config.cpp
  src/forward_model/forward_model.cpp
  src/forward_model/generic_config.cpp
  src/forward_model/emission/emission_model_config.cpp
  src/forward_model/emission/emission_post_process.cpp
  src/forward_model/emission/emission_priors.cpp
  src/forward_model/emission/emission_init.cpp
  src/forward_model/emission/emission.cpp
  src/forward_model/emission/emission_test.cpp
  src/forward_model/secondary_eclipse/secondary_eclipse_model_config.cpp
  src/forward_model/secondary_eclipse/secondary_eclipse_post_process.cpp
  src/forward_model/secondary_eclipse/secondary_eclipse_priors.cpp
  src/forward_model/secondary_eclipse/secondary_eclipse_init.cpp
  src/forward_model/secondary_eclipse/secondary_eclipse.cpp
  src/forward_model/secondary_eclipse/secondary_eclipse_test.cpp
  src/forward_model/secondary_eclipse_bb/secondary_eclipse_bb_model_config.cpp
  src/forward_model/secondary_eclipse_bb/secondary_eclipse_bb_post_process.cpp
  src/forward_model/secondary_eclipse_bb/secondary_eclipse_bb_priors.cpp
  src/forward_model/secondary_eclipse_bb/secondary_eclipse_bb_init.cpp
  src/forward_model/secondary_eclipse_bb/secondary_eclipse_bb.cpp
  src/forward_model/secondary_eclipse_bb/secondary_eclipse_bb_test.cpp
  src/forward_model/transmission/transmission_model_config.cpp
  src/forward_model/transmission/transmission_post_process.cpp
  src/forward_model/transmission/transmission_priors.cpp
  src/forward_model/transmission/transmission_init.cpp
  src/forward_model/transmission/transmission.cpp
  src/forward_model/transmission/transmission_spectrum.cpp
  src/forward_model/transmission/transmission_test.cpp
  src/forward_model/flat_line/flat_line_priors.cpp
  src/forward_model/flat_line/flat_line.cpp
  src/forward_model/flat_line/flat_line_post_process.cpp
  src/forward_model/flat_line/flat_line_test.cpp
  src/forward_model/atmosphere/atmosphere.cpp
  src/forward_model/stellar_spectrum/sampled_stellar_spectrum.cpp
  src/forward_model/stellar_spectrum/spectrum_file.cpp
  src/forward_model/stellar_spectrum/star_blackbody.cpp
  src/forward_model/stellar_spectrum/star_file_spectrum.cpp
  src/forward_model/stellar_spectrum/stellar_spectrum_grid.cpp
  src/forward_model/modules/stellar_contamination/stellar_contamination.cpp
  src/forward_model/modules/velocity_broadening/velocity_broadening.cpp
  src/forward_model/modules/phase_resolved_broadening/phase_resolved_broadening.cpp
  src/observations/highres_observation.cpp
  src/observations/highres_observation_load_data.cpp
  src/cloud_model/grey_cloud_model.cpp
  src/cloud_model/fixed_cloud_model.cpp
  src/cloud_model/kh_cloud_model.cpp
  src/cloud_model/power_law_cloud_model.cpp
  src/model_main/model_main.cpp
  src/observations/observations.cpp
  src/observations/observations_load_data.cpp
  src/observations/observations_filter_response.cpp
  src/radiative_transfer/disort/discrete_ordinate.cpp
  src/radiative_transfer/adding_doubling/adding_doubling.cpp
  src/radiative_transfer/short_characteristics/short_characteristics.cpp
  src/retrieval/multinest_parameter.cpp
  src/retrieval/multinest_loglike.cpp
  src/retrieval/post_process_spectra.cpp
  src/retrieval/post_process.cpp
  src/retrieval/select_forward_model.cpp
  src/retrieval/retrieval_load_observations.cpp
  src/retrieval/retrieval.cpp
  src/retrieval/priors.cpp
  src/spectral_grid/convert.cpp
  src/spectral_grid/spectral_band_convolve.cpp
  src/spectral_grid/spectral_band_integrate.cpp
  src/spectral_grid/spectral_band.cpp
  src/spectral_grid/spectral_grid.cpp
  src/spectral_grid/spectral_grid_interpolate.cpp
  src/transport_coeff/cross_section_file.cpp
  src/transport_coeff/opacity_species.cpp
  src/transport_coeff/sampled_data.cpp
  src/transport_coeff/species_rayleigh_cross_sections.cpp
  src/transport_coeff/transport_coeff.cpp
  src/transport_coeff/h_m.cpp
  src/chemistry/isoprofile_chemistry.cpp
  src/chemistry/background_chemistry.cpp
  src/chemistry/isoprofile_clr_chemistry.cpp
  src/chemistry/fastchem_chemistry.cpp
  src/chemistry/free_chemistry.cpp
  src/chemistry/free_cbspline_chemistry.cpp
  src/chemistry/fixed_chemistry.cpp
  src/chemistry/step_function_chemistry.cpp
  src/temperature/piecewise_poly_temperature.cpp
  src/temperature/cubic_b_spline_temperature.cpp
  src/temperature/milne_solution_temperature.cpp
  src/temperature/guillot_temperature.cpp
  src/temperature/adiabate_cubic_spline.cpp
)

set(SRC_CUDA
  src/CUDA_kernels/band_integration_kernels.cu
  src/CUDA_kernels/convolution_kernels.cu
  src/CUDA_kernels/filter_response_kernels.cu
  src/CUDA_kernels/cross_section_kernels.cu
  src/CUDA_kernels/data_management_kernels.cu
  src/CUDA_kernels/log_like_kernels.cu
  src/CUDA_kernels/secondary_eclipse_kernels.cu
  src/CUDA_kernels/secondary_eclipse_bb.cu
  src/CUDA_kernels/transmission_spectrum_kernels.cu
  src/radiative_transfer/short_characteristics/short_characteristics_kernels.cu
  src/radiative_transfer/adding_doubling/adding_doubling_kernels.cu
  src/CUDA_kernels/rayleigh_scattering_kernels.cu
  src/CUDA_kernels/cloud_model_kernels.cu
  src/CUDA_kernels/contribution_function_kernels.cu
  src/CUDA_kernels/star_blackbody.cu
  src/CUDA_kernels/star_file_spectrum.cu
  src/CUDA_kernels/stellar_contamination.cu
  src/CUDA_kernels/stellar_spectrum_grid.cu
  src/CUDA_kernels/observation_kernels.cu
  src/CUDA_kernels/radiative_transfer_kernels.cu
  src/CUDA_kernels/batched_opacity_kernels.cu
  src/CUDA_kernels/doppler_shift_kernels.cu
  src/CUDA_kernels/highres_loglike_kernels.cu
  src/forward_model/modules/velocity_broadening/convolution_kernels_highres.cu
  src/forward_model/modules/phase_resolved_broadening/convolution_kernels_phase_resolved.cu)


set(SRC_PYTHON_WRAPPER
  python/src/bear_python_wrapper.cpp)


#check for OpenMP
find_package(OpenMP REQUIRED C CXX)



include(FetchContent)

#download the Boost Math library if necessary
FetchContent_Declare(
  boost_math
  GIT_REPOSITORY https://github.com/boostorg/math
  GIT_TAG        ed01dae24893bb69c02c6d599acc74bdb8f46bda
)

FetchContent_GetProperties(boost_math)

if (NOT boost_math_POPULATED)
  FetchContent_Populate(boost_math)
  add_subdirectory(${boost_math_SOURCE_DIR} ${boost_math_BINARY_DIR})
endif()



#download FastChem 3, if necessary
FetchContent_Declare(
  fastchem
  GIT_REPOSITORY https://github.com/newstrangeworlds/fastchem
  GIT_TAG        0979ebe637de51241c77e73071b0a3da5a16c8f2
)

FetchContent_GetProperties(fastchem)
set(FASTCHEM_ONLY_LIBRARY ON)

if (NOT fastchem_POPULATED)
  FetchContent_Populate(fastchem)
  add_subdirectory(${fastchem_SOURCE_DIR} ${fastchem_BINARY_DIR})
endif()



#download Eigen (shared dependency of DisORT and Adding-Doubling RT)
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        3.4.0
  GIT_SHALLOW    TRUE
)

set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

FetchContent_GetProperties(eigen)

if (NOT eigen_POPULATED)
  FetchContent_Populate(eigen)
  add_subdirectory(${eigen_SOURCE_DIR} ${eigen_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()



#download DisORT C++ library
FetchContent_Declare(
  disortpp
  GIT_REPOSITORY https://github.com/NewStrangeWorlds/DisORT
  GIT_TAG        main
  GIT_SHALLOW    TRUE
)

set(BUILD_PYTHON_BINDINGS OFF CACHE BOOL "" FORCE)
set(DISORTPP_MARCH_NATIVE ON CACHE BOOL "" FORCE)

FetchContent_GetProperties(disortpp)

if (NOT disortpp_POPULATED)
  FetchContent_Populate(disortpp)
  add_subdirectory(${disortpp_SOURCE_DIR} ${disortpp_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()



#download Adding-Doubling RT library
FetchContent_Declare(
  adding_doubling
  GIT_REPOSITORY https://github.com/NewStrangeWorlds/AddingDoublingRT
  GIT_TAG        main
  GIT_SHALLOW    TRUE
)

set(ADRT_ENABLE_CUDA ON CACHE BOOL "" FORCE)

FetchContent_GetProperties(adding_doubling)

if (NOT adding_doubling_POPULATED)
  FetchContent_Populate(adding_doubling)
  add_subdirectory(${adding_doubling_SOURCE_DIR} ${adding_doubling_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()



#download Multibnest if necessary
FetchContent_Declare(
  multinest
  GIT_REPOSITORY https://github.com/farhanferoz/MultiNest
  GIT_TAG        cc616c17ed440fbff8d35937997c86256ed4d744
)

FetchContent_GetProperties(multinest)

if (NOT multinest_POPULATED)
  FetchContent_Populate(multinest)
endif()


set(MultiNest_SOURCE
  ${multinest_SOURCE_DIR}/MultiNest_v3.12_CMake/multinest/src/cwrapper.f90
  ${multinest_SOURCE_DIR}/MultiNest_v3.12_CMake/multinest/src/kmeans_clstr.f90
  ${multinest_SOURCE_DIR}/MultiNest_v3.12_CMake/multinest/src/nested.F90
  ${multinest_SOURCE_DIR}/MultiNest_v3.12_CMake/multinest/src/priors.f90
  ${multinest_SOURCE_DIR}/MultiNest_v3.12_CMake/multinest/src/posterior.F90
  ${multinest_SOURCE_DIR}/MultiNest_v3.12_CMake/multinest/src/utils.f90
  ${multinest_SOURCE_DIR}/MultiNest_v3.12_CMake/multinest/src/utils1.f90
  ${multinest_SOURCE_DIR}/MultiNest_v3.12_CMake/multinest/src/xmeans_clstr.f90)
  

# If we are compiling with gfortran, add -ffree-line-length-none
foreach(lang C CXX Fortran)
    if(CMAKE_${lang}_COMPILER_ID STREQUAL GNU)
        MESSAGE(STATUS "Detected gfortran, adding -ffree-line-length-none compiler flag.")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
	 if(CMAKE_${lang}_COMPILER_VERSION VERSION_GREATER 10)
            set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -std=legacy")
	    endif()
	break()
    endif()
endforeach()

# libmultinest.so
add_library(multinest_shared SHARED ${MultiNest_SOURCE})
SET_TARGET_PROPERTIES(multinest_shared PROPERTIES LINKER_LANGUAGE Fortran)


#check for LAPACK library (required by MultiNest)
find_package(LAPACK REQUIRED)



#compilation target for the CUDA files
#CUDA requires different compiler options, so it's compiled separately
add_library(bear_cuda ${SRC_CUDA})
target_include_directories(bear_cuda PUBLIC ${boost_math_SOURCE_DIR}/include/)
target_link_libraries(bear_cuda PUBLIC adding_doubling_cuda)

if(DEBUG_CUDA MATCHES ON)
  target_compile_definitions(bear_cuda PRIVATE DEBUG_CUDA)
endif()


if (USE_PYTHON MATCHES ON)
    #find the python version
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    
    set(PYBIND11_FINDPYTHON ON)

    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG        7c33cdc2d39c7b99a122579f53bc94c8eb3332ff
    )

    FetchContent_GetProperties(pybind11)
    if(NOT pybind11_POPULATED)
        FetchContent_Populate(pybind11)
        add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
    endif()
    
  

  #compilation target for BeAR's Python library
  add_library(bear_for_py ${SRC})
  target_include_directories(bear_for_py PUBLIC ${boost_math_SOURCE_DIR}/include/)
  target_include_directories(bear_for_py PUBLIC ${eigen_SOURCE_DIR})
  
  add_library(bear_cuda_for_py ${SRC_CUDA})
  target_include_directories(bear_cuda_for_py PUBLIC ${boost_math_SOURCE_DIR}/include/)
  target_link_libraries(bear_cuda_for_py PUBLIC adding_doubling_cuda)

  if(DEBUG_CUDA MATCHES ON)
    target_compile_definitions(bear_cuda_for_py PRIVATE DEBUG_CUDA)
  endif()
  
  
  target_compile_options(bear_for_py PRIVATE "${OpenMP_CXX_FLAGS}" PRIVATE "-O3" PRIVATE "-fPIC")
  target_compile_options(bear_cuda_for_py PRIVATE "-fPIC")

  
  #compilation target for BeAR's Python wrapper
  pybind11_add_module(pybear ${SRC_PYTHON_WRAPPER})
  
  set_target_properties(pybear PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/lib")
  
  target_link_libraries(pybear PUBLIC fastchem_lib_for_py bear_cuda_for_py bear_for_py multinest_shared disortpp adding_doubling OpenMP::OpenMP_CXX pybind11::module ${LAPACK_LIBRARIES})
endif()






# build application for BeAR
add_executable(bear ${SRC})
target_include_directories(bear PUBLIC ${boost_math_SOURCE_DIR}/include/)
target_include_directories(bear PUBLIC ${eigen_SOURCE_DIR})


#link files for the final BeAR code
target_link_libraries(bear PUBLIC fastchem_lib bear_cuda multinest_shared disortpp adding_doubling OpenMP::OpenMP_CXX ${LAPACK_LIBRARIES})
