
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

cmake_minimum_required (VERSION 3.10 FATAL_ERROR)

# check command line args, user can specify SM to override manual setting
set(SM "0" CACHE STRING "GPU SM value")


project(helios-r C CXX CUDA Fortran)


#some C++/C flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -pedantic -MMD")
set(CMAKE_C_FLAGS "-Wall -pedantic -MMD")


set(CMAKE_BUILD_TYPE Release)



if(SM MATCHES "0")
  set(CUDA_ARCH_FLAGS "-arch=sm_61")
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
string(APPEND CMAKE_CUDA_FLAGS " -std=c++11 -lineinfo -Xptxas -v")

#string (APPEND CMAKE_CUDA_FLAGS " -cudart shared" )



#define the source files 
set(SRC
  helios_src/additional/aux_functions.cpp
  helios_src/additional/piecewise_poly.cpp
  helios_src/config/global_config.cpp
  helios_src/forward_model/emission/emission_model_config.cpp
  helios_src/forward_model/emission/emission_post_process_temperatures.cpp
  helios_src/forward_model/emission/emission_post_process.cpp
  helios_src/forward_model/emission/emission_priors.cpp
  helios_src/forward_model/emission/emission_init.cpp
  helios_src/forward_model/emission/emission.cpp
  helios_src/forward_model/emission/emission_test.cpp
  helios_src/forward_model/secondary_eclipse/secondary_eclipse_model_config.cpp
  helios_src/forward_model/secondary_eclipse/secondary_eclipse_post_process_temperatures.cpp
  helios_src/forward_model/secondary_eclipse/secondary_eclipse_post_process.cpp
  helios_src/forward_model/secondary_eclipse/secondary_eclipse_priors.cpp
  helios_src/forward_model/secondary_eclipse/secondary_eclipse_init.cpp
  helios_src/forward_model/secondary_eclipse/secondary_eclipse.cpp
  helios_src/forward_model/secondary_eclipse/secondary_eclipse_test.cpp
  helios_src/forward_model/transmission/transmission_model_config.cpp
  helios_src/forward_model/transmission/transmission_post_process.cpp
  helios_src/forward_model/transmission/transmission_priors.cpp
  helios_src/forward_model/transmission/transmission_init.cpp
  helios_src/forward_model/transmission/transmission.cpp
  helios_src/forward_model/transmission/transmission_spectrum.cpp
  helios_src/forward_model/transmission/transmission_test.cpp
  helios_src/forward_model/atmosphere/atmosphere.cpp
  helios_src/cloud_model/grey_cloud_model.cpp
  helios_src/cloud_model/kh_cloud_model.cpp
  helios_src/model_main/model_main.cpp
  helios_src/observations/observations.cpp
  helios_src/observations/observations_load_data.cpp
  helios_src/observations/observations_filter_response.cpp
  helios_src/radiative_transfer/discrete_ordinate.cpp
  helios_src/radiative_transfer/short_characteristics.cpp
  helios_src/retrieval/multinest_parameter.cpp
  helios_src/retrieval/multinest_loglike.cpp
  helios_src/retrieval/post_process_spectra.cpp
  helios_src/retrieval/post_process.cpp
  helios_src/retrieval/select_forward_model.cpp
  helios_src/retrieval/retrieval_load_observations.cpp
  helios_src/retrieval/retrieval.cpp
  helios_src/retrieval/priors.cpp
  helios_src/spectral_grid/convert.cpp
  helios_src/spectral_grid/spectral_band_convolve.cpp
  helios_src/spectral_grid/spectral_band_integrate.cpp
  helios_src/spectral_grid/spectral_band.cpp
  helios_src/spectral_grid/spectral_grid.cpp
  helios_src/spectral_grid/spectral_grid_interpolate.cpp
  helios_src/transport_coeff/cross_section_file.cpp
  helios_src/transport_coeff/opacity_species.cpp
  helios_src/transport_coeff/sampled_data.cpp
  helios_src/transport_coeff/species_rayleigh_cross_sections.cpp
  helios_src/transport_coeff/transport_coeff.cpp
  helios_src/chemistry/isoprofile_chemistry.cpp
  helios_src/chemistry/fastchem_chemistry.cpp
  helios_src/chemistry/free_chemistry.cpp
  helios_src/temperature/piecewise_poly_temperature.cpp
  helios_src/radiative_transfer/cdisort_src/cdisort.c
  helios_src/radiative_transfer/cdisort_src/locate.c)
  

set(SRC_CUDA
  helios_src/CUDA_kernels/band_integration_kernels.cu
  helios_src/CUDA_kernels/convolution_kernels.cu
  helios_src/CUDA_kernels/filter_response_kernels.cu
  helios_src/CUDA_kernels/cross_section_kernels.cu
  helios_src/CUDA_kernels/data_management_kernels.cu
  helios_src/CUDA_kernels/log_like_kernels.cu
  helios_src/CUDA_kernels/secondary_eclipse_kernels.cu
  helios_src/CUDA_kernels/transmission_spectrum_kernels.cu
  helios_src/CUDA_kernels/short_characteristics_kernels.cu
  helios_src/CUDA_kernels/rayleigh_scattering_kernels.cu
  helios_src/CUDA_kernels/cloud_model_kernels.cu
  helios_src/CUDA_kernels/contribution_function_kernels.cu)


#check for OpenMP
find_package(OpenMP REQUIRED)



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



#download the FastChem 2.1 if necessary
FetchContent_Declare(
  fastchem
  GIT_REPOSITORY https://github.com/exoclime/fastchem
  GIT_TAG        e8dbeb7ec7719dbd0021d672a0b5d330e924a764
)

FetchContent_GetProperties(fastchem)
set(FASTCHEM_ONLY_LIBRARY ON)

if (NOT fastchem_POPULATED)
  FetchContent_Populate(fastchem)
  add_subdirectory(${fastchem_SOURCE_DIR} ${fastchem_BINARY_DIR})
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
add_library(helios_r_cuda ${SRC_CUDA})
target_include_directories(helios_r_cuda PUBLIC ${boost_math_SOURCE_DIR}/include/)


# build application for Helios-r
add_executable(helios-r ${SRC})
target_include_directories(helios-r PUBLIC ${boost_math_SOURCE_DIR}/include/)


#link files for the final helios-r code
target_link_libraries(helios-r PUBLIC fastchem_lib helios_r_cuda multinest_shared OpenMP::OpenMP_CXX ${LAPACK_LIBRARIES})
