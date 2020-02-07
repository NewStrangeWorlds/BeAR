
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


project(helios-r C CXX CUDA)


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
  helios_src/forward_model/brown_dwarf_model_config.cpp
  helios_src/forward_model/brown_dwarf_post_process_temperatures.cpp
  helios_src/forward_model/brown_dwarf_post_process.cpp
  helios_src/forward_model/brown_dwarf_priors.cpp
  helios_src/forward_model/brown_dwarf_init.cpp
  helios_src/forward_model/brown_dwarf.cpp
  helios_src/model_main/model_main.cpp
  helios_src/observations/observations.cpp
  helios_src/radiative_transfer/discrete_ordinate.cpp
  helios_src/radiative_transfer/short_characteristics.cpp
  helios_src/retrieval/multinest_parameter.cpp
  helios_src/retrieval/multinest_loglike.cpp
  helios_src/retrieval/post_process_spectra.cpp
  helios_src/retrieval/post_process.cpp
  helios_src/retrieval/retrieval_set_priors.cpp
  helios_src/retrieval/retrieval_load_observations.cpp
  helios_src/retrieval/retrieval.cpp
  helios_src/spectral_grid/convert.cpp
  helios_src/spectral_grid/spectral_band_convolve.cpp
  helios_src/spectral_grid/spectral_band_integrate.cpp
  helios_src/spectral_grid/spectral_band.cpp
  helios_src/spectral_grid/spectral_grid.cpp
  helios_src/transport_coeff/cia_coefficients.cpp
  helios_src/transport_coeff/h2.cpp
  helios_src/transport_coeff/sampled_cross_section_data.cpp
  helios_src/transport_coeff/transport_coeff_single_species.cpp
  helios_src/transport_coeff/transport_coeff.cpp
  helios_src/chemistry/isoprofile_chemistry.cpp
  helios_src/chemistry/fastchem_chemistry.cpp
  helios_src/temperature/piecewise_poly_temperature.cpp
  cdisort_src/cdisort.c
  cdisort_src/locate.c)

  

set(SRC_FASTCHEM
  fastchem_src/bisection.cpp
  fastchem_src/calc_densities_ph.cpp
  fastchem_src/calc_densities.cpp
  fastchem_src/calc_mean_mol_weight.cpp
  fastchem_src/calc_total_h_density.cpp
  fastchem_src/check.cpp
  fastchem_src/fastchem.cpp
  fastchem_src/get.cpp
  fastchem_src/init_add_species.cpp
  fastchem_src/init_read_files.cpp
  fastchem_src/init_solver.cpp
  fastchem_src/init.cpp
  fastchem_src/mass_action_constant.cpp
  fastchem_src/set.cpp
  fastchem_src/solve_fastchem.cpp
  fastchem_src/solver_linsol_quadsol.cpp
  fastchem_src/solver_nelder_mead_single.cpp
  fastchem_src/solver_newtsol.cpp
  fastchem_src/solver_scaling_factor.cpp
  fastchem_src/solver.cpp)
  

set(SRC_CUDA
  helios_src/CUDA_kernels/band_integration_kernels.cu
  helios_src/CUDA_kernels/convolution_kernels.cu
  helios_src/CUDA_kernels/cross_section_kernels.cu
  helios_src/CUDA_kernels/data_management_kernels.cu
  helios_src/CUDA_kernels/log_like_kernels.cu
  helios_src/CUDA_kernels/short_characteristics_kernels.cu)



#compilation target for the CUDA files
#CUDA requires different compiler options, so it's compiled separately
add_library(helios_r_cuda ${SRC_CUDA})


#compilation target for FastChem
add_library(fastchem ${SRC_FASTCHEM})
target_compile_options(fastchem PRIVATE "-Ofast")



# build application for Helios-r
add_executable(helios-r ${SRC})
  
  
#link the Boost library 
#note that Helios-r only needs the Boost header files
find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()


#check for OpenMP
find_package(OpenMP REQUIRED)


#check for LAPACK library (required by MultiNest)
find_package(LAPACK REQUIRED)


#link gfortran (MultiNest is written in Fortran)
find_library(GFORTRAN_LIBRARY NAMES gfortran)


#add the static MultiNest library
set(MULTINEST_LOCATION ${CMAKE_SOURCE_DIR})
string(APPEND MULTINEST_LOCATION "/multinest/libmultinest.a")
add_library(multinest STATIC IMPORTED)
set_property(TARGET multinest PROPERTY IMPORTED_LOCATION ${MULTINEST_LOCATION})



#link files for the final helios-r code
target_link_libraries(helios-r PUBLIC fastchem helios_r_cuda multinest OpenMP::OpenMP_CXX ${GFORTRAN_LIBRARY} ${LAPACK_LIBRARIES})

