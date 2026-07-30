

#include "contribution_function_kernels.h"


#include <iostream>
#include <vector>
#include "math.h"
#include <stdio.h>
#include <new>

#include "data_management_kernels.h"
#include "../additional/physical_const.h"


#include "error_check.h"
#include "planck_function.h"

namespace bear{



__global__ 
void contributionFunctionDevice(
  float* contribution_function_gpu,
  const float* absorption_coeff_dev, 
  const double* wavenumber_list_dev,
  const double* temperature_dev, 
  const double* vertical_grid_dev,
  const int nb_spectral_points, 
  const int nb_grid_points)
{
  //tid is the wavenumber index
  //int tid = blockIdx.x * blockDim.x + threadIdx.x;
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {

    double cumulative_transmission = 1.0;

    //the absorption coefficient of the lower level of each layer is carried
    //over to the next iteration to avoid reading it from global memory twice
    double abs_upper = absorption_coeff_dev[(nb_grid_points-1)*nb_spectral_points + tid];

    for (int i=nb_grid_points-1; i>0; i--)
    {
      const double abs_lower = absorption_coeff_dev[(i-1)*nb_spectral_points + tid];

      const double delta_z = (vertical_grid_dev[i] - vertical_grid_dev[i-1]);
      const double optical_depth_layer = delta_z * (abs_upper + abs_lower)/2.;

      const double layer_transmission = exp(-optical_depth_layer);

      contribution_function_gpu[i*nb_spectral_points + tid] = static_cast<float>(2 * constants::pi * planckFunction(temperature_dev[i], wavenumber_list_dev[tid]) * (1.0 - layer_transmission) * cumulative_transmission);

      cumulative_transmission *= layer_transmission;

      abs_upper = abs_lower;
    }
  }


}






__host__ void contributionFunctionGPU(
  float* contribution_function_dev,
  float* absorption_coeff_dev, 
  double* wavenumber_list_dev,
  std::vector<double>& temperature, 
  std::vector<double>& vertical_grid,
  const size_t& nb_spectral_points)
{
  size_t nb_grid_points = temperature.size();

  const int bytes = nb_grid_points*sizeof(double);



  double* temperature_dev = nullptr;

  cudaMalloc(&temperature_dev, bytes);
  cudaMemcpy(temperature_dev, &temperature[0], bytes, cudaMemcpyHostToDevice);


  double* vertical_grid_dev = nullptr;

  cudaMalloc(&vertical_grid_dev, bytes);
  cudaMemcpy(vertical_grid_dev, &vertical_grid[0], bytes, cudaMemcpyHostToDevice);


  int threads = 256;
  //int blocks = min(( int(nb_spectral_points)+ threads-1)/threads, 2048);
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;



  contributionFunctionDevice<<<blocks,threads>>>(contribution_function_dev,
                                                 absorption_coeff_dev, wavenumber_list_dev,
                                                 temperature_dev, vertical_grid_dev,
                                                 nb_spectral_points, nb_grid_points);


  CUDA_CHECK_AFTER_KERNEL();

  gpuErrchk(cudaFree(temperature_dev));
  gpuErrchk(cudaFree(vertical_grid_dev));
}






}
