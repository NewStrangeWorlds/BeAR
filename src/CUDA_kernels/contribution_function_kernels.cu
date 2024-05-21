

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

namespace helios{



__global__ void contributionFunctionDevice(double* contribution_function_gpu,
                                           const double* absorption_coeff_dev, const double* wavenumber_list_dev,
                                           const double* temperature_dev, const double* vertical_grid_dev,
                                           const int nb_spectral_points, const int nb_grid_points)
{


  //tid is the wavenumber index
  //int tid = blockIdx.x * blockDim.x + threadIdx.x;
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {

    double cumulative_optical_depth = 0;
    double cumulative_transmission = 1.0;


    for (int i=nb_grid_points-1; i>0; i--)
    {
      const double delta_z = (vertical_grid_dev[i] - vertical_grid_dev[i-1]);
      const double optical_depth_layer = delta_z * ( absorption_coeff_dev[i*nb_spectral_points + tid ] + absorption_coeff_dev[(i-1)*nb_spectral_points + tid])/2.;

      const double layer_transmission = exp(-optical_depth_layer);

      contribution_function_gpu[i*nb_spectral_points + tid] = 2 * constants::pi * planckFunction(temperature_dev[i], wavenumber_list_dev[tid]) * (1.0 - layer_transmission) * cumulative_transmission;

      cumulative_transmission *= layer_transmission;
      //printf("%d  %d  %f  %f  %f  %f\n", tid, i, optical_depth_layer, layer_transmission, cumulative_transmission, contribution_function_gpu[i*nb_spectral_points + tid]);

      /*const double delta_z = (vertical_grid_dev[i] - vertical_grid_dev[i-1]);

      const double optical_depth_layer = delta_z * ( absorption_coeff_dev[i*nb_spectral_points + tid ] + absorption_coeff_dev[(i-1)*nb_spectral_points + tid])/2.;

      const double cumulative_optical_depth_layer = cumulative_optical_depth + optical_depth_layer;

      const double weighting_function = (-exp(-cumulative_optical_depth_layer) + exp(-cumulative_optical_depth)) / delta_z;

      cumulative_optical_depth = cumulative_optical_depth_layer;
      
      contribution_function_gpu[i*nb_spectral_points + tid] = planckFunction(temperature_dev[i], wavenumber_list_dev[tid]) * weighting_function;*/
    }
  }


}






__host__ void contributionFunctionGPU(double* contribution_function_dev,
                                      double* absorption_coeff_dev, double* wavenumber_list_dev,
                                      std::vector<double>& temperature, std::vector<double>& vertical_grid,
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


  cudaThreadSynchronize();


  int threads = 256;
  //int blocks = min(( int(nb_spectral_points)+ threads-1)/threads, 2048);
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;



  contributionFunctionDevice<<<blocks,threads>>>(contribution_function_dev,
                                                 absorption_coeff_dev, wavenumber_list_dev,
                                                 temperature_dev, vertical_grid_dev,
                                                 nb_spectral_points, nb_grid_points);


  cudaThreadSynchronize();

  cudaFree(temperature_dev);
  cudaFree(vertical_grid_dev);

  gpuErrchk( cudaPeekAtLastError() ); 
}






}
