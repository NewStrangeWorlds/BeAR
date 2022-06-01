/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
*
* Helios-r2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Helios-r2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* Helios-r2 directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#include "log_like_kernels.h"

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "error_check.h"
#include "reduce_kernels.h"
#include "../additional/physical_const.h"


namespace helios{


//computes the log likelihood
//each thread calculates one data point, the results are then summed by each block via a blockReduce method and finally all collected by thread 0 
__global__ void logLikeDevice(double* observation, double* error, double* model, const int nb_spectral_points, const double error_inflation_coefficient, double* d_log_like)
{
  double d_log_like_sum = 0;


  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_spectral_points; i += blockDim.x * gridDim.x)
  {
    //Eq. 22 from Paper I
    const double error_square = error[i]*error[i] + error_inflation_coefficient;

    //Eq. 23 from Paper I
    d_log_like_sum += - 0.5 * log(error_square * 2.0 * constants::pi) - 0.5 * (observation[i] - model[i])*(observation[i] - model[i]) / error_square;
  }


  d_log_like_sum = blockReduceSum(d_log_like_sum);


  if (threadIdx.x == 0)
    atomicAdd(d_log_like, d_log_like_sum);

}





__host__ double logLikeHost(double* observation, double* observation_error, double* model_spectrum, const size_t nb_spectral_points, const double error_inflation_coefficient)
{
  double h_log_like = 0;

  double* d_log_like = nullptr;

  cudaMalloc(&d_log_like, sizeof(double));
  cudaMemset(d_log_like, 0, sizeof(double));

 

  int threads = 256;

  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;


  logLikeDevice<<<blocks,threads>>>(observation,observation_error,model_spectrum, nb_spectral_points, error_inflation_coefficient, d_log_like);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );


  cudaMemcpy(&h_log_like, d_log_like, sizeof(double), cudaMemcpyDeviceToHost);

  cudaDeviceSynchronize();
  
  cudaFree(d_log_like);
  
  cudaDeviceSynchronize();


  return h_log_like;
}


}
