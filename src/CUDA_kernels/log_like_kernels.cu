/*
* This file is part of the BeAR code (https://github.com/newstrangeworlds/BeAR).
* Copyright (C) 2024 Daniel Kitzmann
*
* BeAR is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* BeAR is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* BeAR directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <cmath>
#include <stdio.h>

#include "../retrieval/retrieval.h"

#include "error_check.h"
#include "reduce_kernels.h"
#include "../additional/physical_const.h"


namespace bear{


//computes the log likelihood
//each thread calculates one data point, the results are then summed by each block via a blockReduce method and finally all collected by thread 0 
__global__ void logLikeDevice(
  double* observation,
  double* error,
  double* loglike_weight,
  double* model,
  const int nb_spectral_points,
  const double error_inflation_coefficient,
  double* d_log_like)
{
  double d_log_like_sum = 0;


  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_spectral_points; i += blockDim.x * gridDim.x)
  {
    //Eq. 22 from Paper I
    const double error_square = error[i]*error[i] + error_inflation_coefficient;

    //Eq. 23 from Paper I
    d_log_like_sum += (- 0.5 * log(error_square * 2.0 * constants::pi) - 0.5 * (observation[i] - model[i])*(observation[i] - model[i]) / error_square) * loglike_weight[i];
  }


  d_log_like_sum = blockReduceSum(d_log_like_sum);


  if (threadIdx.x == 0)
    atomicAdd(d_log_like, d_log_like_sum);
}



__host__ double Retrieval::logLikeDev(
  std::vector<double*> model_spectrum,
  const double error_inflation_coefficient)
{
  double log_like = 0;

  for (size_t i=0; i<nb_observations; ++i)
  {
    double* d_log_like = nullptr;

    cudaMalloc(&d_log_like, sizeof(double));
    cudaMemset(d_log_like, 0, sizeof(double));

    const int threads = 128;
    const int nb_points = observations[i].nbPoints();

    int blocks = nb_points / threads;
    if (nb_points % threads) blocks++;


    logLikeDevice<<<blocks,threads>>>(
      observations[i].data_gpu,
      observations[i].data_error_gpu,
      observations[i].likelihood_weight_gpu,
      model_spectrum[i],
      nb_points,
      error_inflation_coefficient,
      d_log_like);


    cudaDeviceSynchronize();
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    double h_log_like = 0;
    cudaMemcpy(&h_log_like, d_log_like, sizeof(double), cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();
    cudaFree(d_log_like);
    cudaDeviceSynchronize();

    log_like += h_log_like;
  }

  return log_like;
}


}
