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


#include "cross_section_kernels.h"

#include <iostream>
#include <vector>
#include "math.h"
#include <stdio.h>
#include <new>

#include <fstream>
#include <string>
#include <iomanip>

#include "error_check.h"
#include "data_management_kernels.h"
#include <cstring>
#include "../additional/physical_const.h"
#include "../transport_coeff/species_definition.h"
#include "../transport_coeff/opacity_species.h"


namespace bear{


// H- bound-free coefficients (constant memory)
__constant__ double d_C_n[7] = {0.0, 152.519, 49.534, -118.858, 92.536, -34.194, 4.982};

// H- free-free coefficients for lambda > 0.3645 micron (constant memory)
__constant__ double d_A_n1[7] = {0.0, 0.0, 2483.3460, -3449.8890, 2200.0400, -696.2710, 88.2830};
__constant__ double d_B_n1[7] = {0.0, 0.0, 285.8270, -1158.3820, 2427.7190, -1841.4000, 444.5170};
__constant__ double d_C_n1[7] = {0.0, 0.0, -2054.2910, 8746.5230, -13651.1050, 8624.9700, -1863.8650};
__constant__ double d_D_n1[7] = {0.0, 0.0, 2827.7760, -11485.6320, 16755.5240, -10051.5300, 2095.2880};
__constant__ double d_E_n1[7] = {0.0, 0.0, -1341.5370, 5303.6090, -7510.4940, 4400.0670, -901.7880};
__constant__ double d_F_n1[7] = {0.0, 0.0, 208.9520, -812.9390, 1132.7380, -655.0200, 132.9850};

// H- free-free coefficients for 0.1823 <= lambda <= 0.3645 micron (constant memory)
__constant__ double d_A_n2[7] = {0.0, 518.1021, 473.2636, -482.2089, 115.5291, 0.0, 0.0};
__constant__ double d_B_n2[7] = {0.0, -734.8666, 1443.4137, -737.1616, 169.6374, 0.0, 0.0};
__constant__ double d_C_n2[7] = {0.0, 1021.1775, -1977.3395, 1096.8827, -245.6490, 0.0, 0.0};
__constant__ double d_D_n2[7] = {0.0, -479.0721, 922.3575, -521.1341, 114.2430, 0.0, 0.0};
__constant__ double d_E_n2[7] = {0.0, 93.1373, -178.9275, 101.7963, -21.9972, 0.0, 0.0};
__constant__ double d_F_n2[7] = {0.0, -6.4285, 12.3600, -7.0571, 1.5097, 0.0, 0.0};


//tabulate the H- bound-free cross-section based on the formula from John (1988) 
//and the coefficients from Bell & Berrington (1987)
__global__ void calcHmbfCrossSectionDevice(
  const int nb_spectral_points,
  double* wavelengths_dev,
  float* __restrict__ cross_section_dev)
{
  const double lambda_0 = 1.6419; //photo-detachment threshold

  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < nb_spectral_points)
  {
    const double lambda = wavelengths_dev[tid];

    //outside the photo-detachment range the cross-section is zero;
    //without this, the caller would read uninitialised device memory
    cross_section_dev[tid] = 0.0f;

    if (lambda <= lambda_0 && lambda >= 0.125)
    {
      const double x = 1.0/lambda - 1.0/lambda_0;
      const double sqrt_x = sqrt(x);

      //powers of sqrt(x): x^0, x^0.5, x^1, x^1.5, x^2, x^2.5
      //computed incrementally instead of calling pow() 6 times
      //f = C_n[1]*x^0 + C_n[2]*x^0.5 + C_n[3]*x^1 + C_n[4]*x^1.5 + C_n[5]*x^2 + C_n[6]*x^2.5
      double x_power = 1.0;  // x^0
      double f = d_C_n[1];   // i=1: C_n[1] * x^0

      x_power *= sqrt_x;     // x^0.5
      f += d_C_n[2] * x_power;

      x_power *= sqrt_x;     // x^1
      f += d_C_n[3] * x_power;

      x_power *= sqrt_x;     // x^1.5
      f += d_C_n[4] * x_power;

      x_power *= sqrt_x;     // x^2
      f += d_C_n[5] * x_power;

      x_power *= sqrt_x;     // x^2.5
      f += d_C_n[6] * x_power;

      const double lambda3 = lambda * lambda * lambda;
      const double x15 = x * sqrt_x;  // x^1.5
      const double sigma = 1e-18 * lambda3 * x15 * f;

      cross_section_dev[tid] = sigma;
    }
  }
}



//fused H- continuum for all grid points in one launch:
//blockIdx.y indexes a grid-point work item; bound-free and free-free
//contributions are combined into a single write per (grid point, wavelength)
__global__ void calcHmContinuumBatchedDevice(
  const float* __restrict__ bf_cross_sections,
  const double* __restrict__ wavelengths_dev,
  const double* __restrict__ hm_number_densities,
  const double* __restrict__ h_e_products,        //h number density * e- pressure
  const double* __restrict__ temperatures,
  const int* __restrict__ grid_points,
  const int nb_spectral_points,
  float* __restrict__ absorption_coeff_device)
{
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int wi = blockIdx.y;

  if (tid >= nb_spectral_points) return;

  //bound-free contribution from the pre-tabulated cross-section
  double sigma = bf_cross_sections[tid] * hm_number_densities[wi];

  //free-free contribution (John 1988 / Bell & Berrington 1987)
  const double lambda = wavelengths_dev[tid];

  if (lambda >= 0.1823)
  {
    const double inv_l = 1.0 / lambda;
    const double l2 = lambda * lambda;
    const double inv_l2 = inv_l * inv_l;
    const double inv_l3 = inv_l2 * inv_l;
    const double inv_l4 = inv_l2 * inv_l2;

    const double theta = 5040.0 / temperatures[wi];
    const double sqrt_theta = sqrt(theta);
    double theta_power = theta;

    const bool use_set1 = (lambda > 0.3645);

    double ff_sigma = 0;

    #pragma unroll
    for (unsigned int i=1; i<7; ++i)
    {
      const double A_i = use_set1 ? d_A_n1[i] : d_A_n2[i];
      const double B_i = use_set1 ? d_B_n1[i] : d_B_n2[i];
      const double C_i = use_set1 ? d_C_n1[i] : d_C_n2[i];
      const double D_i = use_set1 ? d_D_n1[i] : d_D_n2[i];
      const double E_i = use_set1 ? d_E_n1[i] : d_E_n2[i];
      const double F_i = use_set1 ? d_F_n1[i] : d_F_n2[i];

      const double lambda_poly = l2 * A_i + B_i + C_i * inv_l + D_i * inv_l2
                                 + E_i * inv_l3 + F_i * inv_l4;

      ff_sigma += theta_power * lambda_poly;
      theta_power *= sqrt_theta;
    }

    ff_sigma *= 1e-29;
    if (ff_sigma < 0) ff_sigma = 0;

    sigma += ff_sigma * h_e_products[wi];
  }

  //each grid point appears at most once in the work list,
  //so this write needs no atomics
  absorption_coeff_device[grid_points[wi]*nb_spectral_points + tid] += sigma;
}


//set all absorption coefficients to 0
__global__ void initCrossSectionsDevice(
  const int nb_points, 
  float* absorption_coeff_device)
{
  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < nb_points)
    absorption_coeff_device[tid] = 0.0;
}



//prints the calculated cross-sections
//this is just for debug purposes
__global__ void printAC(const size_t nb_spectral_points, const size_t nb_grid_points, const size_t grid_point,
                        float* absorption_coeff_device)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid == 0)
  {
    for (int i=0; i<nb_grid_points; ++i)
    {
       //double* absorption_coeff = absorption_coeff_device[i]; //pointer to the absorption coefficient per wavenumber
       printf("%d %1.6e\n",i,absorption_coeff_device[i*nb_spectral_points + 1]);
    }

  }

}


//checks the calculated cross-sections for too large or small values
//this is just for debug purposes
__global__ 
void checkCrossSections(
  const int nb_points, 
  float* absorption_coeff_dev)
{
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_points; tid += blockDim.x * gridDim.x)
  {
    if (absorption_coeff_dev[tid] != absorption_coeff_dev[tid] || absorption_coeff_dev[tid] > 1e10 || absorption_coeff_dev[tid] < 0)
      printf("abs coeff are wrong: %d %1.6e\n",tid,absorption_coeff_dev[tid]);
  }

}



//sets all cross-sections back to zero
__host__ void initCrossSectionsHost(
  const size_t nb_points, 
  float* absorption_coeff_device)
{
  cudaMemset(absorption_coeff_device, 0, nb_points * sizeof(float));

  gpuErrchk( cudaPeekAtLastError() );
}



//Checks computed absorption coefficients for unrealistic values
//This is just a debug function that is normally only called to validate calculations
__host__ void checkCrossSectionsHost(
  const size_t nb_spectral_points, 
  const size_t nb_grid_points,
  float* absorption_coeff_device)
{
  int threads = 256;
  //int blocks = min((int (nb_spectral_points)/4+threads-1)/threads, 2048);
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  checkCrossSections<<<blocks,threads>>>(nb_spectral_points*nb_grid_points, absorption_coeff_device);

  CUDA_CHECK_AFTER_KERNEL();
}



__host__ void GasHm::calcContinuumBatchedGPU(
  const std::vector<double>& temperature,
  const std::vector<std::vector<double>>& number_densities,
  float* absorption_coeff_device)
{
  const int nb_spectral_points = spectral_grid->nbSpectralPoints();

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  //pre-tabulate the bound-free cross-section if not already done
  if (bound_free_sigma_dev == nullptr)
  {
    allocateOnDevice(bound_free_sigma_dev, nb_spectral_points);

    calcHmbfCrossSectionDevice<<<blocks,threads>>>(
      nb_spectral_points,
      spectral_grid->wavelength_list_gpu,
      bound_free_sigma_dev);

    CUDA_CHECK_AFTER_KERNEL();
  }

  //collect the per-layer scalars, skipping layers without any H-/H/e- content
  const size_t nb_grid_points = temperature.size();

  std::vector<double> hm_nd, h_e, temps;
  std::vector<int> gps;
  hm_nd.reserve(nb_grid_points); h_e.reserve(nb_grid_points);
  temps.reserve(nb_grid_points); gps.reserve(nb_grid_points);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    const double hm = number_densities[i][_Hm];
    const double he = number_densities[i][_H] * number_densities[i][_e]
                      * constants::boltzmann_k * temperature[i];

    if (hm == 0 && he == 0) continue;

    hm_nd.push_back(hm);
    h_e.push_back(he);
    temps.push_back(temperature[i]);
    gps.push_back(static_cast<int>(i));
  }

  const size_t n = gps.size();
  if (n == 0) return;

  //grow the cached staging buffers on demand
  constexpr size_t item_bytes = 3*sizeof(double) + sizeof(int);

  if (hm_batch_capacity < n)
  {
    freeHmBatchBuffers();
    gpuErrchk(cudaMalloc(&hm_batch_dev, n * item_bytes));
    gpuErrchk(cudaMallocHost(&hm_batch_host, n * item_bytes));
    hm_batch_capacity = n;
  }

  //pack everything into the pinned blob and upload with a single memcpy
  const size_t off_he = n * sizeof(double);
  const size_t off_t  = 2 * n * sizeof(double);
  const size_t off_gp = 3 * n * sizeof(double);

  std::memcpy(hm_batch_host,          hm_nd.data(), n * sizeof(double));
  std::memcpy(hm_batch_host + off_he, h_e.data(),   n * sizeof(double));
  std::memcpy(hm_batch_host + off_t,  temps.data(), n * sizeof(double));
  std::memcpy(hm_batch_host + off_gp, gps.data(),   n * sizeof(int));

  gpuErrchk(cudaMemcpy(hm_batch_dev, hm_batch_host,
    n * item_bytes, cudaMemcpyHostToDevice));

  dim3 blocks2(blocks, n);

  calcHmContinuumBatchedDevice<<<blocks2,threads>>>(
    bound_free_sigma_dev,
    spectral_grid->wavelength_list_gpu,
    (double*)hm_batch_dev,
    (double*)(hm_batch_dev + off_he),
    (double*)(hm_batch_dev + off_t),
    (int*)(hm_batch_dev + off_gp),
    nb_spectral_points,
    absorption_coeff_device);

  CUDA_CHECK_AFTER_KERNEL();
}


__host__ void GasHm::freeHmBatchBuffers()
{
  if (hm_batch_dev != nullptr)
  {
    gpuErrchk(cudaFree(hm_batch_dev));
    hm_batch_dev = nullptr;
  }

  if (hm_batch_host != nullptr)
  {
    gpuErrchk(cudaFreeHost(hm_batch_host));
    hm_batch_host = nullptr;
  }

  hm_batch_capacity = 0;
}







}
