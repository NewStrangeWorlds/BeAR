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


#include "batched_opacity_kernels.h"
#include "error_check.h"
#include "data_management_kernels.h"

#include <cstdio>


namespace bear{


__global__ void batchedCrossSectionsDevice(
  float** __restrict__ cs1,
  float** __restrict__ cs2,
  float** __restrict__ cs3,
  float** __restrict__ cs4,
  const float* __restrict__ temp_factors,
  const float* __restrict__ pres_factors,
  const float* __restrict__ log_number_densities,
  const int* __restrict__ grid_points,
  const int nb_spectral_points,
  float* __restrict__ absorption_coeff)
{
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int wi = blockIdx.y;

  if (tid >= nb_spectral_points) return;

  const float tf = temp_factors[wi];
  const float pf = pres_factors[wi];

  float c1 = cs1[wi][tid];
  float c2 = cs2[wi][tid];
  c1 = c1 + (c2 - c1) * pf;
  c2 = cs3[wi][tid] + (cs4[wi][tid] - cs3[wi][tid]) * pf;
  float sigma = c1 + (c2 - c1) * tf;

  atomicAdd(
    &absorption_coeff[grid_points[wi] * nb_spectral_points + tid],
    __exp10f(sigma + log_number_densities[wi]));
}


__global__ void batchedRayleighDevice(
  float** __restrict__ rayleigh_cs,
  const double* __restrict__ number_densities,
  const int* __restrict__ grid_points,
  const int nb_spectral_points,
  float* __restrict__ scattering_coeff)
{
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int wi = blockIdx.y;

  if (tid >= nb_spectral_points) return;

  atomicAdd(
    &scattering_coeff[grid_points[wi] * nb_spectral_points + tid],
    rayleigh_cs[wi][tid] * static_cast<float>(number_densities[wi]));
}


void allocateBatchBuffers(BatchedDeviceBuffers& buffers, size_t capacity)
{
  if (buffers.capacity >= capacity) return;

  freeBatchBuffers(buffers);

  allocateOnDevice(buffers.cs1_ptrs_dev, capacity);
  allocateOnDevice(buffers.cs2_ptrs_dev, capacity);
  allocateOnDevice(buffers.cs3_ptrs_dev, capacity);
  allocateOnDevice(buffers.cs4_ptrs_dev, capacity);
  allocateOnDevice(buffers.temp_factors_dev, capacity);
  allocateOnDevice(buffers.pres_factors_dev, capacity);
  allocateOnDevice(buffers.cs_log_number_densities_dev, capacity);
  allocateOnDevice(buffers.cs_grid_points_dev, capacity);

  allocateOnDevice(buffers.ray_ptrs_dev, capacity);
  allocateOnDevice(buffers.ray_number_densities_dev, capacity);
  allocateOnDevice(buffers.ray_grid_points_dev, capacity);

  buffers.capacity = capacity;
}


void freeBatchBuffers(BatchedDeviceBuffers& buffers)
{
  if (buffers.capacity == 0) return;

  deleteFromDevice(buffers.cs1_ptrs_dev);
  deleteFromDevice(buffers.cs2_ptrs_dev);
  deleteFromDevice(buffers.cs3_ptrs_dev);
  deleteFromDevice(buffers.cs4_ptrs_dev);
  deleteFromDevice(buffers.temp_factors_dev);
  deleteFromDevice(buffers.pres_factors_dev);
  deleteFromDevice(buffers.cs_log_number_densities_dev);
  deleteFromDevice(buffers.cs_grid_points_dev);

  deleteFromDevice(buffers.ray_ptrs_dev);
  deleteFromDevice(buffers.ray_number_densities_dev);
  deleteFromDevice(buffers.ray_grid_points_dev);

  buffers.capacity = 0;
}


void launchBatchedCrossSections(
  const std::vector<float*>& cs1_ptrs,
  const std::vector<float*>& cs2_ptrs,
  const std::vector<float*>& cs3_ptrs,
  const std::vector<float*>& cs4_ptrs,
  const std::vector<float>& temp_factors,
  const std::vector<float>& pres_factors,
  const std::vector<float>& log_number_densities,
  const std::vector<int>& grid_points,
  BatchedDeviceBuffers& buffers,
  int nb_spectral_points,
  float* absorption_coeff_device)
{
  const size_t nb_work_items = grid_points.size();
  if (nb_work_items == 0) return;

  gpuErrchk(cudaMemcpy(buffers.cs1_ptrs_dev, cs1_ptrs.data(),
    nb_work_items * sizeof(float*), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(buffers.cs2_ptrs_dev, cs2_ptrs.data(),
    nb_work_items * sizeof(float*), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(buffers.cs3_ptrs_dev, cs3_ptrs.data(),
    nb_work_items * sizeof(float*), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(buffers.cs4_ptrs_dev, cs4_ptrs.data(),
    nb_work_items * sizeof(float*), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(buffers.temp_factors_dev, temp_factors.data(),
    nb_work_items * sizeof(float), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(buffers.pres_factors_dev, pres_factors.data(),
    nb_work_items * sizeof(float), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(buffers.cs_log_number_densities_dev, log_number_densities.data(),
    nb_work_items * sizeof(float), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(buffers.cs_grid_points_dev, grid_points.data(),
    nb_work_items * sizeof(int), cudaMemcpyHostToDevice));

  const int threads = 256;
  dim3 blocks((nb_spectral_points + threads - 1) / threads, nb_work_items);

  batchedCrossSectionsDevice<<<blocks, threads>>>(
    buffers.cs1_ptrs_dev,
    buffers.cs2_ptrs_dev,
    buffers.cs3_ptrs_dev,
    buffers.cs4_ptrs_dev,
    buffers.temp_factors_dev,
    buffers.pres_factors_dev,
    buffers.cs_log_number_densities_dev,
    buffers.cs_grid_points_dev,
    nb_spectral_points,
    absorption_coeff_device);

  CUDA_CHECK_AFTER_KERNEL();
}


void launchBatchedRayleigh(
  const std::vector<float*>& rayleigh_ptrs,
  const std::vector<double>& number_densities,
  const std::vector<int>& grid_points,
  BatchedDeviceBuffers& buffers,
  int nb_spectral_points,
  float* scattering_coeff_device)
{
  const size_t nb_work_items = grid_points.size();
  if (nb_work_items == 0) return;

  gpuErrchk(cudaMemcpy(buffers.ray_ptrs_dev, rayleigh_ptrs.data(),
    nb_work_items * sizeof(float*), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(buffers.ray_number_densities_dev, number_densities.data(),
    nb_work_items * sizeof(double), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(buffers.ray_grid_points_dev, grid_points.data(),
    nb_work_items * sizeof(int), cudaMemcpyHostToDevice));

  const int threads = 256;
  dim3 blocks((nb_spectral_points + threads - 1) / threads, nb_work_items);

  batchedRayleighDevice<<<blocks, threads>>>(
    buffers.ray_ptrs_dev,
    buffers.ray_number_densities_dev,
    buffers.ray_grid_points_dev,
    nb_spectral_points,
    scattering_coeff_device);

  CUDA_CHECK_AFTER_KERNEL();
}


}
