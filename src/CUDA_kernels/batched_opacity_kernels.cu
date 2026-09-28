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
#include <cstring>


namespace bear{


//one block column per grid point (blockIdx.y indexes a group of work items that
//share the same grid point); each thread sums all species of that grid point in
//a register and performs a single non-atomic write, which makes the result
//deterministic and avoids atomic traffic
__global__ void batchedCrossSectionsDevice(
  float** __restrict__ cs1,
  float** __restrict__ cs2,
  float** __restrict__ cs3,
  float** __restrict__ cs4,
  const float* __restrict__ temp_factors,
  const float* __restrict__ pres_factors,
  const float* __restrict__ log_number_densities,
  const int* __restrict__ group_offsets,
  const int* __restrict__ group_grid_points,
  const int nb_spectral_points,
  float* __restrict__ absorption_coeff)
{
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int g = blockIdx.y;

  if (tid >= nb_spectral_points) return;

  const int k_start = group_offsets[g];
  const int k_end   = group_offsets[g + 1];

  float sum = 0.0f;

  for (int k = k_start; k < k_end; ++k)
  {
    const float tf = temp_factors[k];
    const float pf = pres_factors[k];

    float c1 = cs1[k][tid];
    float c2 = cs2[k][tid];
    c1 = c1 + (c2 - c1) * pf;
    c2 = cs3[k][tid] + (cs4[k][tid] - cs3[k][tid]) * pf;
    const float sigma = c1 + (c2 - c1) * tf;

    sum += __exp10f(sigma + log_number_densities[k]);
  }

  absorption_coeff[group_grid_points[g] * nb_spectral_points + tid] += sum;
}


__global__ void batchedRayleighDevice(
  float** __restrict__ rayleigh_cs,
  const double* __restrict__ number_densities,
  const int* __restrict__ group_offsets,
  const int* __restrict__ group_grid_points,
  const int nb_spectral_points,
  float* __restrict__ scattering_coeff)
{
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int g = blockIdx.y;

  if (tid >= nb_spectral_points) return;

  const int k_start = group_offsets[g];
  const int k_end   = group_offsets[g + 1];

  float sum = 0.0f;

  for (int k = k_start; k < k_end; ++k)
    sum += rayleigh_cs[k][tid] * static_cast<float>(number_densities[k]);

  scattering_coeff[group_grid_points[g] * nb_spectral_points + tid] += sum;
}


void allocateBatchBuffers(BatchedDeviceBuffers& buffers, size_t capacity)
{
  if (buffers.capacity >= capacity) return;

  freeBatchBuffers(buffers);

  gpuErrchk(cudaMalloc(&buffers.cs_blob_dev,
    capacity * batch_cs_item_bytes + batch_blob_extra_bytes));
  gpuErrchk(cudaMallocHost(&buffers.cs_blob_host,
    capacity * batch_cs_item_bytes + batch_blob_extra_bytes));

  gpuErrchk(cudaMalloc(&buffers.ray_blob_dev,
    capacity * batch_ray_item_bytes + batch_blob_extra_bytes));
  gpuErrchk(cudaMallocHost(&buffers.ray_blob_host,
    capacity * batch_ray_item_bytes + batch_blob_extra_bytes));

  buffers.capacity = capacity;
}


void freeBatchBuffers(BatchedDeviceBuffers& buffers)
{
  if (buffers.capacity == 0) return;

  gpuErrchk(cudaFree(buffers.cs_blob_dev));
  gpuErrchk(cudaFreeHost(buffers.cs_blob_host));
  gpuErrchk(cudaFree(buffers.ray_blob_dev));
  gpuErrchk(cudaFreeHost(buffers.ray_blob_host));

  buffers.cs_blob_dev = nullptr;
  buffers.cs_blob_host = nullptr;
  buffers.ray_blob_dev = nullptr;
  buffers.ray_blob_host = nullptr;

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
  const size_t n = grid_points.size();
  if (n == 0) return;

  //group consecutive work items that share a grid point
  //(prepareBatchedGPU emits items with the grid point as the outer loop,
  //so equal grid points are always contiguous)
  std::vector<int> group_offsets, group_grid;
  group_offsets.reserve(n + 1);
  group_grid.reserve(n);

  for (size_t k = 0; k < n; ++k)
    if (k == 0 || grid_points[k] != grid_points[k - 1])
    {
      group_offsets.push_back(static_cast<int>(k));
      group_grid.push_back(grid_points[k]);
    }

  group_offsets.push_back(static_cast<int>(n));
  const size_t nb_groups = group_grid.size();

  //pack all metadata tightly into the pinned staging blob
  //(pointer arrays first to keep 8-byte alignment)
  const size_t off_cs2  = n * sizeof(float*);
  const size_t off_cs3  = 2 * off_cs2;
  const size_t off_cs4  = 3 * off_cs2;
  const size_t off_temp = 4 * off_cs2;
  const size_t off_pres = off_temp + n * sizeof(float);
  const size_t off_logn = off_temp + 2 * n * sizeof(float);
  const size_t off_goff = off_temp + 3 * n * sizeof(float);
  const size_t off_ggrid = off_goff + (nb_groups + 1) * sizeof(int);
  const size_t total_bytes = off_ggrid + nb_groups * sizeof(int);

  char* host = buffers.cs_blob_host;
  std::memcpy(host,            cs1_ptrs.data(), n * sizeof(float*));
  std::memcpy(host + off_cs2,  cs2_ptrs.data(), n * sizeof(float*));
  std::memcpy(host + off_cs3,  cs3_ptrs.data(), n * sizeof(float*));
  std::memcpy(host + off_cs4,  cs4_ptrs.data(), n * sizeof(float*));
  std::memcpy(host + off_temp, temp_factors.data(), n * sizeof(float));
  std::memcpy(host + off_pres, pres_factors.data(), n * sizeof(float));
  std::memcpy(host + off_logn, log_number_densities.data(), n * sizeof(float));
  std::memcpy(host + off_goff, group_offsets.data(), (nb_groups + 1) * sizeof(int));
  std::memcpy(host + off_ggrid, group_grid.data(), nb_groups * sizeof(int));

  //single upload of the whole batch
  gpuErrchk(cudaMemcpy(buffers.cs_blob_dev, host,
    total_bytes, cudaMemcpyHostToDevice));

  char* dev = buffers.cs_blob_dev;

  const int threads = 256;
  dim3 blocks((nb_spectral_points + threads - 1) / threads, nb_groups);

  batchedCrossSectionsDevice<<<blocks, threads>>>(
    (float**)dev,
    (float**)(dev + off_cs2),
    (float**)(dev + off_cs3),
    (float**)(dev + off_cs4),
    (float*)(dev + off_temp),
    (float*)(dev + off_pres),
    (float*)(dev + off_logn),
    (int*)(dev + off_goff),
    (int*)(dev + off_ggrid),
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
  const size_t n = grid_points.size();
  if (n == 0) return;

  //group consecutive work items that share a grid point (see cross-section launcher)
  std::vector<int> group_offsets, group_grid;
  group_offsets.reserve(n + 1);
  group_grid.reserve(n);

  for (size_t k = 0; k < n; ++k)
    if (k == 0 || grid_points[k] != grid_points[k - 1])
    {
      group_offsets.push_back(static_cast<int>(k));
      group_grid.push_back(grid_points[k]);
    }

  group_offsets.push_back(static_cast<int>(n));
  const size_t nb_groups = group_grid.size();

  //pack all metadata tightly into the pinned staging blob
  //(8-byte members first to keep alignment)
  const size_t off_dens = n * sizeof(float*);
  const size_t off_goff = off_dens + n * sizeof(double);
  const size_t off_ggrid = off_goff + (nb_groups + 1) * sizeof(int);
  const size_t total_bytes = off_ggrid + nb_groups * sizeof(int);

  char* host = buffers.ray_blob_host;
  std::memcpy(host,            rayleigh_ptrs.data(), n * sizeof(float*));
  std::memcpy(host + off_dens, number_densities.data(), n * sizeof(double));
  std::memcpy(host + off_goff, group_offsets.data(), (nb_groups + 1) * sizeof(int));
  std::memcpy(host + off_ggrid, group_grid.data(), nb_groups * sizeof(int));

  //single upload of the whole batch
  gpuErrchk(cudaMemcpy(buffers.ray_blob_dev, host,
    total_bytes, cudaMemcpyHostToDevice));

  char* dev = buffers.ray_blob_dev;

  const int threads = 256;
  dim3 blocks((nb_spectral_points + threads - 1) / threads, nb_groups);

  batchedRayleighDevice<<<blocks, threads>>>(
    (float**)dev,
    (double*)(dev + off_dens),
    (int*)(dev + off_goff),
    (int*)(dev + off_ggrid),
    nb_spectral_points,
    scattering_coeff_device);

  CUDA_CHECK_AFTER_KERNEL();
}


}
