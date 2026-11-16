// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "batched_opacity_kernels.h"
#include "error_check.h"
#include "data_management_kernels.h"

#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <numeric>


namespace bear{


namespace {

//threads per block of the cross-section kernel (spectral points per block), the
//number of consecutive grid points a block accumulates in registers and the number
//of distinct tables a run of work items may reference
constexpr int batch_cs_threads = 256;
constexpr int batch_cs_group = 8;
constexpr int batch_cs_run_tables = 16;

//logarithm of the number density of a padding item, makes its opacity vanish
constexpr float batch_cs_padding_log_density = -1e30f;

//packed work item: the four table slots (4 bits each) in .x, the temperature and
//pressure interpolation factors and the logarithm of the number density in .y, .z, .w
__host__ __device__ __forceinline__
unsigned int packSlots(const int slot1, const int slot2, const int slot3, const int slot4)
{
  return static_cast<unsigned int>(slot1) | (static_cast<unsigned int>(slot2) << 4)
       | (static_cast<unsigned int>(slot3) << 8) | (static_cast<unsigned int>(slot4) << 12);
}

}


//predicated table read: the value at the spectral point of this thread if the
//table is used, otherwise zero (as inline PTX so that the read is issued without
//a branch and all reads of a run are in flight at the same time)
__device__ __forceinline__
float tableLoad(const bool used, const float* __restrict__ q, const int tid)
{
  float value = 0.0f;

  asm("{\n\t.reg .pred p;\n\tsetp.ne.u32 p, %1, 0;\n\t@p ld.global.nc.f32 %0, [%2];\n\t}"
      : "+f"(value) : "r"(static_cast<unsigned int>(used)), "l"(q + tid));

  return value;
}



//a block handles batch_cs_threads consecutive spectral points for a group of
//batch_cs_group consecutive grid points; its work items come in runs of one species
//over the grid points of the group (padded with vanishing items where a species
//has no cross sections at a grid point); neighbouring grid points of one species
//mostly share their four tabulated tables, so a run references only a few distinct
//tables: these are read once per thread, with all reads of the run in flight at
//the same time, into the thread's own column of shared memory, from which the
//items interpolate; the sums per grid point stay in registers, their order
//(species by species) is deterministic
__global__ void batchedCrossSectionsDevice(
  const float* const* __restrict__ run_tables,
  const int* __restrict__ run_nb_tables,
  const float4* __restrict__ items,
  const int* __restrict__ group_run_offsets,
  const int nb_spectral_points,
  const int nb_grid_points,
  float* __restrict__ absorption_coeff)
{
  __shared__ float s_tables[batch_cs_run_tables * batch_cs_threads];

  const int tid = blockIdx.x * batch_cs_threads + threadIdx.x;
  const int group = blockIdx.y;

  if (tid >= nb_spectral_points) return;

  float* table = s_tables + threadIdx.x;

  float sum[batch_cs_group];

  #pragma unroll
  for (int j = 0; j < batch_cs_group; ++j)
    sum[j] = 0.0f;

  const int run_end = group_run_offsets[group + 1];

  for (int run = group_run_offsets[group]; run < run_end; ++run)
  {
    const float* const* tables = run_tables + run * batch_cs_run_tables;
    const int nb_tables = run_nb_tables[run];

    float loaded[batch_cs_run_tables];

    #pragma unroll
    for (int u = 0; u < batch_cs_run_tables; ++u)
      loaded[u] = tableLoad(u < nb_tables, tables[u], tid);

    #pragma unroll
    for (int u = 0; u < batch_cs_run_tables; ++u)
      table[u * batch_cs_threads] = loaded[u];

    #pragma unroll
    for (int j = 0; j < batch_cs_group; ++j)
    {
      const float4 item = items[run * batch_cs_group + j];
      const unsigned int slots = __float_as_uint(item.x);

      const float w1 = table[(slots & 15u) * batch_cs_threads];
      const float w2 = table[((slots >> 4) & 15u) * batch_cs_threads];
      const float w3 = table[((slots >> 8) & 15u) * batch_cs_threads];
      const float w4 = table[(slots >> 12) * batch_cs_threads];

      const float c1 = w1 + (w2 - w1) * item.z;
      const float c2 = w3 + (w4 - w3) * item.z;
      const float sigma = c1 + (c2 - c1) * item.y;

      sum[j] += __exp10f(sigma + item.w);
    }
  }

  const int gp_start = group * batch_cs_group;

  #pragma unroll
  for (int j = 0; j < batch_cs_group; ++j)
    if (gp_start + j < nb_grid_points)
      absorption_coeff[(gp_start + j) * nb_spectral_points + tid] += sum[j];
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
  int nb_grid_points,
  float* absorption_coeff_device)
{
  const size_t n = grid_points.size();
  if (n == 0) return;

  //the items arrive species by species with increasing grid point; the species
  //boundaries are where the grid point does not increase
  std::vector<size_t> species_starts;

  for (size_t k = 0; k < n; ++k)
    if (k == 0 || grid_points[k] <= grid_points[k-1])
      species_starts.push_back(k);

  species_starts.push_back(n);

  const size_t nb_species = species_starts.size() - 1;
  const int nb_groups = (nb_grid_points + batch_cs_group - 1) / batch_cs_group;

  //runs of one species over the grid points of a group, padded to batch_cs_group
  //items with vanishing items; a run references at most batch_cs_run_tables
  //distinct tables, a species is split into further runs when it needs more
  std::vector<const float*> run_tables;
  std::vector<int> run_nb_tables;
  std::vector<float4> items;
  std::vector<int> group_run_offsets(nb_groups + 1, 0);

  std::vector<size_t> species_pos(nb_species);
  for (size_t s = 0; s < nb_species; ++s) species_pos[s] = species_starts[s];

  const float4 padding_item = {0.0f, 0.0f, 0.0f, batch_cs_padding_log_density};

  auto openRun = [&]()
  {
    run_tables.resize(run_tables.size() + batch_cs_run_tables, nullptr);
    run_nb_tables.push_back(0);
    items.resize(items.size() + batch_cs_group, padding_item);
  };

  for (int g = 0; g < nb_groups; ++g)
  {
    const int gp_start = g * batch_cs_group;
    const int gp_end = gp_start + batch_cs_group;

    for (size_t s = 0; s < nb_species; ++s)
    {
      //items of this species inside the group (grid points increase within a species)
      const size_t first = species_pos[s];
      size_t last = first;

      while (last < species_starts[s+1] && grid_points[last] < gp_end) ++last;

      species_pos[s] = last;

      if (last == first) continue;

      openRun();
      group_run_offsets[g+1]++;

      for (size_t k = first; k < last; ++k)
      {
        const float* item_tables[4] = {cs1_ptrs[k], cs2_ptrs[k], cs3_ptrs[k], cs4_ptrs[k]};

        //distinct tables of the item not yet in the run
        int nb_new = 0;

        for (int j = 0; j < 4; ++j)
        {
          bool found = false;

          for (int i = 0; i < j; ++i)
            if (item_tables[i] == item_tables[j]) found = true;

          const float* const* tables = run_tables.data() + run_tables.size() - batch_cs_run_tables;

          for (int u = 0; u < run_nb_tables.back() && !found; ++u)
            if (tables[u] == item_tables[j]) found = true;

          if (!found) nb_new++;
        }

        if (run_nb_tables.back() + nb_new > batch_cs_run_tables)
        {
          openRun();
          group_run_offsets[g+1]++;
        }

        const float** tables = run_tables.data() + run_tables.size() - batch_cs_run_tables;
        int& nb_tables = run_nb_tables.back();
        int slots[4];

        for (int j = 0; j < 4; ++j)
        {
          slots[j] = -1;

          for (int u = 0; u < nb_tables; ++u)
            if (tables[u] == item_tables[j]) { slots[j] = u; break; }

          if (slots[j] < 0)
          {
            slots[j] = nb_tables;
            tables[nb_tables++] = item_tables[j];
          }
        }

        float4 item;
        const unsigned int packed = packSlots(slots[0], slots[1], slots[2], slots[3]);
        std::memcpy(&item.x, &packed, sizeof(float));
        item.y = temp_factors[k];
        item.z = pres_factors[k];
        item.w = log_number_densities[k];

        items[items.size() - batch_cs_group + (grid_points[k] - gp_start)] = item;
      }
    }
  }

  for (int g = 0; g < nb_groups; ++g)
    group_run_offsets[g+1] += group_run_offsets[g];

  const size_t nb_runs = run_nb_tables.size();

  //unused table slots of a run point to its first table (they are never read)
  for (size_t r = 0; r < nb_runs; ++r)
    for (int u = run_nb_tables[r]; u < batch_cs_run_tables; ++u)
      run_tables[r * batch_cs_run_tables + u] = run_tables[r * batch_cs_run_tables];

  allocateBatchBuffers(buffers, nb_runs + nb_groups + 2);

  //pack all metadata tightly into the pinned staging blob
  //(the 16-byte items first, then the pointer array, to keep the alignment)
  const size_t off_tab  = items.size() * sizeof(float4);
  const size_t off_ntab = off_tab + run_tables.size() * sizeof(float*);
  const size_t off_goff = off_ntab + nb_runs * sizeof(int);
  const size_t total_bytes = off_goff + (nb_groups + 1) * sizeof(int);

  char* host = buffers.cs_blob_host;
  std::memcpy(host,            items.data(), items.size() * sizeof(float4));
  std::memcpy(host + off_tab,  run_tables.data(), run_tables.size() * sizeof(float*));
  std::memcpy(host + off_ntab, run_nb_tables.data(), nb_runs * sizeof(int));
  std::memcpy(host + off_goff, group_run_offsets.data(), (nb_groups + 1) * sizeof(int));

  //single upload of the whole batch
  gpuErrchk(cudaMemcpy(buffers.cs_blob_dev, host,
    total_bytes, cudaMemcpyHostToDevice));

  const char* dev = buffers.cs_blob_dev;

  const dim3 blocks((nb_spectral_points + batch_cs_threads - 1) / batch_cs_threads, nb_groups);

  batchedCrossSectionsDevice<<<blocks, batch_cs_threads>>>(
    (const float* const*)(dev + off_tab),
    (const int*)(dev + off_ntab),
    (const float4*)dev,
    (const int*)(dev + off_goff),
    nb_spectral_points,
    nb_grid_points,
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

  //group the work items by grid point (stable sort keeps the species order
  //within a grid point, so the summation order is independent of the input order)
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(),
    [&](const int a, const int b) { return grid_points[a] < grid_points[b]; });

  std::vector<int> group_offsets, group_grid;
  group_offsets.reserve(n + 1);
  group_grid.reserve(n);

  for (size_t i = 0; i < n; ++i)
    if (i == 0 || grid_points[order[i]] != grid_points[order[i-1]])
    {
      group_offsets.push_back(static_cast<int>(i));
      group_grid.push_back(grid_points[order[i]]);
    }

  group_offsets.push_back(static_cast<int>(n));
  const size_t nb_groups = group_grid.size();

  allocateBatchBuffers(buffers, n);

  //pack all metadata tightly into the pinned staging blob
  //(8-byte members first to keep alignment)
  const size_t off_dens = n * sizeof(float*);
  const size_t off_goff = off_dens + n * sizeof(double);
  const size_t off_ggrid = off_goff + (nb_groups + 1) * sizeof(int);
  const size_t total_bytes = off_ggrid + nb_groups * sizeof(int);

  char* host = buffers.ray_blob_host;

  float** h_ptrs = reinterpret_cast<float**>(host);
  double* h_dens = reinterpret_cast<double*>(host + off_dens);

  for (size_t i = 0; i < n; ++i)
  {
    h_ptrs[i] = rayleigh_ptrs[order[i]];
    h_dens[i] = number_densities[order[i]];
  }

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
