// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef BATCHED_OPACITY_KERNELS_H
#define BATCHED_OPACITY_KERNELS_H

#include <vector>
#include <cstddef>


namespace bear{


//the per-work-item metadata of a batch is packed into a single blob so that it
//can be uploaded with one memcpy from pinned host memory instead of many small
//pageable transfers; the launchers compute the typed array pointers per call
struct BatchedDeviceBuffers {
  char* cs_blob_dev = nullptr;
  char* cs_blob_host = nullptr;   //pinned host staging buffer

  char* ray_blob_dev = nullptr;
  char* ray_blob_host = nullptr;  //pinned host staging buffer

  size_t capacity = 0;
};

//bytes per work item in the packed blobs; each item contributes its per-item
//metadata plus at most 2 ints of grouping data (group offset + group grid point),
//with one extra int overall for the trailing offsets entry
constexpr size_t batch_cs_item_bytes = 4*sizeof(float*) + 3*sizeof(float) + 2*sizeof(int);
constexpr size_t batch_ray_item_bytes = sizeof(float*) + sizeof(double) + 2*sizeof(int);
constexpr size_t batch_blob_extra_bytes = sizeof(int);


void allocateBatchBuffers(BatchedDeviceBuffers& buffers, size_t capacity);
void freeBatchBuffers(BatchedDeviceBuffers& buffers);

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
  float* absorption_coeff_device);

void launchBatchedRayleigh(
  const std::vector<float*>& rayleigh_ptrs,
  const std::vector<double>& number_densities,
  const std::vector<int>& grid_points,
  BatchedDeviceBuffers& buffers,
  int nb_spectral_points,
  float* scattering_coeff_device);


}


#endif
