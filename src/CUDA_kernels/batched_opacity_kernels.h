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

//bytes per work item in the packed blobs
constexpr size_t batch_cs_item_bytes = 4*sizeof(float*) + 3*sizeof(float) + sizeof(int);
constexpr size_t batch_ray_item_bytes = sizeof(float*) + sizeof(double) + sizeof(int);


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
