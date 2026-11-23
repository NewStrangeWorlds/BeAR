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
//pageable transfers; the launcher computes the typed array pointers per call
struct BatchedDeviceBuffers {
  char* blob_dev = nullptr;
  char* blob_host = nullptr;   //pinned host staging buffer

  size_t capacity = 0;         //in bytes
};


void allocateBatchBuffers(BatchedDeviceBuffers& buffers, size_t nb_bytes);
void freeBatchBuffers(BatchedDeviceBuffers& buffers);

//gas absorption and Rayleigh scattering coefficients of all grid points in a
//single kernel launch; the kernel writes every element of the output arrays, 
//they do not need to be initialised
//the cross-section and Rayleigh work items are expected species by species with 
//increasing grid points
//if scattering_coeff_device is a nullptr, the scattering coefficients are added 
//to the absorption coefficients instead (total gas extinction)
void launchBatchedOpacities(
  const std::vector<float*>& cs1_ptrs,
  const std::vector<float*>& cs2_ptrs,
  const std::vector<float*>& cs3_ptrs,
  const std::vector<float*>& cs4_ptrs,
  const std::vector<float>& temp_factors,
  const std::vector<float>& pres_factors,
  const std::vector<float>& log_number_densities,
  const std::vector<int>& grid_points,
  const std::vector<float*>& rayleigh_ptrs,
  const std::vector<double>& rayleigh_number_densities,
  const std::vector<int>& rayleigh_grid_points,
  BatchedDeviceBuffers& buffers,
  int nb_spectral_points,
  int nb_grid_points,
  float* absorption_coeff_device,
  float* scattering_coeff_device);


}


#endif
