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


struct BatchedDeviceBuffers {
  float** cs1_ptrs_dev = nullptr;
  float** cs2_ptrs_dev = nullptr;
  float** cs3_ptrs_dev = nullptr;
  float** cs4_ptrs_dev = nullptr;
  float* temp_factors_dev = nullptr;
  float* pres_factors_dev = nullptr;
  float* cs_log_number_densities_dev = nullptr;
  int* cs_grid_points_dev = nullptr;

  float** ray_ptrs_dev = nullptr;
  double* ray_number_densities_dev = nullptr;
  int* ray_grid_points_dev = nullptr;

  size_t capacity = 0;
};


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
