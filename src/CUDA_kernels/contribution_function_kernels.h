// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _contribution_function_kernels_h
#define _contribution_function_kernels_h


#include <vector>



namespace bear{


void contributionFunctionGPU(
  float* contribution_function_dev,
  float* absorption_coeff_device, 
  double* wavenumber_list_device,
  std::vector<double>& temperature, 
  std::vector<double>& vertical_grid,
  const size_t& nb_spectral_points);

}


#endif

