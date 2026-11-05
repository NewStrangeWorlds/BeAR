// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _cross_section_kernels_h
#define _cross_section_kernels_h

#include <vector>

namespace bear{

void initCrossSectionsHost(
  const size_t nb_points, 
  float* absorption_coeff_device);

void checkCrossSectionsHost(
  const size_t nb_spectral_points, 
  const size_t nb_grid_points,
  float* absorption_coeff_device);

}

#endif
