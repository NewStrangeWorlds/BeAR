// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _planck_function_kernel_h
#define _planck_function_kernel_h

#include "../additional/physical_const.h"

namespace bear{


__forceinline__ __device__ 
double planckFunction(const double temperature, const double wavenumber)
{

  return 2. * constants::planck_h * constants::light_c * constants::light_c * wavenumber*wavenumber*wavenumber 
         / ( exp(constants::planck_h * wavenumber * constants::light_c / constants::boltzmann_k / temperature) - 1.0);

}


__forceinline__ __device__ 
float planckFunction(
  const float temperature, 
  const float wavenumber_cube, 
  const float wavenumber)
{  
  constexpr float c1 = 2. * constants::planck_h * constants::light_c * constants::light_c;
  constexpr float c2 = constants::planck_h * constants::light_c / constants::boltzmann_k;
  
  // If temperature is near zero, avoid division by zero/overflow
  if (temperature < 1e-5) return 0.0f;
  
  // Use pre-calculated constants to reduce operations to: 
  // 1 division, 1 exp, 3 multiplications
  return (c1 * wavenumber_cube) / (__expf(c2 * wavenumber / temperature) - 1.0f);
}


}

#endif