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
double planckFunctionOld(
  const double temperature, 
  const double wavenumber_cube, 
  const double wavenumber)
{  
  constexpr double c1 = 2. * constants::planck_h * constants::light_c * constants::light_c;
  constexpr double c2 = constants::planck_h * constants::light_c / constants::boltzmann_k;
  
  // If temperature is near zero, avoid division by zero/overflow
  if (temperature < 1e-5) return 0.0;
  
  // Use pre-calculated constants to reduce operations to: 
  // 1 division, 1 exp, 3 multiplications
  return (c1 * wavenumber_cube) / (exp(c2 * wavenumber / temperature) - 1.0);
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