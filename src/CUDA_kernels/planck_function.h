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


//single-precision Planck function with a full-precision exponential; used by the
//scattering solvers, where the fast __expf (error growing with its argument)
//would limit the agreement with the double-precision CPU path
__forceinline__ __device__ 
float planckFunctionAccurate(
  const float temperature, 
  const float wavenumber_cube, 
  const float wavenumber)
{  
  constexpr float c1 = 2. * constants::planck_h * constants::light_c * constants::light_c;
  constexpr float c2 = constants::planck_h * constants::light_c / constants::boltzmann_k;

  if (temperature < 1e-5) return 0.0f;

  return (c1 * wavenumber_cube) / expm1f(c2 * wavenumber / temperature);
}


//difference B(temperature_1) - B(temperature_2) of the single-precision Planck function
//with full relative accuracy, also when the two values are close: with u = c2 nu / T,
//  B_1 - B_2 = c1 nu^3 expm1(u_2 - u_1) / ((1 - exp(-u_1)) expm1(u_2))
//the direct difference is used where the second term is negligible anyway
__forceinline__ __device__ 
float planckDifferenceAccurate(
  const float temperature_1,
  const float temperature_2,
  const float wavenumber_cube, 
  const float wavenumber)
{  
  constexpr float c1 = 2. * constants::planck_h * constants::light_c * constants::light_c;
  constexpr float c2 = constants::planck_h * constants::light_c / constants::boltzmann_k;

  if (temperature_1 < 1e-5 || temperature_2 < 1e-5)
    return planckFunctionAccurate(temperature_1, wavenumber_cube, wavenumber) 
         - planckFunctionAccurate(temperature_2, wavenumber_cube, wavenumber);

  const float u1 = c2 * wavenumber / temperature_1;
  const float u2 = c2 * wavenumber / temperature_2;
  const float du = c2 * wavenumber * (temperature_1 - temperature_2) / (temperature_1 * temperature_2);

  if (fabsf(du) > 60.0f || u2 > 80.0f)
    return planckFunctionAccurate(temperature_1, wavenumber_cube, wavenumber) 
         - planckFunctionAccurate(temperature_2, wavenumber_cube, wavenumber);

  return c1 * wavenumber_cube * expm1f(du) / ((1.0f - expf(-u1)) * expm1f(u2));
}


}

#endif