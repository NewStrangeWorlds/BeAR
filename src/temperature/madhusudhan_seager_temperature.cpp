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


#include "madhusudhan_seager_temperature.h"
#include "../additional/exceptions.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <limits>


namespace bear {


MadhusudhanSeagerTemperature::MadhusudhanSeagerTemperature()
{
  nb_parameters = 6;

  std::cout << "\n- Temperature profile: Madhusudhan & Seager (2009)\n\n";
}


// Madhusudhan & Seager (2009) three-layer T-P profile.
// Matches the TP_MS() function in CHIMERA (Line et al. 2021).
//
// In BeAR, pressure[0] is the highest pressure (bottom of atmosphere) and
// pressure.back() is the lowest pressure (top), so P0 = pressure.back().
bool MadhusudhanSeagerTemperature::calcProfile(
  const std::vector<double>& parameters,
  const double /*surface_gravity*/,
  const std::vector<double>& pressure,
  std::vector<double>& temperature)
{
  const double T0     = parameters[0];
  const double log_P1 = parameters[1];
  const double log_P2 = parameters[2];
  const double log_P3 = parameters[3];
  const double alpha1 = parameters[4];
  const double alpha2 = parameters[5];

  const double P1 = std::pow(10.0, log_P1);
  const double P2 = std::pow(10.0, log_P2);
  const double P3 = std::pow(10.0, log_P3);

  // Top-of-atmosphere reference pressure
  const double P0 = pressure.back();

  // Continuity condition at P1: T2 such that layer-2 formula matches layer-1 at P1
  const double ln_P1_P0_over_a1 = std::log(P1 / P0) / alpha1;
  const double ln_P1_P2_over_a2 = std::log(P1 / P2) / alpha2;
  const double T2 = ln_P1_P0_over_a1 * ln_P1_P0_over_a1 + T0
                  - ln_P1_P2_over_a2 * ln_P1_P2_over_a2;

  // Isothermal temperature for layer 3: value of layer-2 formula at P3
  const double ln_P3_P2_over_a2 = std::log(P3 / P2) / alpha2;
  const double T3 = ln_P3_P2_over_a2 * ln_P3_P2_over_a2 + T2;

  temperature.assign(pressure.size(), 0.0);

  for (size_t i = 0; i < pressure.size(); ++i)
  {
    const double P = pressure[i];

    if (P < P1)
    {
      // Layer 1: upper atmosphere
      const double x = std::log(P / P0) / alpha1;
      temperature[i] = x * x + T0;
    }
    else if (P < P3)
    {
      // Layer 2
      const double x = std::log(P / P2) / alpha2;
      temperature[i] = x * x + T2;
    }
    else
    {
      // Layer 3: isothermal deep atmosphere
      temperature[i] = T3;
    }
  }

  return checkProfile(temperature);
}


}
