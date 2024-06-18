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


#include "milne_solution_temperature.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"

#include <algorithm>
#include <vector>
#include <cmath>


namespace helios {



//calculate the temperature based on the analytical Milne solution
bool MilneTemperature::calcProfile(
  const std::vector<double>& parameters, 
  const double surface_gravity,
  const std::vector<double>& pressure, 
  std::vector<double>& temperature)
{
  temperature.assign(pressure.size(), 0);

  const double kappa_ross = parameters[0];
  const double t_eff = parameters[1];

  for (size_t i=0; i<temperature.size(); ++i)
  {
    const double tau_ross = kappa_ross * pressure[i] * 1e6 / surface_gravity;
    
    temperature[i] = 3.0/4.0 * std::pow(t_eff, 4) * (hopfFunction(tau_ross) + tau_ross);
    temperature[i] = std::pow(temperature[i], 0.25);
  }

  //neglect models with too low temperatures
  bool neglect_model = false;
  
  for (auto & i : temperature)
    if (i < 50) {i = 50; neglect_model = true;}


  return neglect_model;
}


double MilneTemperature::hopfFunction(const double optical_depth)
{
  const double log_optical_depth = std::log10(optical_depth);

  if (optical_depth < 0.01)
    return 0.577351 + (optical_depth - 0.0) * (0.588236 - 0.577351)/(0.01 - 0.0);

  if (optical_depth > 5)
    return  0.710398 + (log_optical_depth - log10(5.0)) * (0.710446 - 0.710398)/(10000.0 - log10(5.0));


  return (  fit_p[0] * std::pow(log_optical_depth, 4)
          + fit_p[1] * std::pow(log_optical_depth, 3)
          + fit_p[2] * std::pow(log_optical_depth, 2)
          + fit_p[3] * log_optical_depth
          + fit_p[4]) /
           (std::pow(log_optical_depth, 4)
            + fit_q[0] * std::pow(log_optical_depth, 3)
            + fit_q[1] * std::pow(log_optical_depth, 2)
            + fit_q[2] * log_optical_depth
            + fit_q[3]);
}


}

