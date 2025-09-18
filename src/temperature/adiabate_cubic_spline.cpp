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


#include "adiabate_cubic_spline.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <string>

#include "../../_deps/boost_math-src/include/boost/math/interpolators/cardinal_cubic_b_spline.hpp"


namespace bear {


AdiabateSplineTemperature::AdiabateSplineTemperature(const size_t nb_control_points_)
 : nb_control_points{nb_control_points_}
{
  if (nb_control_points < 5)
  {
    std::string error_message = "Adiabatate & Cubic B spline temperature profile requires at least 5 control points!";
    throw InvalidInput(std::string ("AdiabateSplineTemperature::AdiabateSplineTemperature"), error_message);
  }

  nb_parameters = nb_control_points + 2;
}



//calculate the temperature by using an adiabate and a cubic B spline
//uses the spline routine from the Boost library
bool AdiabateSplineTemperature::calcProfile(
  const std::vector<double>& parameters,
  const double surface_gravity,
  const std::vector<double>& pressure,
  std::vector<double>& temperature)
{
  if (parameters.size() != nb_parameters)
    std::cout << "The number of free parameters is not equal to the number of control points + 2!\n";
  
  temperature.assign(pressure.size(), 0);
  
  double pressure_rcb = parameters[0];
  
  size_t rcb_idx = findRadiavativeConvectiveBoundary(pressure_rcb, pressure);
  
  const double gamma = parameters[1];
  const double temperature_bottom = parameters[2];
  
  double temperature_rcb = temperature_bottom;

  if (rcb_idx > 0)
  { 
    std::vector<double> pressure_adiabate(pressure.begin(), pressure.begin() + rcb_idx + 1);

    pressure_adiabate.back() = pressure_rcb;
    
    addAdiabate(temperature_bottom, rcb_idx, gamma, pressure_adiabate, temperature);

    temperature_rcb = temperature[rcb_idx];
  }


  if (rcb_idx < pressure.size() - 1)
  { 
    const std::vector<double> spline_parameters(parameters.begin() + 3, parameters.end());

    addCubicSpline(
      spline_parameters, 
      rcb_idx, pressure, 
      temperature, 
      pressure_rcb, 
      temperature_rcb);
  }
  
  return checkProfile(temperature);
}



int AdiabateSplineTemperature::findRadiavativeConvectiveBoundary(
  const double rcb_pressure,
  const std::vector<double>& pressure)
{
  if (rcb_pressure <= pressure.back())
    return pressure.size() - 1;

  if (rcb_pressure >= pressure.front())
    return 0;

  for (size_t i=0; i<pressure.size(); ++i)
    if (pressure[i] <= rcb_pressure)
      return i;

  return 0;
}



void AdiabateSplineTemperature::addAdiabate(
  const double temperature_bottom,
  const unsigned int rcb_idx,
  const double gamma,
  const std::vector<double>& pressure,
  std::vector<double>& temperature)
{
  temperature[0] = temperature_bottom;
  const double kappa = 1.0 - 1.0/gamma;

  for (size_t i=1; i<pressure.size(); ++i)
  {
    // const double delta_ln_p = std::log(pressure[i]*1e6) - std::log(pressure[i-1]*1e6);
    // const double delta_ln_t = gamma * delta_ln_p;
    // const double ln_t = std::log(temperature[i-1]) + delta_ln_t;
    
    // temperature[i] = std::exp(ln_t);
    temperature[i] = temperature[0] * std::pow(pressure[i]/pressure[0], kappa);
  }
}


void AdiabateSplineTemperature::addCubicSpline(
  const std::vector<double>& parameters,
  const unsigned int rcb_idx,
  const std::vector<double>& pressure,
  std::vector<double>& temperature,
  const double pressure_rcb,
  const double temperature_rcb)
{
  double control_points_step = (std::log10(pressure_rcb) - std::log10(pressure.back())) / (nb_control_points - 1.0);

  std::vector<double> temperature_control_point(nb_control_points, 0.0);

  temperature_control_point[0] = temperature_rcb;

  for (size_t i=1; i<nb_control_points; ++i)
    temperature_control_point[i] = temperature_control_point[i-1] * parameters[i-1];

  std::reverse(temperature_control_point.begin(), temperature_control_point.end());

  boost::math::interpolators::cardinal_cubic_b_spline<double> temperature_profile (
    temperature_control_point.data(), 
    temperature_control_point.size(), 
    std::log10(pressure.back()), 
    control_points_step);

  for (size_t i=rcb_idx; i<pressure.size(); ++i)
    temperature[i] = temperature_profile(std::log10(pressure[i]));
}


}

