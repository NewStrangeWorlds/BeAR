// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "cubic_b_spline_temperature.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <string>

#include "../../_deps/boost_math-src/include/boost/math/interpolators/cardinal_cubic_b_spline.hpp"


namespace bear {


CubicBSplineTemperature::CubicBSplineTemperature(
  const size_t nb_control_points_,
  const Parametrisation parametrisation_)
 : nb_control_points{nb_control_points_}
{
  if (nb_control_points < 5)
  {
    std::string error_message = "Cubic B spline temperature profile requires at least 5 control points!";
    throw InvalidInput(std::string ("CubicBSplineTemperature::CubicBSplineTemperature"), error_message);
  }

  parametrisation = parametrisation_;

  //parameter_names is the source of truth for the parameter count
  setControlPointNames(nb_control_points);
}



//calculate the temperature by using a cubic B spline
//uses the spline routine from the Boost library
bool CubicBSplineTemperature::calcProfile(
  const std::vector<double>& parameters,
  const double surface_gravity,
  const std::vector<double>& pressure,
  std::vector<double>& temperature)
{
  if (parameters.size() != nb_control_points)
    std::cout << "The number of free parameters is not equal to the number of control points!\n";

  double control_points_step = (std::log10(pressure[0]) - std::log10(pressure.back())) / (nb_control_points - 1.0);


  std::vector<double> temperature_control_point = controlTemperatures(parameters);

  std::reverse(temperature_control_point.begin(), temperature_control_point.end());

  boost::math::interpolators::cardinal_cubic_b_spline<double> temperature_profile (
    temperature_control_point.data(),
    temperature_control_point.size(),
    std::log10(pressure.back()),
    control_points_step);

  temperature.assign(pressure.size(), 0);

  for (size_t i=0; i<pressure.size(); ++i)
    temperature[i] = temperature_profile(std::log10(pressure[i]));


  //neglect models with too low temperatures
  bool neglect_model = false;
  
  for (auto & i : temperature)
    if (i < 50) {i = 50; neglect_model = true;}


  return neglect_model;
}


}

