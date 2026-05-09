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


#include "pchip_temperature.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <string>

#include "../../_deps/boost_math-src/include/boost/math/interpolators/pchip.hpp"


namespace bear {


PchipTemperature::PchipTemperature(const size_t nb_control_points_)
 : nb_control_points{nb_control_points_}
{
  if (nb_control_points < 4)
  {
    std::string error_message = "PCHIP temperature profile requires at least 4 control points!";
    throw InvalidInput(std::string ("PchipTemperature::PchipTemperature"), error_message);
  }

  nb_parameters = nb_control_points;
}



bool PchipTemperature::calcProfile(
  const std::vector<double>& parameters,
  const double surface_gravity,
  const std::vector<double>& pressure,
  std::vector<double>& temperature)
{
  if (parameters.size() != nb_control_points)
    std::cout << "The number of free parameters is not equal to the number of control points!\n";

  double control_points_step = (std::log10(pressure[0]) - std::log10(pressure.back())) / (nb_control_points - 1.0);

  std::vector<double> temperature_control_point(nb_control_points, 0.0);

  for (size_t i = 0; i < nb_parameters; ++i)
    temperature_control_point[i] = parameters[i];

  std::reverse(temperature_control_point.begin(), temperature_control_point.end());

  std::vector<double> x_knots(nb_control_points);
  for (size_t i = 0; i < nb_control_points; ++i)
    x_knots[i] = std::log10(pressure.back()) + i * control_points_step;

  auto temperature_profile = boost::math::interpolators::pchip<std::vector<double>>(
    std::move(x_knots), std::move(temperature_control_point));

  temperature.assign(pressure.size(), 0);

  for (size_t i = 0; i < pressure.size(); ++i)
    temperature[i] = temperature_profile(std::log10(pressure[i]));

  bool neglect_model = false;

  for (auto & t : temperature)
    if (t < 50) {t = 50; neglect_model = true;}

  return neglect_model;
}


}
