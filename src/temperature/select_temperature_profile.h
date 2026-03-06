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


#ifndef _select_temperature_h
#define _select_temperature_h

#include <vector>
#include <string>
#include <algorithm>
#include <memory>

#include "temperature.h"

#include "../config/global_config.h"
#include "../additional/exceptions.h"

#include "piecewise_poly_temperature.h"
#include "milne_solution_temperature.h"
#include "constant_temperature.h"
#include "cubic_b_spline_temperature.h"
#include "guillot_temperature.h"
#include "adiabate_cubic_spline.h"


namespace bear {

//definition of the different chemistry modules with an
//identifier, a keyword to be located in the config file and a short version of the keyword
namespace temp_profile_modules{
  enum id {poly, milne, constant, cubicbspline, guillot, adspline};
  const std::vector<std::string> description {"poly", "milne", "const", "cubicbspline", "guillot", "adiabate_spline"};
}



inline std::unique_ptr<Temperature> selectTemperatureProfile(
  const std::string profile_type,
  const std::vector<std::string>& parameters,
  const std::vector<double>& atmos_boundaries)
{
  //find the corresponding radiative transfer module to the supplied type string
  auto it = std::find(
    temp_profile_modules::description.begin(),
    temp_profile_modules::description.end(),
    profile_type);


  //no module is found
  if (it == temp_profile_modules::description.end())
  {
    std::string error_message = "Temperature profile type " + profile_type + " unknown!\n";
    throw InvalidInput(std::string ("forward_model.config"), error_message);
  }


  //get the id of the chosen module
  temp_profile_modules::id module_id = static_cast<temp_profile_modules::id>(
    std::distance(temp_profile_modules::description.begin(), it));


  //create the temperature profile object based on the chosen module
  switch (module_id)
  {
    case temp_profile_modules::poly :
      if (parameters.size() != 2)
      {
        std::string error_message =
          "Piesewise polynomial temperature profile requires exactly two parameters!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);
      }
      return std::make_unique<PiecewisePolynomialTemperature>(
          std::stoi(parameters[0]),
          std::stoi(parameters[1]),
          atmos_boundaries);

    case temp_profile_modules::milne :
      return std::make_unique<MilneTemperature>();

    case temp_profile_modules::cubicbspline :
      if (parameters.size() != 1)
      {
        std::string error_message =
          "Cubic B spline temperature profile requires exactly one parameter!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);
      }
      return std::make_unique<CubicBSplineTemperature>(std::stoi(parameters[0]));

    case temp_profile_modules::guillot :
      if (parameters.size() != 1)
      {
        std::string error_message =
          "Guillot temperature profile requires exactly one parameter!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);
      }
      return std::make_unique<GuillotTemperature>(parameters[0]);

    case temp_profile_modules::adspline :
      if (parameters.size() != 1)
      {
        std::string error_message =
          "Adiabate Cubic B spline temperature profile requires exactly one parameter!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);
      }
      return std::make_unique<AdiabateSplineTemperature>(std::stoi(parameters[0]));

    case temp_profile_modules::constant :
      return std::make_unique<ConstantTemperature>();
  }


  return nullptr;
}


}
#endif
