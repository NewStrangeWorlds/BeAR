/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2022 Daniel Kitzmann
*
* Helios-r2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Helios-r2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* Helios-r2 directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#ifndef _select_temperature_h
#define _select_temperature_h


#include "temperature.h"

#include "piecewise_poly_temperature.h"


#include "../config/global_config.h"
#include "../additional/exceptions.h"

#include <vector>
#include <string>
#include <algorithm>


namespace helios {

//definition of the different chemistry modules with an
//identifier, a keyword to be located in the config file and a short version of the keyword
namespace temp_profile_modules{
  enum id {poly, constant}; 
  const std::vector<std::string> description {"poly", "const"};
}



inline Temperature* selectTemperatureProfile(const std::string profile_type, const std::vector<std::string>& parameters, 
                                             const double atmos_boundaries [2])
{
  //find the corresponding radiative transfer module to the supplied type string
  auto it = std::find(temp_profile_modules::description.begin(), temp_profile_modules::description.end(), profile_type);


  //no module is found
  if (it == temp_profile_modules::description.end())
  {
    std::string error_message = "Temperature profile type " + profile_type + " unknown!\n";
    throw ExceptionInvalidInput(std::string ("forward_model.config"), error_message);
  }


  //get the id of the chosen module
  temp_profile_modules::id module_id = static_cast<temp_profile_modules::id>(std::distance(temp_profile_modules::description.begin(), it));


  //create the radiative transfer object based on the chosen module
  Temperature* temperature_profile = nullptr;

  switch (module_id)
  {
    case temp_profile_modules::poly :
      if (parameters.size() != 2)
      {
        std::string error_message = "Piesewise polynomial temperature profile requires exactly two parameters!\n";
        throw ExceptionInvalidInput(std::string ("forward_model.config"), error_message);
      }
      {
        PiecewisePolynomialTemperature* temp = new PiecewisePolynomialTemperature(std::stoi(parameters[0]),
                                                                                  std::stoi(parameters[1]),
                                                                                  atmos_boundaries);
        temperature_profile = temp;  
      }
      break;

    case temp_profile_modules::constant :
      {
        //DiscreteOrdinates* disort = new DiscreteOrdinates(spectral_grid, std::stoi(parameters[0]), nb_grid_points, config->use_gpu); 
        //radiative_transfer = disort;
      }
      break;
  }


  return temperature_profile;
}


}
#endif

