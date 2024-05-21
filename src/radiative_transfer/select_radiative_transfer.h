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


#ifndef _select_radiative_transfer_h
#define _select_radiative_transfer_h

#include <vector>
#include <string>
#include <algorithm>

#include "radiative_transfer.h"

#include "short_characteristics.h"
#include "discrete_ordinate.h"
#include "../config/global_config.h"
#include "../additional/exceptions.h"


namespace helios {

//definition of the different chemistry modules with an
//identifier, a keyword to be located in the config file and a short version of the keyword
namespace rt_modules{
  enum id {scm, disort}; 
  const std::vector<std::string> description {"scm", "disort"};
}



inline RadiativeTransfer* selectRadiativeTransfer(
  const std::string rt_type,
  const std::vector<std::string>& parameters,
  const size_t nb_grid_points, 
  GlobalConfig* config,
  SpectralGrid* spectral_grid)
{
  //find the corresponding radiative transfer module to the supplied "type" string
  auto it = std::find(
    rt_modules::description.begin(),
    rt_modules::description.end(),
    rt_type);


  //no module is found
  if (it == rt_modules::description.end())
  {
    std::string error_message = "Radiative transfer type " + rt_type + " unknown!\n";
    throw InvalidInput(std::string ("forward_model.config"), error_message);
  }


  //get the id of the chosen module
  rt_modules::id module_id = static_cast<rt_modules::id>(
    std::distance(rt_modules::description.begin(),
    it));


  //create the radiative transfer object based on the chosen module
  RadiativeTransfer* radiative_transfer = nullptr;

  switch (module_id)
  {
    case rt_modules::scm :
      {
        ShortCharacteristics* scm = new ShortCharacteristics(spectral_grid);
        radiative_transfer = scm; 
      }
      break;

    case rt_modules::disort :
      if (parameters.size() != 1)
      {
        std::string error_message = "Discrete ordinate radiative transfer requires exactly one parameter (number of streams)!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);
      }
      {
        DiscreteOrdinates* disort = new DiscreteOrdinates(
          spectral_grid, 
          std::stoi(parameters[0]), 
          nb_grid_points, 
          config->use_gpu); 
        radiative_transfer = disort;
      }
      break;
  }


  return radiative_transfer;
}


}
#endif

