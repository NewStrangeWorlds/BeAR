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


#ifndef _select_module_h
#define _select_module_h

#include <vector>
#include <string>
#include <algorithm>
#include <memory>

#include "module.h"

#include "../../config/global_config.h"
#include "../../additional/exceptions.h"
#include "../../spectral_grid/spectral_grid.h"

#include "stellar_contamination/stellar_contamination.h"
#include "velocity_broadening/velocity_broadening.h"



namespace bear {

//definition of the different chemistry modules with an
//identifier, a keyword to be located in the config file and a short version of the keyword
namespace modules{
  enum id {stellar_contamination, velocity_broadening};
  const std::vector<std::string> description {"stellar_contamination", "velocity_broadening"};
}



inline std::unique_ptr<Module> selectModule(
  const std::string model_type,
  const std::vector<std::string>& parameters,
  SpectralGrid* spectral_grid)
{
  //find the corresponding radiative transfer module to the supplied type string
  auto it = std::find(
    modules::description.begin(),
    modules::description.end(),
    model_type);


  //no module is found
  if (it == modules::description.end())
  {
    std::string error_message = "Module " + model_type + " unknown!\n";
    throw InvalidInput(std::string ("forward_model.config"), error_message);
  }


  //get the id of the chosen module
  modules::id module_id = static_cast<modules::id>(
    std::distance(modules::description.begin(), it));


  //create the module object based on the chosen module
  switch (module_id)
  {
    case modules::stellar_contamination :
      if (parameters.size() < 1)
      {
        std::string error_message =
          "Stellar activity module requires at least one parameter!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);
      }
      return std::make_unique<StellarContamination>(
          parameters,
          spectral_grid);
    
    case modules::velocity_broadening :
      if (parameters.size() != 0)
      {
        std::string error_message =
          "Velocity broadening module does not require any parameters!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);
      }
      return std::make_unique<VelocityBroadening>(
          parameters,
          spectral_grid);
  }

  return nullptr;
}


}
#endif
