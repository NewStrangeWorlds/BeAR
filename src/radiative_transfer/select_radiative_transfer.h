// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _select_radiative_transfer_h
#define _select_radiative_transfer_h

#include <vector>
#include <string>
#include <algorithm>
#include <memory>

#include "radiative_transfer.h"

#include "short_characteristics/short_characteristics.h"
#include "disort/discrete_ordinate.h"
#include "adding_doubling/adding_doubling.h"
#include "../config/global_config.h"
#include "../additional/exceptions.h"


namespace bear {

//definition of the different chemistry modules with an
//identifier, a keyword to be located in the config file and a short version of the keyword
namespace rt_modules{
  enum id {scm, disort, adding_doubling};
  const std::vector<std::string> description {"scm", "disort", "adding_doubling"};
}



inline std::unique_ptr<RadiativeTransfer> selectRadiativeTransfer(
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
    throw InvalidInput(std::string ("forward_model.toml"), error_message);
  }


  //get the id of the chosen module
  rt_modules::id module_id = static_cast<rt_modules::id>(
    std::distance(rt_modules::description.begin(),
    it));


  //create the radiative transfer object based on the chosen module
  switch (module_id)
  {
    case rt_modules::scm :
      return std::make_unique<ShortCharacteristics>(spectral_grid);

    case rt_modules::disort :
      if (parameters.size() != 1)
      {
        std::string error_message = "Discrete ordinate radiative transfer requires exactly one parameter (number of streams)!\n";
        throw InvalidInput(std::string ("forward_model.toml"), error_message);
      }
      return std::make_unique<DiscreteOrdinates>(
          spectral_grid,
          std::stoi(parameters[0]),
          nb_grid_points,
          config->use_gpu);

    case rt_modules::adding_doubling :
      if (parameters.size() != 1)
      {
        std::string error_message = "Adding-doubling radiative transfer requires exactly one parameter (number of quadrature points)!\n";
        throw InvalidInput(std::string ("forward_model.toml"), error_message);
      }
      return std::make_unique<AddingDoubling>(
          spectral_grid,
          std::stoi(parameters[0]),
          nb_grid_points,
          config->use_gpu);
  }


  return nullptr;
}


}
#endif
