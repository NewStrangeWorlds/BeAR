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
#include "toon_source_function/toon_source_function.h"
#include "../config/global_config.h"
#include "../additional/exceptions.h"


namespace bear {

//definition of the different chemistry modules with an
//identifier, a keyword to be located in the config file and a short version of the keyword
namespace rt_modules{
  enum id {scm, disort, adding_doubling, toon};
  const std::vector<std::string> description {"scm", "disort", "adding_doubling", "toon"};
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
    {
      if (parameters.size() < 1 || parameters.size() > 2)
      {
        std::string error_message = "Adding-doubling radiative transfer requires one or two parameters (number of quadrature points, optional delta_m/no_delta_m)!\n";
        throw InvalidInput(std::string ("forward_model.toml"), error_message);
      }

      //delta-M scaling of strongly forward-peaked phase functions is on by default
      bool use_delta_m = true;

      if (parameters.size() == 2)
      {
        if (parameters[1] == "no_delta_m")
          use_delta_m = false;
        else if (parameters[1] != "delta_m")
        {
          std::string error_message = "Unknown adding-doubling option \"" + parameters[1] + "\" (allowed: delta_m, no_delta_m)!\n";
          throw InvalidInput(std::string ("forward_model.toml"), error_message);
        }
      }

      return std::make_unique<AddingDoubling>(
          spectral_grid,
          std::stoi(parameters[0]),
          nb_grid_points,
          config->use_gpu,
          use_delta_m);
    }

    case rt_modules::toon :
    {
      if (parameters.size() > 1)
      {
        std::string error_message = "Toon source function radiative transfer takes at most one parameter (number of angles)!\n";
        throw InvalidInput(std::string ("forward_model.toml"), error_message);
      }

      //Gauss-Legendre angles for the intensity integration, two angles as in the scm scheme by default
      const size_t nb_angles = parameters.size() == 1 ? std::stoi(parameters[0]) : 2;

      return std::make_unique<ToonSourceFunction>(
          spectral_grid,
          nb_angles,
          nb_grid_points,
          config->use_gpu);
    }
  }


  return nullptr;
}


}
#endif
