// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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
#include "phase_resolved_broadening/phase_resolved_broadening.h"



namespace bear {

//definition of the different chemistry modules with an
//identifier, a keyword to be located in the config file and a short version of the keyword
namespace modules{
  enum id {stellar_contamination, velocity_broadening, phase_resolved_broadening};
  const std::vector<std::string> description {"stellar_contamination", "velocity_broadening", "phase_resolved_broadening"};
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
    throw InvalidInput(std::string ("forward_model.toml"), error_message);
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
        throw InvalidInput(std::string ("forward_model.toml"), error_message);
      }
      return std::make_unique<StellarContamination>(
          parameters,
          spectral_grid);
    
    case modules::velocity_broadening :
      if (parameters.size() != 0)
      {
        std::string error_message =
          "Velocity broadening module does not require any parameters!\n";
        throw InvalidInput(std::string ("forward_model.toml"), error_message);
      }
      return std::make_unique<VelocityBroadening>(
          parameters,
          spectral_grid);

    case modules::phase_resolved_broadening :
      if (parameters.size() != 0)
      {
        std::string error_message =
          "Phase-resolved broadening module does not require any parameters!\n";
        throw InvalidInput(std::string ("forward_model.toml"), error_message);
      }
      return std::make_unique<PhaseResolvedBroadening>(
          parameters,
          spectral_grid);
  }

  return nullptr;
}


}
#endif
