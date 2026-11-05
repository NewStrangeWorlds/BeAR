// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _select_stellar_model_h
#define _select_stellar_model_h

#include <vector>
#include <string>
#include <algorithm>
#include <memory>

#include "stellar_spectrum.h"

#include "../../config/global_config.h"
#include "../../additional/exceptions.h"
#include "../../spectral_grid/spectral_grid.h"

#include "star_blackbody.h"
#include "star_file_spectrum.h"
#include "stellar_spectrum_grid.h"


namespace bear {

//definition of the different chemistry modules with an
//identifier, a keyword to be located in the config file and a short version of the keyword
namespace stellar_modules{
  enum id {blackbody, file, grid};
  const std::vector<std::string> description {"blackbody", "file", "grid"};
}



inline std::unique_ptr<StellarSpectrumModel> selectStellarModel(
  const std::string model_type,
  const std::vector<std::string>& parameters,
  SpectralGrid* spectral_grid)
{
  //find the corresponding radiative transfer module to the supplied type string
  auto it = std::find(
    stellar_modules::description.begin(),
    stellar_modules::description.end(),
    model_type);


  //no module is found
  if (it == stellar_modules::description.end())
  {
    std::string error_message = "Stellar model " + model_type + " unknown!\n";
    throw InvalidInput(std::string ("forward_model.toml"), error_message);
  }


  //get the id of the chosen module
  stellar_modules::id module_id = static_cast<stellar_modules::id>(
    std::distance(stellar_modules::description.begin(), it));


  //create the stellar spectrum object based on the chosen module
  switch (module_id)
  {
    case stellar_modules::blackbody :
      if (parameters.size() != 0)
      {
        std::string error_message =
          "Stellar blackbody model requires no parameters!\n";
        throw InvalidInput(std::string ("forward_model.toml"), error_message);
      }
      return std::make_unique<StarBlackBody>(spectral_grid);

    case stellar_modules::file :
      if (parameters.size() != 1)
      {
        std::string error_message =
          "Stellar model requires one parameter (the file path)!\n";
        throw InvalidInput(std::string ("forward_model.toml"), error_message);
      }
      return std::make_unique<StarSpectrumFile>(
          parameters[0],
          spectral_grid);

    case stellar_modules::grid :
      if (parameters.size() != 1)
      {
        std::string error_message =
          "Stellar grid model requires one parameter (the parameter file path)!\n";
        throw InvalidInput(std::string ("forward_model.toml"), error_message);
      }
      return std::make_unique<StellarSpectrumGrid>(
          parameters[0],
          spectral_grid);
  }


  return nullptr;
}


}
#endif
