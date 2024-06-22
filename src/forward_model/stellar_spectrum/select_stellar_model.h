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


#ifndef _select_stellar_model_h
#define _select_stellar_model_h

#include <vector>
#include <string>
#include <algorithm>

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



inline StellarSpectrumModel* selectStellarModel(
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
    throw InvalidInput(std::string ("forward_model.config"), error_message);
  }


  //get the id of the chosen module
  stellar_modules::id module_id = static_cast<stellar_modules::id>(
    std::distance(stellar_modules::description.begin(), it));


  //create the temperature profile object based on the chosen module
  StellarSpectrumModel* stellar_spectrum_model = nullptr;

  switch (module_id)
  {
    case stellar_modules::blackbody :
      if (parameters.size() != 0)
      {
        std::string error_message = 
          "Stellar blackbody model requires no parameters!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);
      }
      {
        StarBlackBody* star = new StarBlackBody(spectral_grid);
        stellar_spectrum_model = star;  
      }
      break;

    case stellar_modules::file :
      if (parameters.size() != 1)
      {
        std::string error_message = 
          "Stellar model requires one parameter (the file path)!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);
      }
      {
        StarSpectrumFile* star = new StarSpectrumFile(
          parameters[0],
          spectral_grid);
        stellar_spectrum_model = star;  
      }
      break;

    case stellar_modules::grid :
      if (parameters.size() != 1)
      {
        std::string error_message = 
          "Stellar grid model requires one parameter (the parameter file path)!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);
      }
      {
        StellarSpectrumGrid* star = new StellarSpectrumGrid(
          parameters[0],
          spectral_grid);
        stellar_spectrum_model = star;  
      }
      break;

  }


  return stellar_spectrum_model;
}


}
#endif

