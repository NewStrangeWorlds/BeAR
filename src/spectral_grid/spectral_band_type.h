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


#ifndef _spectral_band_type_h
#define _spectral_band_type_h

#include <string>
#include <algorithm>
#include <vector>

#include "../additional/exceptions.h"


namespace bear {


namespace band_type{
  enum class id {
    spectroscopy, 
    photometry, 
    band_spectroscopy}; 
  const std::vector<std::string> description {
    "spectroscopy", 
    "photometry", 
    "band-spectroscopy"};
  const std::vector<std::string> description_alt {
    "Spectroscopy", 
    "Photometry", 
    "Band-spectroscopy"};
}



inline band_type::id selectBandType(
  const std::string band_input,
  const std::string observation_name)
{
  auto it = std::find(
    band_type::description.begin(),
    band_type::description.end(),
    band_input);

  auto it_alt = std::find(
    band_type::description_alt.begin(),
    band_type::description_alt.end(),
    band_input);

  if (it == band_type::description.end() && it_alt == band_type::description_alt.end())
  {
    std::string error_message = 
      "Unsupported band type *" 
      + band_input 
      + "* for observation "
      + observation_name + "\n";
    throw InvalidInput(std::string ("selectBandType"), error_message);
  }
  
  if (it != band_type::description.end())
    return static_cast<band_type::id>(std::distance(band_type::description.begin(), it));
  else
    return static_cast<band_type::id>(std::distance(band_type::description_alt.begin(), it_alt));
}

}


#endif