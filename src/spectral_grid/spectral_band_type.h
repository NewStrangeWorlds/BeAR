// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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