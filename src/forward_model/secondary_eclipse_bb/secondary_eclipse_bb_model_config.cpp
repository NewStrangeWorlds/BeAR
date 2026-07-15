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


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "secondary_eclipse_bb.h"

#include "../../additional/exceptions.h"


namespace bear{


OccultationBlackBodyConfig::OccultationBlackBodyConfig (
  const std::string& folder_path,
  const std::string& file_name)
{
  readConfigFile(folder_path + file_name);
}


OccultationBlackBodyConfig::OccultationBlackBodyConfig (
  const std::string stellar_spectrum_model_,
  const std::vector<std::string>& stellar_model_parameters_)
{
  stellar_spectrum_model = stellar_spectrum_model_;
  stellar_model_parameters = stellar_model_parameters_;
}



void OccultationBlackBodyConfig::readConfigFile(const std::string& file_name)
{
  std::cout << "Parameters read from " << file_name << " :\n";

  toml::table cfg = parseConfigFile(file_name);

  readModelBlock(cfg, "stellar_spectrum",
    stellar_spectrum_model, stellar_model_parameters, "Stellar spectrum model");
}



}

