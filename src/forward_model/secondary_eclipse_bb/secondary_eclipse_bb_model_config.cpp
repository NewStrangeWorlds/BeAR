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


OccultationBlackBodyConfig::OccultationBlackBodyConfig (const std::string& folder_path)
{
  const std::string config_file_name = folder_path + "forward_model.config";

  readConfigFile(config_file_name);
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
  std::fstream file;
  file.open(file_name.c_str(), std::ios::in);

  
  if (file.fail())  
    throw FileNotFound(std::string ("OccultationBlackBodyConfig::readConfigFile"), file_name);

  
  std::string line;
  std::string input;

  //the stellar spectrum model
  std::getline(file, line);
  std::getline(file, line);
  std::istringstream stellar_input(line);

  stellar_input >> stellar_spectrum_model;

  while (stellar_input >> input)
    stellar_model_parameters.push_back(input);

  std::cout << "- Stellar spectrum model: " << stellar_spectrum_model;
  for (auto & i : stellar_model_parameters) std::cout << "  " << i;
  std::cout << "\n";

  file.close();
}



}

