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

#include "emission.h"

#include "../../additional/exceptions.h"


namespace bear{


EmissionModelConfig::EmissionModelConfig (const std::string& folder_path)
{
  const std::string config_file_name = folder_path + "forward_model.config";

  readConfigFile(config_file_name);
}



void EmissionModelConfig::readConfigFile(const std::string& file_name)
{
  std::fstream file;
  file.open(file_name.c_str(), std::ios::in);

  if (file.fail())
    throw FileNotFound(std::string ("EmissionModelConfig::readConfigFile"), file_name);

  std::string line;
  std::string input;

  std::vector<double> pressure_boundaries;
  
  readAtmosphereConfig(file, nb_grid_points, pressure_boundaries);
  atmos_boundaries[0] = pressure_boundaries[0];
  atmos_boundaries[1] = pressure_boundaries[1];

  readTemperatureConfig(file, temperature_profile_model, temperature_profile_parameters);

  readCloudConfig(file, cloud_model, cloud_model_parameters);

  if (cloud_model.size() == 0) 
    use_cloud_model = false;
  else
    use_cloud_model = true;

  //the radiative transfer input
  std::getline(file, line);
  std::getline(file, line);
  
  std::istringstream input_stream(line);
  input_stream.str(line); input_stream.clear();

  input_stream >> radiative_transfer_model;

  while (input_stream >> input)
    radiative_transfer_parameters.push_back(input);

  std::cout << "- Radiative transfer model: " << radiative_transfer_model;
  for (auto & i : radiative_transfer_parameters) std::cout << "  " << i;
  std::cout << "\n";

  std::getline(file, line);

  readChemistryConfig(file, chemistry_model, chemistry_parameters);
  
  readOpacityConfig(file, opacity_species_symbol, opacity_species_folder);

  file.close();
}


}

