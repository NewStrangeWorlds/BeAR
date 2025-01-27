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

#include "transmission.h"

#include "../../additional/exceptions.h"


namespace bear{


TransmissionModelConfig::TransmissionModelConfig (const std::string& folder_path)
{
  const std::string config_file_name = folder_path + "forward_model.config";

  readConfigFile(config_file_name);
}



void TransmissionModelConfig::readConfigFile(const std::string& file_name)
{
  std::fstream file;
  file.open(file_name.c_str(), std::ios::in);

  if (file.fail())
    throw FileNotFound(std::string ("TransmissionModelConfig::readConfigFile"), file_name);

  std::cout << "Parameters read from " << file_name << " :\n";


  std::vector<double> pressure_boundaries;
  
  readAtmosphereConfig(file, nb_grid_points, pressure_boundaries);
  atmos_boundaries[0] = pressure_boundaries[0];
  atmos_boundaries[1] = pressure_boundaries[1];

  std::string fit_mode = readParameter(
    file,
    std::string("Fit for mean molecular weight or scale height"),
    std::vector<std::string>({"mmw", "sh", "no", "No"}));
  
  if (fit_mode == "mmw")
  {
    fit_mean_molecular_weight = true;
    std::cout << "- Fit for mean molecular weight: yes\n";
  }

  if (fit_mode == "sh")
  {
    fit_scale_height = true;
    std::cout << "- Fit for scale height: yes\n";
  }
  
  std::string var_gravity = readParameter(
    file,
    std::string("Use variable gravity"),
    std::vector<std::string>({"Yes", "yes", "no", "No"}));
  
  if (var_gravity == "Yes" || var_gravity == "yes")
  {
    use_variable_gravity = true;
    std::cout << "- Use variable gravity: yes\n";
  }
  else
  {
    use_variable_gravity = false;
    std::cout << "- Use variable gravity: no\n";
  }

  readTemperatureConfig(file, temperature_profile_model, temperature_profile_parameters);

  readCloudConfig(file, cloud_model, cloud_model_parameters);

  if (cloud_model.size() == 0) 
    use_cloud_model = false;
  else
    use_cloud_model = true;

  readModuleConfig(file, modules, modules_parameters);

  if (modules.size() == 0) 
    use_optional_modules = false;
  else
    use_optional_modules = true;
  
  readChemistryConfig(file, chemistry_model, chemistry_parameters);
  
  readOpacityConfig(file, opacity_species_symbol, opacity_species_folder);


  file.close();
}




}

