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


TransmissionModelConfig::TransmissionModelConfig(
  const std::string& folder_path,
  const std::string& file_name)
{
  readConfigFile(folder_path + file_name);
}



TransmissionModelConfig::TransmissionModelConfig(
  const int nb_grid_points_,
  const double atmos_bottom_pressure_,
  const double atmos_top_pressure_,
  const std::string& temperature_profile_model_,
  const std::vector<std::string>& temperature_profile_parameters_,
  const std::vector<std::string>& chemistry_model_,
  const std::vector<std::vector<std::string>>& chemistry_parameters_,
  const std::vector<std::string>& opacity_species_symbol_,
  const std::vector<std::string>& opacity_species_folder_)
  : TransmissionModelConfig(
      false,
      false,
      false,
      nb_grid_points_,
      atmos_bottom_pressure_,
      atmos_top_pressure_,
      temperature_profile_model_,
      temperature_profile_parameters_,
      chemistry_model_,
      chemistry_parameters_,
      opacity_species_symbol_,
      opacity_species_folder_,
      std::vector<std::string>(),
      std::vector<std::vector<std::string>>(),
      std::vector<std::string>(),
      std::vector<std::vector<std::string>>())
{
  
}


TransmissionModelConfig::TransmissionModelConfig(
  const bool fit_mean_molecular_weight_, 
  const bool fit_scale_height_, 
  const bool use_variable_gravity_,
  const int nb_grid_points_,
  const double atmos_bottom_pressure_,
  const double atmos_top_pressure_,
  const std::string& temperature_profile_model_,
  const std::vector<std::string>& temperature_profile_parameters_,
  const std::vector<std::string>& chemistry_model_,
  const std::vector<std::vector<std::string>>& chemistry_parameters_,
  const std::vector<std::string>& opacity_species_symbol_,
  const std::vector<std::string>& opacity_species_folder_,
  const std::vector<std::string>& cloud_model_,
  const std::vector<std::vector<std::string>>& cloud_model_parameters_,
  const std::vector<std::string>& modules_,
  const std::vector<std::vector<std::string>>& modules_parameters_)
{
  fit_mean_molecular_weight = fit_mean_molecular_weight_;
  fit_scale_height = fit_scale_height_;
  use_variable_gravity = use_variable_gravity_;
  nb_grid_points = nb_grid_points_;
  atmos_boundaries[0] = atmos_bottom_pressure_;
  atmos_boundaries[1] = atmos_top_pressure_;
  temperature_profile_model = temperature_profile_model_;
  temperature_profile_parameters = temperature_profile_parameters_;
  chemistry_model = chemistry_model_;
  chemistry_parameters = chemistry_parameters_;
  cloud_model = cloud_model_;
  cloud_model_parameters = cloud_model_parameters_;
  modules = modules_;
  modules_parameters = modules_parameters_;
  opacity_species_symbol = opacity_species_symbol_;
  opacity_species_folder = opacity_species_folder_;
}



void TransmissionModelConfig::readConfigFile(const std::string& file_name)
{
  std::cout << "Parameters read from " << file_name << " :\n";

  toml::table cfg = parseConfigFile(file_name);

  std::vector<double> pressure_boundaries;

  readAtmosphereConfig(cfg, nb_grid_points, pressure_boundaries);
  atmos_boundaries[0] = pressure_boundaries[0];
  atmos_boundaries[1] = pressure_boundaries[1];

  const std::string fit_mode = readParameter(cfg, "fit_mode",
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

  use_variable_gravity = readBooleanParameter(cfg, "use_variable_gravity", false);

  readModelBlock(cfg, "temperature",
    temperature_profile_model, temperature_profile_parameters, "Temperature profile");

  readModelList(cfg, "clouds",
    cloud_model, cloud_model_parameters, /*skip_none=*/true, /*required=*/false, "Cloud model");

  readModelList(cfg, "modules",
    modules, modules_parameters, /*skip_none=*/true, /*required=*/false, "Optional modules");

  readModelList(cfg, "chemistry",
    chemistry_model, chemistry_parameters, /*skip_none=*/false, /*required=*/true, "Chemistry model");

  readOpacityConfig(cfg, "opacity",
    opacity_species_symbol, opacity_species_folder, /*required=*/true);
}


}

