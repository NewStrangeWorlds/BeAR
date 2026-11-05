// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "emission.h"

#include "../../additional/exceptions.h"


namespace bear{


EmissionModelConfig::EmissionModelConfig (
  const std::string& folder_path,
  const std::string& file_name)
{
  readConfigFile(folder_path + file_name);
}


EmissionModelConfig::EmissionModelConfig (
  const int nb_grid_points_,
  const double atmos_bottom_pressure_,
  const double atmos_top_pressure_,
  const std::string& temperature_profile_model_,
  const std::vector<std::string>& temperature_profile_parameters_,
  const std::string radiative_transfer_model_,
  const std::vector<std::string>& radiative_transfer_parameters_,
  const std::vector<std::string>& chemistry_model_,
  const std::vector<std::vector<std::string>>& chemistry_parameters_,
  const std::vector<std::string>& opacity_species_symbol_,
  const std::vector<std::string>& opacity_species_folder_)
  : EmissionModelConfig(
      nb_grid_points_,
      atmos_bottom_pressure_,
      atmos_top_pressure_,
      temperature_profile_model_,
      temperature_profile_parameters_,
      radiative_transfer_model_,
      radiative_transfer_parameters_,
      chemistry_model_,
      chemistry_parameters_,
      opacity_species_symbol_,
      opacity_species_folder_,
      std::vector<std::string>(),
      std::vector<std::vector<std::string>>())
{

}


EmissionModelConfig::EmissionModelConfig (
  const int nb_grid_points_,
  const double atmos_bottom_pressure_,
  const double atmos_top_pressure_,
  const std::string& temperature_profile_model_,
  const std::vector<std::string>& temperature_profile_parameters_,
  const std::string radiative_transfer_model_,
  const std::vector<std::string>& radiative_transfer_parameters_,
  const std::vector<std::string>& chemistry_model_,
  const std::vector<std::vector<std::string>>& chemistry_parameters_,
  const std::vector<std::string>& opacity_species_symbol_,
  const std::vector<std::string>& opacity_species_folder_,
  const std::vector<std::string>& cloud_model_,
  const std::vector<std::vector<std::string>>& cloud_model_parameters_)
{
  nb_grid_points = nb_grid_points_;
  atmos_boundaries[0] = atmos_bottom_pressure_;
  atmos_boundaries[1] = atmos_top_pressure_;
  temperature_profile_model = temperature_profile_model_;
  temperature_profile_parameters = temperature_profile_parameters_;
  radiative_transfer_model = radiative_transfer_model_;
  radiative_transfer_parameters = radiative_transfer_parameters_;
  chemistry_model = chemistry_model_;
  chemistry_parameters = chemistry_parameters_;
  opacity_species_symbol = opacity_species_symbol_;
  opacity_species_folder = opacity_species_folder_;

  cloud_model = cloud_model_;
  cloud_model_parameters = cloud_model_parameters_;
}




void EmissionModelConfig::readConfigFile(const std::string& file_name)
{
  std::cout << "Parameters read from " << file_name << " :\n";

  toml::table cfg = parseConfigFile(file_name);

  std::vector<double> pressure_boundaries;

  readAtmosphereConfig(cfg, nb_grid_points, pressure_boundaries);
  atmos_boundaries[0] = pressure_boundaries[0];
  atmos_boundaries[1] = pressure_boundaries[1];

  readModelBlock(cfg, "temperature",
    temperature_profile_model, temperature_profile_parameters, "Temperature profile");

  readModelList(cfg, "clouds",
    cloud_model, cloud_model_parameters, /*skip_none=*/true, /*required=*/false, "Cloud model");

  readModelBlock(cfg, "radiative_transfer",
    radiative_transfer_model, radiative_transfer_parameters, "Radiative transfer model");

  readModelList(cfg, "chemistry",
    chemistry_model, chemistry_parameters, /*skip_none=*/false, /*required=*/true, "Chemistry model");

  readOpacityConfig(cfg, "opacity",
    opacity_species_symbol, opacity_species_folder, /*required=*/true);
}


}

