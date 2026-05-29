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

#include "phase_curve.h"

#include "../../additional/exceptions.h"


namespace bear{


PhaseCurveConfig::PhaseCurveConfig (const std::string& folder_path)
{
  const std::string config_file_name = folder_path + "forward_model.config";

  readConfigFile(config_file_name);
}


PhaseCurveConfig::PhaseCurveConfig (
  const int nb_grid_points_,
  const double atmos_bottom_pressure_,
  const double atmos_top_pressure_,
  const std::string temperature_profile_model_,
  const std::vector<std::string>& temperature_profile_parameters_,
  const std::string radiative_transfer_model_,
  const std::vector<std::string>& radiative_transfer_parameters_,
  const std::vector<std::string>& chemistry_model_,
  const std::vector<std::vector<std::string>>& chemistry_parameters_,
  const std::vector<std::string>& opacity_species_symbol_,
  const std::vector<std::string>& opacity_species_folder_)
  : PhaseCurveConfig(
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



PhaseCurveConfig::PhaseCurveConfig (
  const int nb_grid_points_,
  const double atmos_bottom_pressure_,
  const double atmos_top_pressure_,
  const std::string temperature_profile_model_,
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

  radiative_transfer_model = radiative_transfer_model_;
  radiative_transfer_parameters = radiative_transfer_parameters_;

  temperature_profile_model = temperature_profile_model_;
  temperature_profile_parameters = temperature_profile_parameters_;

  chemistry_model = chemistry_model_;
  chemistry_parameters = chemistry_parameters_;

  opacity_species_symbol = opacity_species_symbol_;
  opacity_species_folder = opacity_species_folder_;

  cloud_model = cloud_model_;
  cloud_model_parameters = cloud_model_parameters_;
}


void PhaseCurveConfig::readConfigFile(const std::string& file_name)
{
  std::fstream file;
  file.open(file_name.c_str(), std::ios::in);

  if (file.fail())
    throw FileNotFound(std::string ("PhaseCurveConfig::readConfigFile"), file_name);

  std::string line;
  std::string input;

  std::vector<double> pressure_boundaries;

  readAtmosphereConfig(file, nb_grid_points, pressure_boundaries);
  atmos_boundaries[0] = pressure_boundaries[0];
  atmos_boundaries[1] = pressure_boundaries[1];


  readTemperatureConfig(file, temperature_profile_model, temperature_profile_parameters);

  // Stellar spectrum model
  std::getline(file, line);
  std::getline(file, line);
  {
    std::istringstream stellar_input(line);
    stellar_input >> stellar_spectrum_model;
    while (stellar_input >> input)
      stellar_model_parameters.push_back(input);
    std::cout << "- Stellar spectrum model: " << stellar_spectrum_model;
    for (auto& p : stellar_model_parameters) std::cout << "  " << p;
    std::cout << "\n";
  }
  // Optional: highres_stellar_smooth_sigma <value>
  std::getline(file, line);
  {
    std::istringstream check(line);
    std::string keyword;
    if (check >> keyword && keyword == "highres_stellar_smooth_sigma")
    {
      check >> highres_stellar_smooth_sigma;
      std::cout << "- High-res stellar smoothing sigma: "
                << highres_stellar_smooth_sigma << " px\n";
      std::getline(file, line);  // read blank separator
    }
    // if not the keyword, treat line as the blank separator (no action needed)
  }

  readCloudConfig(file, cloud_model, cloud_model_parameters);

  std::getline(file, line);
  std::getline(file, line);

  std::istringstream input_stream(line);

  input_stream >> radiative_transfer_model;

  while (input_stream >> input)
    radiative_transfer_parameters.push_back(input);

  std::cout << "- Radiative transfer model: " << radiative_transfer_model;
  for (auto & i : radiative_transfer_parameters) std::cout << "  " << i;
  std::cout << "\n";

  std::getline(file, line);

  readModuleConfig(file, modules, modules_parameters);

  readChemistryConfig(file, chemistry_model, chemistry_parameters);

  readOpacityConfig(file, opacity_species_symbol, opacity_species_folder);

  // Optional high-res opacity section. readOpacityConfig consumes the section
  // header itself, so call it directly; if we're at EOF it loads nothing.
  readOpacityConfig(file, opacity_species_symbol_highres, opacity_species_folder_highres);

  if (opacity_species_symbol_highres.empty())
  {
    opacity_species_symbol_highres = opacity_species_symbol;
    opacity_species_folder_highres = opacity_species_folder;
  }

  file.close();
}


}
