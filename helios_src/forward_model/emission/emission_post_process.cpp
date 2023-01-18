/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2022 Daniel Kitzmann
*
* Helios-r2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Helios-r2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* Helios-r2 directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>

#include "emission.h"

#include "../../retrieval/retrieval.h"
#include "../../chemistry/chem_species.h"
#include "../atmosphere/atmosphere.h"


namespace helios{


//calls the model specific posterior calculations
void EmissionModel::postProcess(
  const std::vector< std::vector<double> >& model_parameter, 
  const std::vector< std::vector<double> >& model_spectrum_bands, 
  const size_t best_fit_model)
{
  const size_t nb_models = model_parameter.size();

  //data structures for post process
  std::vector<double> effective_temperatures(nb_models, 0);
  std::vector<std::vector<double>> temperature_profiles(nb_models, std::vector<double>(nb_grid_points, 0));

  std::vector<chemical_species_id> postprocess_species {_H2O, _K, _NH3, _CH4};
  std::vector<std::vector<std::vector<double>>> mixing_ratios(nb_models, std::vector<std::vector<double>>(constants::species_data.size(), std::vector<double>(nb_grid_points,0)));


  for (size_t i=0; i<nb_models; ++i)
    postProcessModel(model_parameter[i], model_spectrum_bands[i], temperature_profiles[i], effective_temperatures[i], mixing_ratios[i]);


  for (auto & i : postprocess_species)
    savePostProcessChemistry(mixing_ratios, i);


  savePostProcessEffectiveTemperatures(effective_temperatures);
  savePostProcessTemperatures(temperature_profiles);
}




void EmissionModel::postProcessModel(
  const std::vector<double>& model_parameter, 
  const std::vector<double>& model_spectrum_bands, 
  std::vector<double>& temperature_profile, 
  double& effective_temperature,
  std::vector<std::vector<double>>& mixing_ratios)
{
  calcAtmosphereStructure(model_parameter);

  
  for (auto & i : constants::species_data)
  { 
    for (size_t j=0; j<nb_grid_points; ++j)
      mixing_ratios[i.id][j] = atmosphere.number_densities[j][i.id]/atmosphere.number_densities[j][_TOTAL];
  }


  temperature_profile = atmosphere.temperature;

  const double radius_distance_scaling = radiusDistanceScaling(model_parameter);

  effective_temperature = postProcessEffectiveTemperature(model_spectrum_bands, radius_distance_scaling);
}



void EmissionModel::savePostProcessChemistry(
  const std::vector<std::vector<std::vector<double>>>& mixing_ratios, 
  const unsigned int species)
{
  std::fstream file;
  std::string file_name = config->retrieval_folder_path + "/chem_";
  
  file_name += constants::species_data[species].symbol;
  file_name += ".dat";

  file.open(file_name.c_str(), std::ios::out);

  
  const size_t nb_models = mixing_ratios.size();

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    file << std::setprecision(10) << std::scientific << atmosphere.pressure[i];

    for (size_t j=0; j<nb_models; ++j)
      file << "\t" << mixing_ratios[j][species][i];

    file << "\n";
  }

}




void EmissionModel::savePostProcessTemperatures(
  const std::vector<std::vector<double>>& temperature_profiles)
{
  std::fstream file;
  std::string file_name = config->retrieval_folder_path + "/temperature_structures.dat";
  file.open(file_name.c_str(), std::ios::out);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    file << std::setprecision(10) << std::scientific << atmosphere.pressure[i];

    for(size_t j=0; j<temperature_profiles.size(); ++j)
      file << "\t" << temperature_profiles[j][i];

    file << "\n";
  }

}





void EmissionModel::savePostProcessEffectiveTemperatures(
  const std::vector<double>& effective_temperatures)
{
  std::string file_name = config->retrieval_folder_path + "/effective_temperatures.dat";
  
  std::fstream file(file_name.c_str(), std::ios::out);

  for (size_t i=0; i<effective_temperatures.size(); ++i)
    file << std::setprecision(10) << std::scientific << effective_temperatures[i] << "\n";
}



std::vector<double> EmissionModel::convertSpectrumToModel(const std::vector<double>& spectrum)
{
  std::vector<double> model_spectrum = spectrum;
  
  //convert from W m-2 cm to W m-2 micron-1
  for (size_t i=0; i<spectral_grid->nbSpectralPoints(); ++i)
    model_spectrum[i] = model_spectrum[i]/spectral_grid->wavelength_list[i]/spectral_grid->wavelength_list[i]*10000.0;

  return model_spectrum;
}


}

