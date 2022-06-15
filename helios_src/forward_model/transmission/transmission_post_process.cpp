/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
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


#include "transmission.h"


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>


#include "../../retrieval/retrieval.h"
#include "../../chemistry/chem_species.h"

#include "../atmosphere/atmosphere.h"


namespace helios{


//calls the model specific posterior calculations
void TransmissionModel::postProcess(const std::vector< std::vector<double> >& model_parameter, const std::vector< std::vector<double> >& model_spectrum_bands, const size_t best_fit_model)
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


  savePostProcessTemperatures(temperature_profiles);
}




void TransmissionModel::postProcessModel(const std::vector<double>& model_parameter, const std::vector<double>& model_spectrum_bands, 
                                       std::vector<double>& temperature_profile, double& effective_temperature,
                                       std::vector<std::vector<double>>& mixing_ratios)
{
  calcAtmosphereStructure(model_parameter);

  
  for (auto & i : constants::species_data)
  { 
    for (size_t j=0; j<nb_grid_points; ++j)
      mixing_ratios[i.id][j] = atmosphere.number_densities[j][i.id]/atmosphere.number_densities[j][_TOTAL];
  }


  temperature_profile = atmosphere.temperature;
}



void TransmissionModel::savePostProcessChemistry(const std::vector<std::vector<std::vector<double>>>& mixing_ratios, const unsigned int species)
{
  std::fstream file;
  std::string file_name = retrieval->config->retrieval_folder_path + "/chem_";
  
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




void TransmissionModel::savePostProcessTemperatures(const std::vector<std::vector<double>>& temperature_profiles)
{
  //save the temperature profiles into a file
  std::fstream file;
  std::string file_name = retrieval->config->retrieval_folder_path + "/temperature_structures.dat";
  file.open(file_name.c_str(), std::ios::out);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    file << std::setprecision(10) << std::scientific << atmosphere.pressure[i];

    for(size_t j=0; j<temperature_profiles.size(); ++j)
      file << "\t" << temperature_profiles[j][i];

    file << "\n";
  }

}



}

