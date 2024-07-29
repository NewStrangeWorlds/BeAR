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
#include <iomanip>

#include "transmission.h"

#include "../../chemistry/chem_species.h"
#include "../atmosphere/atmosphere.h"


namespace bear{


TransmissionPostProcessConfig::TransmissionPostProcessConfig (const std::string& folder_path)
{
  const std::string config_file_name = folder_path + "post_process.config";

  readConfigFile(config_file_name);
}


void TransmissionPostProcessConfig::readConfigFile(const std::string& file_name)
{
  std::fstream file;
  file.open(file_name.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "\n Post-process config file not found. Using default options!\n\n";

    return;
  }

  std::cout << "\nParameters read from " << file_name << " :\n";

  std::string line;
  std::string input;

  std::getline(file, line);
  
  file >> input >> line;
  if (input == "Yes" || input == "yes" || input == "Y" || input == "y")
    delete_sampler_files = true;
  
  std::cout << "  - Delete sampler files: " << delete_sampler_files << "\n";


  std::getline(file, line);

  file >> input >> line;
  if (input == "Yes" || input == "yes" || input == "Y" || input == "y")
    save_spectra = true;
  else
    save_spectra = false;
  
  std::cout << "  - Save posterior spectra: " << save_spectra << "\n";


  std::getline(file, line);

  file >> input >> line;
  if (input == "Yes" || input == "yes" || input == "Y" || input == "y")
    save_temperatures = true;
  else
    save_temperatures = false;
  
  std::cout << "  - Save temperature structures: " << save_temperatures << "\n";


  std::getline(file, line);
  std::getline(file, line);
  std::string species;



  std::istringstream chem_input(line);
  std::vector<std::string> chemical_species;

  while (chem_input >> species)
    chemical_species.push_back(species);

  for (auto & i : chemical_species)
  {
    for (size_t j=0; j<constants::species_data.size(); ++j)
    {
      if (constants::species_data[j].symbol == i)
      {
        species_to_save.push_back(constants::species_data[j].id); 
        break;
      } 
    }
  }
   
  std::cout << "  - Save chemical species profiles: ";
  for (auto & i : species_to_save)
    std::cout << constants::species_data[i].symbol << " ";
  std::cout << "\n";

  file.close();
}



//calls the model specific posterior calculations
void TransmissionModel::postProcess(
  const std::vector< std::vector<double> >& model_parameter, 
  const size_t best_fit_model)
{
  TransmissionPostProcessConfig post_process_config(config->retrieval_folder_path);

  const size_t nb_models = model_parameter.size();
  std::vector< std::vector<double> > model_spectrum_bands;
  
  if (post_process_config.save_spectra)
    calcPostProcessSpectra(model_parameter, best_fit_model, model_spectrum_bands);


  std::vector<std::vector<double>> temperature_profiles(nb_models, std::vector<double>(nb_grid_points, 0));

  std::vector<std::vector<std::vector<double>>> mixing_ratios(
    nb_models, 
    std::vector<std::vector<double>>(constants::species_data.size(), 
    std::vector<double>(nb_grid_points,0)));


  for (size_t i=0; i<nb_models; ++i)
    postProcessModel(
      model_parameter[i], 
      model_spectrum_bands[i], 
      temperature_profiles[i], 
      mixing_ratios[i]);

  if (post_process_config.species_to_save.size() > 0)
    for (auto & i : post_process_config.species_to_save)
      savePostProcessChemistry(mixing_ratios, i);

  if (post_process_config.save_temperatures)
    savePostProcessTemperatures(temperature_profiles);
}



void TransmissionModel::postProcessModel(
  const std::vector<double>& model_parameter, 
  const std::vector<double>& model_spectrum_bands, 
  std::vector<double>& temperature_profile,
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



void TransmissionModel::savePostProcessChemistry(
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



void TransmissionModel::savePostProcessTemperatures(const std::vector<std::vector<double>>& temperature_profiles)
{
  //save the temperature profiles into a file
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


}

