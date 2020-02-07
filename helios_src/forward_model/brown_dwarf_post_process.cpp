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


#include "brown_dwarf.h"


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>


#include "../retrieval/retrieval.h"



namespace helios{


//calls the model specific posterior calculations
void BrownDwarfModel::postProcess(const std::vector< std::vector<double> >& model_parameter, const std::vector< std::vector<double> >& model_spectrum_bands)
{
  const size_t nb_models = model_parameter.size();

  //data structures for post process
  std::vector<double> effective_temperatures(nb_models, 0);
  std::vector<std::vector<double>> temperature_profiles(nb_models, std::vector<double>(nb_grid_points, 0));


  for (size_t i=0; i<nb_models; ++i)
    postProcessModel(model_parameter[i], model_spectrum_bands[i], temperature_profiles[i], effective_temperatures[i]);


  savePostProcessEffectiveTemperatures(effective_temperatures);
  savePostProcessTemperatures(temperature_profiles);
}




void BrownDwarfModel::postProcessModel(const std::vector<double>& model_parameter, const std::vector<double>& model_spectrum_bands, 
                                       std::vector<double>& temperature_profile, double& effective_temperature)
{
  calcAtmosphereStructure(model_parameter);

  temperature_profile = temperature;

  effective_temperature = postProcessEffectiveTemperature(model_spectrum_bands);
}




void BrownDwarfModel::savePostProcessTemperatures(const std::vector<std::vector<double>>& temperature_profiles)
{
  //save the temperature profiles into a file
  std::fstream file;
  std::string file_name = retrieval->config->retrieval_folder_path + "/temperature_structures.dat";
  file.open(file_name.c_str(), std::ios::out);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    file << std::setprecision(10) << std::scientific << pressure[i];

    for(size_t j=0; j<temperature_profiles.size(); ++j)
      file << "\t" << temperature_profiles[j][i];

    file << "\n";
  }

}





void BrownDwarfModel::savePostProcessEffectiveTemperatures(const std::vector<double>& effective_temperatures)
{
  //save the effective temperatures
  std::string file_name = retrieval->config->retrieval_folder_path + "/effective_temperatures.dat";
  
  std::fstream file(file_name.c_str(), std::ios::out);

  for (size_t i=0; i<effective_temperatures.size(); ++i)
    file << std::setprecision(10) << std::scientific << effective_temperatures[i] << "\n";
}



}

