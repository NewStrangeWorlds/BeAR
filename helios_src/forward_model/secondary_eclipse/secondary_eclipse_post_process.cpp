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


#include "secondary_eclipse.h"


#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>



#include "../../retrieval/retrieval.h"
#include "../../chemistry/chem_species.h"

#include "../../CUDA_kernels/data_management_kernels.h"
#include "../../CUDA_kernels/cross_section_kernels.h"
#include "../../CUDA_kernels/contribution_function_kernels.h"

#include "../atmosphere/atmosphere.h"


namespace helios{


//calls the model specific posterior calculations
void SecondaryEclipseModel::postProcess(const std::vector< std::vector<double> >& model_parameter, const std::vector< std::vector<double> >& model_spectrum_bands, const size_t best_fit_model)
{
  const size_t nb_models = model_parameter.size();

  //data structures for post process
  std::vector<double> effective_temperatures(nb_models, 0);
  std::vector<std::vector<double>> temperature_profiles(nb_models, std::vector<double>(nb_grid_points, 0));

  std::vector<chemical_species_id> postprocess_species {_H2O, _Na, _K, _TiO};
  std::vector<std::vector<std::vector<double>>> mixing_ratios(nb_models, std::vector<std::vector<double>>(constants::species_data.size(), std::vector<double>(nb_grid_points,0)));


  for (size_t i=0; i<nb_models; ++i)
  {
    postProcessModel(model_parameter[i], model_spectrum_bands[i], temperature_profiles[i], effective_temperatures[i], mixing_ratios[i]);
    
    if (i == 0)
    //if (i == best_fit_model)
      postProcessContributionFunctions();
  }


  for (auto & i : postprocess_species)
    savePostProcessChemistry(mixing_ratios, i);

  //savePostProcessEffectiveTemperatures(effective_temperatures);
  savePostProcessTemperatures(temperature_profiles);
}




void SecondaryEclipseModel::postProcessModel(const std::vector<double>& model_parameter, const std::vector<double>& model_spectrum_bands, 
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

  //effective_temperature = postProcessEffectiveTemperature(model_spectrum_bands);
}



void SecondaryEclipseModel::savePostProcessChemistry(const std::vector<std::vector<std::vector<double>>>& mixing_ratios, const unsigned int species)
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




void SecondaryEclipseModel::savePostProcessTemperatures(const std::vector<std::vector<double>>& temperature_profiles)
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





void SecondaryEclipseModel::savePostProcessEffectiveTemperatures(const std::vector<double>& effective_temperatures)
{
  //save the effective temperatures
  std::string file_name = retrieval->config->retrieval_folder_path + "/effective_temperatures.dat";
  
  std::fstream file(file_name.c_str(), std::ios::out);

  for (size_t i=0; i<effective_temperatures.size(); ++i)
    file << std::setprecision(10) << std::scientific << effective_temperatures[i] << "\n";
}




void SecondaryEclipseModel::postProcessContributionFunctions()
{
  size_t nb_spectral_points = retrieval->spectral_grid.nbSpectralPoints();


  initCrossSectionsHost(nb_spectral_points*nb_grid_points, absorption_coeff_gpu);

  for (size_t i=0; i<nb_grid_points; ++i)
    transport_coeff.calcTransportCoefficientsGPU(atmosphere.temperature[i], atmosphere.pressure[i], atmosphere.number_densities[i],
                                                 nb_grid_points, i,
                                                 absorption_coeff_gpu, nullptr);

  double* contribution_functions_dev = nullptr;

  //intialise the high-res spectrum on the GPU (set it to 0) 
  allocateOnDevice(contribution_functions_dev, nb_spectral_points*nb_grid_points);
  
  contributionFunctionGPU(contribution_functions_dev, absorption_coeff_gpu,
                          retrieval->spectral_grid.wavenumber_list_gpu,
                          atmosphere.temperature, atmosphere.altitude, nb_spectral_points);


  std::vector<double> contribution_functions_all(nb_spectral_points*nb_grid_points, 0.0);

  moveToHost(contribution_functions_dev, contribution_functions_all);
  deleteFromDevice(contribution_functions_dev);

  std::vector< std::vector<double> > contribution_functions(nb_grid_points, std::vector<double>(nb_spectral_points, 0));


  for (size_t i=0; i<nb_spectral_points; ++i)
    for (size_t j=0; j<nb_grid_points; ++j)
      contribution_functions[j][i] = contribution_functions_all[j*nb_spectral_points + i];


  for (size_t i=0; i<retrieval->observations.size(); ++i)
  {
    std::vector< std::vector<double> > contribution_functions_obs = contribution_functions;

    if (retrieval->observations[i].filter_response.size() != 0)
    {
      for (size_t j=0; j<nb_grid_points; ++j)
        contribution_functions_obs[j] = retrieval->observations[i].applyFilterResponseFunction(contribution_functions[j]);
    }

    std::vector< std::vector<double> > contribution_functions_bands(nb_grid_points, std::vector<double>(retrieval->observations[i].spectral_bands.nbBands()));

    for (size_t j=0; j<nb_grid_points; ++j)
       retrieval->observations[i].spectral_bands.bandIntegrateSpectrum(contribution_functions_obs[j], contribution_functions_bands[j]);

    saveContributionFunctions(contribution_functions_bands, i);
  }

}



void SecondaryEclipseModel::saveContributionFunctions(std::vector< std::vector<double>>& contribution_function, const size_t observation_index)
{
  std::string observation_name = retrieval->observations[observation_index].observationName();
  std::replace(observation_name.begin(), observation_name.end(), ' ', '_'); 
    
  std::string file_name = retrieval->config->retrieval_folder_path + "/contribution_function_" + observation_name + ".dat"; 
  

  std::fstream file(file_name.c_str(), std::ios::out);

  for (size_t j=0; j<nb_grid_points; ++j)
  {
    for (size_t i=0; i<retrieval->observations[observation_index].spectral_bands.nbBands(); ++i)
      file << std::setprecision(10) << std::scientific << contribution_function[j][i] << "\t";
     
    file << "\n";
  }

  file.close();
}



}

