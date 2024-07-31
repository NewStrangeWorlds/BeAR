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

#include "emission.h"

#include "../../retrieval/retrieval.h"
#include "../../chemistry/chem_species.h"
#include "../atmosphere/atmosphere.h"
#include "../../CUDA_kernels/data_management_kernels.h"
#include "../../CUDA_kernels/contribution_function_kernels.h"
#include "../../additional/quadrature.h"


namespace bear{


EmissionPostProcessConfig::EmissionPostProcessConfig (const std::string& folder_path)
{
  const std::string config_file_name = folder_path + "post_process.config";

  readConfigFile(config_file_name);
}



void EmissionPostProcessConfig::readConfigFile(const std::string& file_name)
{
  std::fstream file;
  file.open(file_name.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "\n Post-process config file not found. Using default options!\n\n";

    return;
  }

  std::cout << "\nParameters read from " << file_name << " :\n";

  delete_sampler_files = readBooleanParameter(file, "Delete sampler files");

  save_spectra = readBooleanParameter(file, "Save posterior spectra");

  save_temperatures = readBooleanParameter(file, "Save temperature structures");

  save_effective_temperatures = readBooleanParameter(file, "Save effective temperatures");

  save_contribution_functions = readBooleanParameter(file, "Save contribution functions");
  
  species_to_save = readChemicalSpecies(file, "Save chemical species profiles");

  file.close();
}



//calls the model specific posterior calculations
void EmissionModel::postProcess(
  const std::vector< std::vector<double> >& model_parameter, 
  const size_t best_fit_model,
  bool& delete_unused_files)
{
  EmissionPostProcessConfig post_process_config(config->retrieval_folder_path);

  if (post_process_config.delete_sampler_files)
    delete_unused_files = true;

  const size_t nb_models = model_parameter.size();
  std::vector<double> integrated_flux(nb_models, 0);
  
  if (post_process_config.save_spectra || post_process_config.save_effective_temperatures)
    calcPostProcessSpectra(
      model_parameter, 
      best_fit_model, 
      post_process_config.save_spectra,
      integrated_flux);
  
  //data structures for post process
  std::vector<double> effective_temperatures(nb_models, 0);
  std::vector<std::vector<double>> temperature_profiles(nb_models, std::vector<double>(nb_grid_points, 0));
  std::vector<std::vector<std::vector<double>>> mixing_ratios(
    nb_models, 
    std::vector<std::vector<double>>(constants::species_data.size(), 
    std::vector<double>(nb_grid_points,0)));
  
  for (size_t i=0; i<nb_models; ++i)
  {
    postProcessModel(
      model_parameter[i], 
      integrated_flux[i], 
      temperature_profiles[i], 
      effective_temperatures[i], 
      mixing_ratios[i]);

    if (i == best_fit_model && post_process_config.save_contribution_functions)
      postProcessContributionFunctions(model_parameter[i]);
  }

  if (post_process_config.save_temperatures)
    savePostProcessTemperatures(temperature_profiles);

  if (post_process_config.species_to_save.size() > 0)
    for (auto & i : post_process_config.species_to_save)
      savePostProcessChemistry(mixing_ratios, i);

  if (post_process_config.save_effective_temperatures)
    savePostProcessEffectiveTemperatures(effective_temperatures);
}


void EmissionModel::calcPostProcessSpectra(
  const std::vector< std::vector<double> >& model_parameter,
  const size_t best_fit_model,
  const bool save_spectra,
  std::vector<double>& integrated_flux)
{
  const size_t nb_models = model_parameter.size();
  
  std::vector<std::vector<double>> model_spectrum_bands;
  model_spectrum_bands.resize(nb_models);

  integrated_flux.assign(nb_models, 0);
  
  std::cout << "\n";

  for (size_t i=0; i<nb_models; ++i)
  {
    std::cout << "\rPostprocess spectra, model " << i << " of " << nb_models << std::flush;
    
    std::vector<double> model_spectrum_high_res;

    calcPostProcessSpectrum(
      model_parameter[i],
      model_spectrum_high_res,
      model_spectrum_bands[i]);

    integrated_flux[i] = aux::quadratureTrapezoidal(spectral_grid->wavenumber_list, model_spectrum_high_res);

    if (i == best_fit_model && save_spectra)
      saveBestFitSpectrum(model_spectrum_high_res);
  }

  std::cout << "\n";
  
  if (save_spectra)
    savePostProcessSpectra(model_spectrum_bands);
}



void EmissionModel::postProcessModel(
  const std::vector<double>& model_parameter, 
  const double integrated_flux, 
  std::vector<double>& temperature_profile, 
  double& effective_temperature,
  std::vector<std::vector<double>>& mixing_ratios)
{
  calcAtmosphereStructure(model_parameter);

  for (auto & i : constants::species_data)
  { 
    for (size_t j=0; j<nb_grid_points; ++j)
      mixing_ratios[i.id][j] = atmosphere.number_densities[j][i.id]
                              /atmosphere.number_densities[j][_TOTAL];
  }

  temperature_profile = atmosphere.temperature;

  const double radius_distance_scaling = radiusDistanceScaling(model_parameter);

  effective_temperature = postProcessEffectiveTemperature(integrated_flux, radius_distance_scaling);
}


double EmissionModel::postProcessEffectiveTemperature(
  const double integrated_flux, 
  const double radius_distance_scaling)
{
  const double total_flux_scaled = integrated_flux / radius_distance_scaling;
  
  //use the Stefan-Boltzmann law to convert the total flux to effective temperatures
  //total_flux has units of W m-2 and needs to be converted to cgs
  const double effective_temperature = std::pow(total_flux_scaled*1000.0/constants::stefan_boltzmann, 0.25);

  return effective_temperature;
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



void EmissionModel::postProcessContributionFunctions(
  const std::vector<double>& parameter)
{
  std::vector<double> cloud_parameters(
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param,
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_total_cloud_param);

  opacity_calc.calculateGPU(cloud_models, cloud_parameters);

  double* contribution_functions_dev = nullptr;
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  //intialise the high-res spectrum on the GPU (set it to 0) 
  allocateOnDevice(contribution_functions_dev, nb_spectral_points*nb_grid_points);
  
  contributionFunctionGPU(
    contribution_functions_dev, 
    opacity_calc.absorption_coeff_gpu,
    spectral_grid->wavenumber_list_gpu,
    atmosphere.temperature, 
    atmosphere.altitude,
    nb_spectral_points);

  std::vector<double> contribution_functions_all(nb_spectral_points*nb_grid_points, 0.0);

  moveToHost(contribution_functions_dev, contribution_functions_all);
  deleteFromDevice(contribution_functions_dev);

  std::vector< std::vector<double> > contribution_functions(nb_grid_points, std::vector<double>(nb_spectral_points, 0));


  for (size_t i=0; i<nb_spectral_points; ++i)
    for (size_t j=0; j<nb_grid_points; ++j)
      contribution_functions[j][i] = contribution_functions_all[j*nb_spectral_points + i];


  for (size_t i=0; i<observations.size(); ++i)
  {
    std::vector<std::vector<double>> contribution_functions_band(
      nb_grid_points,
      std::vector<double>(observations[i].nbPoints(), 0));

    const bool is_flux = false;
      
    for (size_t j=0; j<nb_grid_points; ++j)
      contribution_functions_band[j] = 
        observations[i].processModelSpectrum(contribution_functions[j], is_flux);

    saveContributionFunctions(contribution_functions_band, i);
  }
}



void EmissionModel::saveContributionFunctions(
  std::vector< std::vector<double>>& contribution_function, const size_t observation_index)
{
  std::string observation_name = observations[observation_index].observationName();
  std::replace(observation_name.begin(), observation_name.end(), ' ', '_'); 
    
  std::string file_name = config->retrieval_folder_path + "/contribution_function_" + observation_name + ".dat"; 
  

  std::fstream file(file_name.c_str(), std::ios::out);

  for (size_t j=0; j<nb_grid_points; ++j)
  { 
    file << std::setprecision(10) << std::scientific << atmosphere.pressure[j] << "\t";

    for (size_t i=0; i<observations[observation_index].spectral_bands.nbBands(); ++i)
      file << std::setprecision(10) << std::scientific << contribution_function[j][i] << "\t";
     
    file << "\n";
  }

  file.close();
}


//overload of the generic forward model function
//changes units of the computed spectrum from W m-2 cm to W m-2 micron-1
void EmissionModel::saveBestFitSpectrum(const std::vector<double>& spectrum)
{ 
  std::string file_name = config->retrieval_folder_path + "/spectrum_best_fit_hr.dat";

  std::fstream file(file_name.c_str(), std::ios::out);

  for (size_t i=0; i<spectrum.size(); ++i)
    file << std::setprecision(10) << std::scientific
         << spectral_grid->wavelength_list[i] << "\t"
         << spectrum[i]/spectral_grid->wavelength_list[i]/spectral_grid->wavelength_list[i]*10000.0 << "\n";

  file.close();
}


}

