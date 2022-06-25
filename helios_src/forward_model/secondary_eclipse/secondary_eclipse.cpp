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
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>

#include "secondary_eclipse.h"

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../retrieval/priors.h"
#include "../../observations/observations.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"
#include "../atmosphere/atmosphere.h"

#include "../stellar_spectrum.h"

namespace helios{


SecondaryEclipseModel::SecondaryEclipseModel (
  const SecondaryEclipseConfig model_config,
  Priors* priors_,
  GlobalConfig* config_,
  SpectralGrid* spectral_grid_,
  std::vector<Observation>& observations_) 
    : config(config_)
    , spectral_grid(spectral_grid_)
    , atmosphere(
        model_config.nb_grid_points,
        model_config.atmos_boundaries,
        config->use_gpu)
    , opacity_calc(
        config,
        spectral_grid,
        &atmosphere,
        model_config.opacity_species_symbol,
        model_config.opacity_species_folder,
        config->use_gpu,
        model_config.cloud_model != "none")
    , observations(observations_)
    , stellar_spectrum(
        config->retrieval_folder_path + model_config.stellar_spectrum_file,
        config->use_gpu,
        spectral_grid,
        observations)
{
  nb_grid_points = model_config.nb_grid_points;

  std::cout << "Forward model selected: Secondary Eclipse\n\n";

  //this forward model has three free general parameters
  nb_general_param = 3;

  for (auto & i : observations)
    nb_observation_points += i.nbPoints();

  initModules(model_config);

  setPriors(priors_);
}



//determines the basic atmospheric structure (temperature profile, chemistry...) from the free parameters supplied by MultiNest
bool SecondaryEclipseModel::calcAtmosphereStructure(const std::vector<double>& parameter)
{
  const double surface_gravity = std::pow(10,parameter[0]);
 
  bool neglect_model = false;


  //parameters for temperature profile and chemistry
  std::vector<double> temp_parameters(parameter.begin() + nb_general_param + nb_total_chemistry_param, 
                                      parameter.begin() + nb_general_param + nb_total_chemistry_param + temperature_profile->nbParameters());

  std::vector<double> chem_parameters (parameter.begin() + nb_general_param, 
                                       parameter.begin() + nb_general_param + nb_total_chemistry_param);
  

  //determine atmosphere structure
  neglect_model = atmosphere.calcAtmosphereStructure(surface_gravity, temperature_profile, temp_parameters, chemistry, chem_parameters);


  return neglect_model;
}



//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool SecondaryEclipseModel::calcModel(
  const std::vector<double>& parameter, 
  std::vector<double>& spectrum, 
  std::vector<double>& model_spectrum_bands)
{
  bool neglect = calcAtmosphereStructure(parameter);

  std::vector<double> cloud_parameters(
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param,
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_cloud_param);

  std::vector<CloudModel*> cm = {cloud_model};

  opacity_calc.calculate(cm, cloud_parameters);


  spectrum.assign(spectral_grid->nbSpectralPoints(), 0.0);
 
  radiative_transfer->calcSpectrum(
    atmosphere,
    opacity_calc.absorption_coeff, 
    opacity_calc.absorption_coeff, 
    opacity_calc.cloud_optical_depths, 
    opacity_calc.cloud_single_scattering, 
    opacity_calc.cloud_asym_param,
    1.0, 
    spectrum);


  //post-process the planet's high-res emission spectrum and bin it to the observational bands
  std::vector<double> planet_spectrum_bands(nb_observation_points, 0);
  postProcessSpectrum(spectrum, planet_spectrum_bands);


  model_spectrum_bands.assign(nb_observation_points, 0);
  std::vector<double> albedo_contribution(nb_observation_points, 0.0);
  const double radius_ratio = parameter[1];

  //calculate the secondary-eclipse depth with an optional geometric albedo contribution
  for (size_t i=0; i<model_spectrum_bands.size(); ++i)
    model_spectrum_bands[i] = 
      planet_spectrum_bands[i]/stellar_spectrum.flux_bands[i] 
      * radius_ratio*radius_ratio * 1e6 + albedo_contribution[i]*1e6;

  //convert the high-res spectrum to an eclipse depth as well
  for (size_t i=0; i<spectral_grid->nbSpectralPoints(); ++i)
    spectrum[i] = spectrum[i]/stellar_spectrum.flux[i] * radius_ratio*radius_ratio * 1e6;

  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool SecondaryEclipseModel::calcModelGPU(
  const std::vector<double>& parameter, 
  double* model_spectrum_gpu, 
  double* model_spectrum_bands)
{
  const double radius_ratio = parameter[1];

  bool neglect = calcAtmosphereStructure(parameter);


  std::vector<double> cloud_parameters(
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param,
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_cloud_param);

  std::vector<CloudModel*> cm = {cloud_model};

  opacity_calc.calculateGPU(cm, cloud_parameters);


  radiative_transfer->calcSpectrumGPU(
    atmosphere,
    opacity_calc.absorption_coeff_gpu, 
    opacity_calc.scattering_coeff_dev, 
    opacity_calc.cloud_optical_depths_dev,
    opacity_calc.cloud_single_scattering_dev,
    opacity_calc.cloud_asym_param_dev,
    1.0,
    model_spectrum_gpu);


  //post-process the planet's high-res emission spectrum and bin it to the observational bands
  double* planet_spectrum_bands = nullptr;
  allocateOnDevice(planet_spectrum_bands, nb_observation_points);

  postProcessSpectrumGPU(model_spectrum_gpu, planet_spectrum_bands);


  //std::vector<double> albedo_contribution(retrieval->nb_total_bands, geometric_albedo*radius_distance_ratio*radius_distance_ratio);
  std::vector<double> albedo_contribution(nb_observation_points, 0.0);

  double* albedo_contribution_gpu = nullptr;
  //moveToDevice(albedo_contribution_gpu, albedo_contribution);


  calcSecondaryEclipseGPU(
    model_spectrum_bands, 
    planet_spectrum_bands, 
    stellar_spectrum.flux_bands_dev, 
    nb_observation_points,
    radius_ratio, 
    albedo_contribution_gpu);

  deleteFromDevice(albedo_contribution_gpu);
  deleteFromDevice(planet_spectrum_bands);

  //convert the original high-res spectrum also to a secondary eclipse
  calcSecondaryEclipseGPU(
    model_spectrum_gpu, 
    model_spectrum_gpu, 
    stellar_spectrum.flux_dev, 
    spectral_grid->nbSpectralPoints(),
    radius_ratio, 
    albedo_contribution_gpu);

  return neglect;
}



std::vector<double> SecondaryEclipseModel::calcSecondaryEclipse(
  std::vector<double>& planet_spectrum_bands, 
  const double radius_ratio,
  const double geometric_albedo, 
  const double radius_distance_ratio)
{
  std::vector<double> secondary_eclipse(planet_spectrum_bands.size(), 0.0);

  for (size_t i=0; i<secondary_eclipse.size(); ++i)
    secondary_eclipse[i] = planet_spectrum_bands[i]/stellar_spectrum.flux[i] * radius_ratio*radius_ratio
                           + geometric_albedo * radius_distance_ratio*radius_distance_ratio;
  
  return secondary_eclipse;
}




//integrate the high-res spectrum to observational bands
//and convolve if necessary 
void SecondaryEclipseModel::postProcessSpectrum(
  std::vector<double>& model_spectrum, std::vector<double>& model_spectrum_bands)
{
  model_spectrum_bands.assign(nb_observation_points, 0.0);

  std::vector<double>::iterator it = model_spectrum_bands.begin();

  for (size_t i=0; i<observations.size(); ++i)
  { 
    const bool is_flux = true;
    std::vector<double> observation_bands = observations[i].processModelSpectrum(model_spectrum, is_flux);

    //copy the band-integrated values for this observation into the global
    //vector of all band-integrated points, model_spectrum_bands
    std::copy(observation_bands.begin(), observation_bands.end(), it);
    it += observation_bands.size();
  }
}


//integrate the high-res spectrum to observational bands
//and convolve if necessary 
void SecondaryEclipseModel::postProcessSpectrumGPU(
  double* model_spectrum_gpu, 
  double* model_spectrum_bands)
{
  unsigned int start_index = 0;
  for (size_t i=0; i<observations.size(); ++i)
  {
    const bool is_flux = true;
    observations[i].processModelSpectrumGPU(
      model_spectrum_gpu, 
      model_spectrum_bands, 
      start_index, 
      is_flux);

    start_index += observations[i].spectral_bands.nbBands();
  }

}


SecondaryEclipseModel::~SecondaryEclipseModel()
{
  delete radiative_transfer;
  delete temperature_profile;
  delete cloud_model;
  
  for (auto & i : chemistry)
    delete i;
}



}

