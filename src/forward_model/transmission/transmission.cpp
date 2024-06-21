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
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>

#include "transmission.h"

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../retrieval/priors.h"
#include "../../observations/observations.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"
#include "../atmosphere/atmosphere.h"
#include "../../transport_coeff/opacity_calc.h"


namespace helios{


TransmissionModel::TransmissionModel (
  const TransmissionModelConfig model_config,
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
        model_config.use_cloud_model)
    , observations(observations_)
{
  nb_grid_points = model_config.nb_grid_points;
  
  std::cout << "Forward model selected: Transmission\n\n"; 

  //this forward model has three free general parameters
  nb_general_param = 3;

  if (model_config.fit_mean_molecular_weight || model_config.fit_scale_height)
    nb_general_param += 1;

  fit_mean_molecular_weight = model_config.fit_mean_molecular_weight;
  fit_scale_height = model_config.fit_scale_height;

  for (auto & i : observations)
  {
    nb_observation_points += i.nbPoints();
    nb_spectrum_modifier_param += i.nb_modifier_param;
  }

  //select and set up the modules
  initModules(model_config);

  setPriors(priors_);
}



bool TransmissionModel::calcAtmosphereStructure(const std::vector<double>& parameter)
{
  const double surface_gravity = std::pow(10,parameter[0]);

  bool neglect_model = false;

  //parameters for temperature profile and chemistry
  std::vector<double> temp_parameters(
    parameter.begin() + nb_general_param + nb_total_chemistry_param, 
    parameter.begin() + nb_general_param + nb_total_chemistry_param + temperature_profile->nbParameters());

  std::vector<double> chem_parameters (
    parameter.begin() + nb_general_param, 
    parameter.begin() + nb_general_param + nb_total_chemistry_param);


  if (!fit_mean_molecular_weight && !fit_scale_height)
  {
    neglect_model = atmosphere.calcAtmosphereStructure(
      surface_gravity, 
      temperature_profile, 
      temp_parameters, 
      chemistry, 
      chem_parameters);
  }
  else
  {
    const double param = parameter[3];

    if (fit_mean_molecular_weight)
      neglect_model = atmosphere.calcAtmosphereStructure(
        surface_gravity, 
        temperature_profile, 
        temp_parameters, 
        chemistry, 
        chem_parameters,
        param);
    else
      neglect_model = atmosphere.calcAtmosphereStructure(
        surface_gravity, 
        param,
        temperature_profile, 
        temp_parameters, 
        chemistry, 
        chem_parameters);
  }

  return neglect_model;
}



//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool TransmissionModel::calcModel(
  const std::vector<double>& parameter, 
  std::vector<double>& spectrum, 
  std::vector<double>& model_spectrum_bands)
{
  bool neglect = calcAtmosphereStructure(parameter);

  size_t nb_previous_param = nb_general_param + nb_total_chemistry_param + nb_temperature_param;

  std::vector<double> cloud_parameters(
      parameter.begin() + nb_previous_param,
      parameter.begin() + nb_previous_param + nb_total_cloud_param);

  opacity_calc.calculate(cloud_models, cloud_parameters);

  spectrum.assign(spectral_grid->nbSpectralPoints(), 0.0);

  const double bottom_radius = parameter[1] * constants::radius_earth;
  const double star_radius = parameter[2] * constants::radius_sun;
  
  cloud_extinction.assign(
      spectral_grid->nbSpectralPoints(), 
      std::vector<double>(nb_grid_points, 0.0));

  if (cloud_models.size() > 0)
  {
    cloud_models[0]->convertOpticalDepth(
      opacity_calc.cloud_optical_depths, 
      cloud_extinction, 
      atmosphere.altitude);
  }

  calcTransmissionSpectrum(bottom_radius, star_radius, spectrum);


  nb_previous_param = nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_total_cloud_param;

  for (auto & m : modules)
  {
    std::vector<double> module_parameters(
      parameter.begin() + nb_previous_param,
      parameter.begin() + nb_previous_param + m->nbParameters());
  
    m->modifySpectrum(module_parameters, &atmosphere, spectrum);
    
    nb_previous_param += m->nbParameters();
  }


  postProcessSpectrum(spectrum, model_spectrum_bands);


  nb_previous_param = nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_total_modules_param + nb_total_cloud_param;

  std::vector<double> modifier_parameters(
      parameter.begin() + nb_previous_param,
      parameter.begin() + nb_previous_param + nb_spectrum_modifier_param);

  auto param_it = modifier_parameters.begin();

  //apply spectrum modifier if necessary
  size_t start_index = 0;

  for (size_t i=0; i<observations.size(); ++i)
  {
    if (observations[i].nb_modifier_param != 0)
    {
      const double spectrum_modifier = *param_it;

      if (spectrum_modifier != 0)
        for (size_t j=start_index; j<start_index+observations[i].nbPoints(); ++j)
          model_spectrum_bands[j] += spectrum_modifier;

      param_it += observations[i].nb_modifier_param;
    }

    start_index += observations[i].spectral_bands.nbBands();
  }

  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool TransmissionModel::calcModelGPU(
  const std::vector<double>& parameter, 
  double* model_spectrum_gpu, 
  double* model_spectrum_bands)
{
  bool neglect = calcAtmosphereStructure(parameter);

  size_t nb_previous_param = nb_general_param + nb_total_chemistry_param + nb_temperature_param;

  std::vector<double> cloud_parameters(
      parameter.begin() + nb_previous_param,
      parameter.begin() + nb_previous_param + nb_total_cloud_param);


  opacity_calc.calculateGPU(cloud_models, cloud_parameters);


  if (cloud_extinction_gpu == nullptr)
    allocateOnDevice(cloud_extinction_gpu, nb_grid_points*spectral_grid->nbSpectralPoints());

  if (cloud_models.size() > 0)
  { 
    intializeOnDevice(cloud_extinction_gpu, nb_grid_points*spectral_grid->nbSpectralPoints());

    cloud_models[0]->convertOpticalDepthGPU(
      opacity_calc.cloud_optical_depths_dev,
      atmosphere.altitude_dev,
      nb_grid_points,
      spectral_grid->nbSpectralPoints(),
      cloud_extinction_gpu);
  }


  const double bottom_radius = parameter[1] * constants::radius_earth;
  const double star_radius = parameter[2] * constants::radius_sun;

  calcTransitDepthGPU(
    model_spectrum_gpu, 
    opacity_calc.absorption_coeff_gpu, 
    opacity_calc.scattering_coeff_dev, 
    cloud_extinction_gpu,
    atmosphere,
    spectral_grid->nbSpectralPoints(), 
    bottom_radius,
    star_radius);


  nb_previous_param = nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_total_cloud_param;

  for (auto & m : modules)
  {
    std::vector<double> module_parameters(
      parameter.begin() + nb_previous_param,
      parameter.begin() + nb_previous_param + m->nbParameters());
  
    m->modifySpectrumGPU(module_parameters, &atmosphere, model_spectrum_gpu);
    
    nb_previous_param += m->nbParameters();
  }


  postProcessSpectrumGPU(model_spectrum_gpu, model_spectrum_bands);


  nb_previous_param = nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_total_modules_param + nb_total_cloud_param;

  std::vector<double> modifier_parameters(
      parameter.begin() + nb_previous_param,
      parameter.begin() + nb_previous_param + nb_spectrum_modifier_param);

  auto param_it = modifier_parameters.begin();

  //apply spectrum modifier if necessary
  unsigned int start_index = 0;
  
  for (size_t i=0; i<observations.size(); ++i)
  {
    if (observations[i].nb_modifier_param != 0)
    {
      const double spectrum_modifier = *param_it;

      if (spectrum_modifier != 0)
        observations[i].addShiftToSpectrumGPU(
          model_spectrum_bands, 
          start_index, 
          spectrum_modifier);

      param_it += observations[i].nb_modifier_param;
    }

    start_index += observations[i].spectral_bands.nbBands();
  }

  return neglect;
}



void TransmissionModel::postProcessSpectrum(
  std::vector<double>& model_spectrum, 
  std::vector<double>& model_spectrum_bands)
{
  model_spectrum_bands.assign(nb_observation_points, 0.0);
  
  std::vector<double>::iterator it = model_spectrum_bands.begin();

  for (size_t i=0; i<observations.size(); ++i)
  {
    const bool is_flux = false;
    std::vector<double> observation_bands = observations[i].processModelSpectrum(model_spectrum, is_flux);

    //copy the band-integrated values for this observation into the global
    //vector of all band-integrated points, model_spectrum_bands
    std::copy(observation_bands.begin(), observation_bands.end(), it);
    it += observation_bands.size();
  }
}


void TransmissionModel::postProcessSpectrumGPU(
  double* model_spectrum_gpu, 
  double* model_spectrum_bands)
{
  unsigned int start_index = 0;
  for (size_t i=0; i<observations.size(); ++i)
  {
    const bool is_flux = false;
    observations[i].processModelSpectrumGPU(
      model_spectrum_gpu, 
      model_spectrum_bands, 
      start_index, 
      is_flux);

    start_index += observations[i].spectral_bands.nbBands();
  }
}


TransmissionModel::~TransmissionModel()
{
  delete temperature_profile;

  for (auto & i : cloud_models)
    delete i;
  
  for (auto & i : chemistry)
    delete i;

  if (cloud_extinction_gpu != nullptr)
    deleteFromDevice(cloud_extinction_gpu);
}


}

