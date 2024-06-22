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

#include "emission.h"

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../retrieval/priors.h"
#include "../../observations/observations.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"
#include "../atmosphere/atmosphere.h"
#include "../../transport_coeff/opacity_calc.h"


namespace bear{


EmissionModel::EmissionModel (
  const EmissionModelConfig model_config,
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
  
  std::cout << "Forward model selected: Emission\n\n"; 

  //this forward model has three free general parameters
  nb_general_param = 3;

  for (auto & i : observations)
  {
    nb_observation_points += i.nbPoints();
    nb_spectrum_modifier_param += i.nb_modifier_param;
  }

  initModules(model_config);

  setPriors(priors_);
}


double EmissionModel::radiusDistanceScaling(const std::vector<double>& parameter)
{
  const double scaling_factor = parameter[1];
  const double distance = parameter[2] * constants::parsec;

   //we assume a fixed prior radius of 1 Rj
  const double prior_radius = constants::radius_jupiter; 
  
  double scaling = prior_radius/distance;
  scaling = scaling*scaling * scaling_factor;

  return scaling;
}



bool EmissionModel::calcAtmosphereStructure(const std::vector<double>& parameter)
{
  const double surface_gravity = std::pow(10,parameter[0]);
  const double scaling_factor = parameter[1];

  //derived radius in Jupiter radii assuming that the radius prior is 1 Rj
  const double derived_radius = std::sqrt(scaling_factor);  

  //derived mass in Jupiter masses
  const double derived_mass = surface_gravity 
                            * std::pow(derived_radius*constants::radius_jupiter, 2) 
                            / constants::gravitation_const / constants::mass_jupiter;


  bool neglect_model = false;

  //if derived mass is larger than 80 Jupiter masses, 
  //we tell MultiNest to neglect this parameter combination
  if (derived_mass > 80) neglect_model = true;


  //parameters for temperature profile and chemistry
  std::vector<double> temp_parameters(
    parameter.begin() + nb_general_param + nb_total_chemistry_param, 
    parameter.begin() + nb_general_param + nb_total_chemistry_param + temperature_profile->nbParameters());

  std::vector<double> chem_parameters (
    parameter.begin() + nb_general_param, 
    parameter.begin() + nb_general_param + nb_total_chemistry_param);


  neglect_model = atmosphere.calcAtmosphereStructure(
    surface_gravity, 
    temperature_profile, 
    temp_parameters, 
    chemistry, 
    chem_parameters);


  return neglect_model;
}



//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool EmissionModel::calcModel(
  const std::vector<double>& parameter, 
  std::vector<double>& spectrum, 
  std::vector<double>& model_spectrum_bands)
{
  bool neglect = calcAtmosphereStructure(parameter);

  std::vector<double> cloud_parameters(
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param,
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_total_cloud_param);

  opacity_calc.calculate(cloud_models, cloud_parameters);


  spectrum.assign(spectral_grid->nbSpectralPoints(), 0.0);
  
  const double radius_distance_scaling = radiusDistanceScaling(parameter);

  radiative_transfer->calcSpectrum(
    atmosphere,
    opacity_calc.absorption_coeff, 
    opacity_calc.absorption_coeff, 
    opacity_calc.cloud_optical_depths, 
    opacity_calc.cloud_single_scattering, 
    opacity_calc.cloud_asym_param,
    radius_distance_scaling, 
    spectrum);


  postProcessSpectrum(spectrum, model_spectrum_bands);


  size_t nb_previous_param = nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_total_cloud_param;

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
bool EmissionModel::calcModelGPU(
  const std::vector<double>& parameter, 
  double* model_spectrum_gpu, 
  double* model_spectrum_bands)
{ 
  bool neglect = calcAtmosphereStructure(parameter);


  std::vector<double> cloud_parameters(
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param,
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_total_cloud_param);

  opacity_calc.calculateGPU(cloud_models, cloud_parameters);


  const double radius_distance_scaling = radiusDistanceScaling(parameter);

  radiative_transfer->calcSpectrumGPU(
    atmosphere,
    opacity_calc.absorption_coeff_gpu, 
    opacity_calc.scattering_coeff_dev, 
    opacity_calc.cloud_optical_depths_dev,
    opacity_calc.cloud_single_scattering_dev,
    opacity_calc.cloud_asym_param_dev,
    radius_distance_scaling,
    model_spectrum_gpu);


  postProcessSpectrumGPU(model_spectrum_gpu, model_spectrum_bands);


  size_t nb_previous_param = nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_total_cloud_param;

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



//integrate the high-res spectrum to observational bands
//and convolve if necessary 
void EmissionModel::postProcessSpectrum(
  std::vector<double>& model_spectrum, 
  std::vector<double>& model_spectrum_bands)
{
  model_spectrum_bands.assign(nb_observation_points, 0.0);
  
  std::vector<double>::iterator it = model_spectrum_bands.begin();

  for (size_t i=0; i<observations.size(); ++i)
  {
    const bool is_flux = true;

    std::vector<double> observation_bands = 
      observations[i].processModelSpectrum(model_spectrum, is_flux);

    //copy the band-integrated values for this observation into the global
    //vector of all band-integrated points, model_spectrum_bands
    std::copy(observation_bands.begin(), observation_bands.end(), it);
    it += observation_bands.size();
  }
}


//integrate the high-res spectrum to observational bands
//and convolve if necessary 
void EmissionModel::postProcessSpectrumGPU(
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



EmissionModel::~EmissionModel()
{
  delete radiative_transfer;
  delete temperature_profile;
  
  for (auto & i : cloud_models)
    delete i;

  for (auto & i : chemistry)
    delete i;
}



}

