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
    : ForwardModel(config_, spectral_grid_, observations_)
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
{
  nb_grid_points = model_config.nb_grid_points;
  
  std::cout << "Forward model selected: Emission\n\n"; 

  //this forward model has three free general parameters
  nb_general_param = 3;

  initModules(model_config);

  setPriors(priors_);
}


double EmissionModel::radiusDistanceScaling(const std::vector<double>& parameter)
{
  const double scaling_factor = parameter[1];
  const double distance = parameter[2];

   //we assume a fixed prior radius of 1 Rj
  const double prior_radius = constants::radius_jupiter; 
  
  double scaling = prior_radius/distance;
  scaling = scaling*scaling * scaling_factor;

  return scaling;
}


void EmissionModel::extractParameters(
  const std::vector<double>& parameters)
{
  model_parameters = std::vector<double>(
    parameters.begin(), 
    parameters.begin() + nb_general_param);

  size_t nb_previous_param = nb_general_param;
  
  chemistry_parameters = std::vector<double>(
    parameters.begin() + nb_previous_param, 
    parameters.begin() + nb_previous_param + nb_total_chemistry_param);

  nb_previous_param += nb_total_chemistry_param;
  
  temperature_parameters = std::vector<double>(
    parameters.begin() + nb_previous_param, 
    parameters.begin() + nb_previous_param + nb_temperature_param);

  nb_previous_param += nb_temperature_param;

  cloud_parameters = std::vector<double>(
      parameters.begin() + nb_previous_param,
      parameters.begin() + nb_previous_param + nb_total_cloud_param);

  nb_previous_param += nb_total_cloud_param;

  spectrum_modifier_parameters = std::vector<double>(
    parameters.begin() + nb_previous_param, 
    parameters.begin() + nb_previous_param + nb_spectrum_modifier_param);
}




bool EmissionModel::calcAtmosphereStructure(const std::vector<double>& parameter)
{
  const double surface_gravity = std::pow(10,model_parameters[0]);
  const double scaling_factor = model_parameters[1];

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

  neglect_model = atmosphere.calcAtmosphereStructure(
    surface_gravity, 
    1.0,
    false,
    temperature_profile, 
    temperature_parameters, 
    chemistry, 
    chemistry_parameters);


  return neglect_model;
}



//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool EmissionModel::calcModel(
  const std::vector<double>& parameters, 
  std::vector<double>& spectrum, 
  std::vector<std::vector<double>>& spectrum_obs)
{
  extractParameters(parameters);

  bool neglect = calcAtmosphereStructure(parameters);

  opacity_calc.calculate(cloud_models, cloud_parameters);


  spectrum.assign(spectral_grid->nbSpectralPoints(), 0.0);
  
  const double radius_distance_scaling = radiusDistanceScaling(model_parameters);

  radiative_transfer->calcSpectrum(
    atmosphere,
    opacity_calc.absorption_coeff, 
    opacity_calc.absorption_coeff, 
    opacity_calc.cloud_optical_depths, 
    opacity_calc.cloud_single_scattering, 
    opacity_calc.cloud_asym_param,
    radius_distance_scaling, 
    spectrum);


  convertSpectrumToObservation(spectrum, true, spectrum_obs);
  
  applyObservationModifier(spectrum_modifier_parameters, spectrum_obs);

  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool EmissionModel::calcModelGPU(
  const std::vector<double>& parameters, 
  double* spectrum, 
  std::vector<double*>& spectrum_obs)
{ 
  extractParameters(parameters);

  bool neglect = calcAtmosphereStructure(parameters);

  opacity_calc.calculateGPU(cloud_models, cloud_parameters);


  const double radius_distance_scaling = radiusDistanceScaling(model_parameters);

  radiative_transfer->calcSpectrumGPU(
    atmosphere,
    opacity_calc.absorption_coeff_gpu, 
    opacity_calc.scattering_coeff_dev, 
    opacity_calc.cloud_optical_depths_dev,
    opacity_calc.cloud_single_scattering_dev,
    opacity_calc.cloud_asym_param_dev,
    radius_distance_scaling,
    spectrum);

  convertSpectrumToObservationGPU(spectrum, true, spectrum_obs);

  applyObservationModifierGPU(spectrum_modifier_parameters, spectrum_obs);

  return neglect;
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

