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
#include "../../cloud_model/fixed_cloud_model.h"


namespace bear{


TransmissionModel::TransmissionModel (
  const TransmissionModelConfig model_config,
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
        model_config.cloud_model.size() > 0)
{
  nb_grid_points = model_config.nb_grid_points;
  opacity_species_symbol_ = model_config.opacity_species_symbol;
  opacity_species_folder_ = model_config.opacity_species_folder;

  std::cout << "Forward model selected: Transmission\n\n";

  //this forward model has three free general parameters
  nb_general_param = 3;

  if (model_config.fit_mean_molecular_weight || model_config.fit_scale_height)
    nb_general_param += 1;

  fit_mean_molecular_weight = model_config.fit_mean_molecular_weight;
  fit_scale_height = model_config.fit_scale_height;
  use_variable_gravity = model_config.use_variable_gravity;
  
  //select and set up the modules
  initModules(model_config);
}



TransmissionModel::TransmissionModel (
  GlobalConfig* config_, 
  SpectralGrid* spectral_grid_,
  const size_t nb_grid_points_,
  const std::vector<std::string>& opacity_species_symbol,
  const std::vector<std::string>& opacity_species_folder)
    : ForwardModel(config_, spectral_grid_, std::vector<Observation>() = {})
    , atmosphere(
        nb_grid_points_,
        std::vector<double>(2) = {10, 1e-5},
        config->use_gpu)
    , opacity_calc(
        config,
        spectral_grid,
        &atmosphere,
        opacity_species_symbol,
        opacity_species_folder,
        config->use_gpu,
        true)
{
  nb_grid_points = nb_grid_points_;

}



void TransmissionModel::extractParameters(
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

  module_parameters = std::vector<double>(
    parameters.begin() + nb_previous_param, 
    parameters.begin() + nb_previous_param + nb_total_modules_param);

  nb_previous_param += nb_total_modules_param;

  spectrum_modifier_parameters = std::vector<double>(
    parameters.begin() + nb_previous_param, 
    parameters.begin() + nb_previous_param + nb_spectrum_modifier_param);
}



bool TransmissionModel::calcAtmosphereStructure(const std::vector<double>& parameter)
{
  const double surface_gravity = std::pow(10,parameter[0]);
  const double bottom_radius = parameter[1];

  bool neglect_model = false;

  if (!fit_mean_molecular_weight && !fit_scale_height)
  {
    neglect_model = atmosphere.calcAtmosphereStructure(
      surface_gravity,
      bottom_radius,
      use_variable_gravity,
      temperature_profile.get(),
      temperature_parameters,
      chemistry,
      chemistry_parameters);

    return neglect_model;
  }

  //either mean molecular weight or scale height
  const double param = parameter[3];

  if (fit_mean_molecular_weight)
    neglect_model = atmosphere.calcAtmosphereStructure(
      surface_gravity,
      bottom_radius,
      use_variable_gravity,
      temperature_profile.get(),
      temperature_parameters,
      chemistry,
      chemistry_parameters,
      param);
  else
    neglect_model = atmosphere.calcAtmosphereStructure(
      surface_gravity,
      param,
      temperature_profile.get(),
      temperature_parameters,
      chemistry,
      chemistry_parameters);

  return neglect_model;
}



//Helper to get the parameter offset for module at given index
static size_t moduleParamOffset(
  const std::vector<std::unique_ptr<Module>>& modules, size_t idx)
{
  size_t offset = 0;
  for (size_t k = 0; k < idx; ++k)
    offset += modules[k]->nbParameters();
  return offset;
}


//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool TransmissionModel::calcModelCPU(
  const std::vector<double>& parameters,
  std::vector<double>& spectrum,
  std::vector<std::vector<double>>& spectrum_obs)
{
  extractParameters(parameters);

  bool neglect = calcAtmosphereStructure(parameters);

  const double bottom_radius = parameters[1];
  const double star_radius = parameters[2];

  // === Low-res path (skip if no low-res spectral points) ===
  if (spectral_grid->nbSpectralPoints() > 0)
  {
    opacity_calc.calculate(cloud_models, cloud_parameters);

    spectrum.assign(spectral_grid->nbSpectralPoints(), 0.0);

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

    if (modules_lowres_idx.empty() && modules_highres_idx.empty())
    {
      // No module classification (no high-res grid): apply all modules to low-res
      auto param_it = module_parameters.begin();

      for (auto & m : modules)
      {
        std::vector<double> p(param_it, param_it + m->nbParameters());
        m->modifySpectrum(p, &atmosphere, spectrum);
        param_it += m->nbParameters();
      }
    }
    else
    {
      // Apply only low-res modules (e.g. stellar contamination)
      for (size_t idx : modules_lowres_idx)
      {
        size_t offset = moduleParamOffset(modules, idx);
        std::vector<double> p(
          module_parameters.begin() + offset,
          module_parameters.begin() + offset + modules[idx]->nbParameters());
        modules[idx]->modifySpectrum(p, &atmosphere, spectrum);
      }
    }

    convertSpectrumToObservation(spectrum, false, spectrum_obs);
    applyObservationModifier(spectrum_modifier_parameters, spectrum_obs);
  }

  // === High-res path ===
  if (opacity_calc_highres)
  {
    const size_t nb_hr = spectral_grid_highres->nbSpectralPoints();

    opacity_calc_highres->calculate(cloud_models, std::vector<double>{});

    spectrum_highres_.assign(nb_hr, 0.0);
    calcTransmissionSpectrum(
      bottom_radius, star_radius,
      *opacity_calc_highres, nb_hr, spectrum_highres_);

    for (size_t idx : modules_highres_idx)
    {
      size_t offset = moduleParamOffset(modules, idx);
      std::vector<double> p(
        module_parameters.begin() + offset,
        module_parameters.begin() + offset + modules[idx]->nbParameters());
      modules[idx]->modifySpectrum(p, &atmosphere, spectrum_highres_);
    }
  }

  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool TransmissionModel::calcModelGPU(
  const std::vector<double>& parameter,
  float* spectrum,
  std::vector<float*>& spectrum_obs)
{
  extractParameters(parameter);

  bool neglect = calcAtmosphereStructure(parameter);

  const double bottom_radius = parameter[1];
  const double star_radius = parameter[2];

  // === Low-res path (skip if no low-res spectral points) ===
  if (spectral_grid->nbSpectralPoints() > 0)
  {
    opacity_calc.calculateGPU(cloud_models, cloud_parameters);

    if (cloud_models.size() > 0)
    {
      if (cloud_extinction_gpu == nullptr)
        allocateOnDevice(
          cloud_extinction_gpu,
          nb_grid_points*spectral_grid->nbSpectralPoints());

      initializeOnDevice(
        cloud_extinction_gpu,
        nb_grid_points*spectral_grid->nbSpectralPoints());

      cloud_models[0]->convertOpticalDepthGPU(
        opacity_calc.cloud_optical_depths_dev,
        atmosphere.altitude_dev,
        nb_grid_points,
        spectral_grid->nbSpectralPoints(),
        cloud_extinction_gpu);
    }

    calcTransitDepthGPU(
      spectrum,
      opacity_calc.absorption_coeff_gpu,
      opacity_calc.scattering_coeff_dev,
      cloud_extinction_gpu,
      atmosphere,
      spectral_grid->nbSpectralPoints(),
      bottom_radius,
      star_radius);

    if (modules_lowres_idx.empty() && modules_highres_idx.empty())
    {
      // No module classification: apply all modules to low-res
      auto param_it = module_parameters.begin();

      for (auto & m : modules)
      {
        std::vector<double> p(param_it, param_it + m->nbParameters());
        m->modifySpectrumGPU(p, &atmosphere, spectrum);
        param_it += m->nbParameters();
      }
    }
    else
    {
      for (size_t idx : modules_lowres_idx)
      {
        size_t offset = moduleParamOffset(modules, idx);
        std::vector<double> p(
          module_parameters.begin() + offset,
          module_parameters.begin() + offset + modules[idx]->nbParameters());
        modules[idx]->modifySpectrumGPU(p, &atmosphere, spectrum);
      }
    }

    convertSpectrumToObservationGPU(spectrum, false, spectrum_obs);
    applyObservationModifierGPU(spectrum_modifier_parameters, spectrum_obs);
  }

  // === High-res path ===
  if (opacity_calc_highres)
  {
    const size_t nb_hr = spectral_grid_highres->nbSpectralPoints();

    opacity_calc_highres->calculateGPU(cloud_models, std::vector<double>{});

    calcTransitDepthGPU(
      spectrum_highres_gpu_,
      opacity_calc_highres->absorption_coeff_gpu,
      opacity_calc_highres->scattering_coeff_dev,
      nullptr,
      atmosphere,
      nb_hr,
      bottom_radius,
      star_radius);

    // Diagnostic: check spectrum after transit depth, before broadening
    // static bool first_gpu_call = true;
    // if (first_gpu_call)
    // {
    //   std::vector<float> diag_spec(nb_hr);
    //   float* ptr = spectrum_highres_gpu_;
    //   moveToHost(ptr, diag_spec);
    //   float dmin = diag_spec[0], dmax = diag_spec[0];
    //   for (size_t i = 1; i < nb_hr; ++i)
    //   {
    //     if (diag_spec[i] < dmin) dmin = diag_spec[i];
    //     if (diag_spec[i] > dmax) dmax = diag_spec[i];
    //   }
    //   std::cout << "[GPU diag] after calcTransitDepthGPU: nb_hr=" << nb_hr
    //             << " min=" << dmin << " max=" << dmax
    //             << " Rp=" << bottom_radius << " Rs=" << star_radius << "\n";
    //   first_gpu_call = false;
    // }

    for (size_t idx : modules_highres_idx)
    {
      size_t offset = moduleParamOffset(modules, idx);
      std::vector<double> p(
        module_parameters.begin() + offset,
        module_parameters.begin() + offset + modules[idx]->nbParameters());
      modules[idx]->modifySpectrumGPU(p, &atmosphere, spectrum_highres_gpu_);
    }
  }

  return neglect;
}



void TransmissionModel::setCloudProperties(
  const std::vector<std::vector<double>>& cloud_optical_depth)
{
  if (config->use_gpu && cloud_extinction_gpu == nullptr)
    allocateOnDevice(
      cloud_extinction_gpu, 
      nb_grid_points*spectral_grid->nbSpectralPoints());
  else
    cloud_extinction.assign(
      spectral_grid->nbSpectralPoints(), 
      std::vector<double>(nb_grid_points, 0.0));
  
  bool use_cloud = false;

  for (auto & i : cloud_optical_depth)
  {
    double max = *std::max_element(i.begin(), i.end());
    if (max > 0)
    {
      use_cloud = true;
      break;
    }
  }
  
  if (!use_cloud)
    return;

  std::vector<std::vector<double>> single_scattering_albedo(
    nb_grid_points-1, 
    std::vector<double>(spectral_grid->nbSpectralPoints(), 0.0));

  std::vector<std::vector<double>> asymmetry_parameter(
    nb_grid_points-1, 
    std::vector<double>(spectral_grid->nbSpectralPoints(), 0.0));

  cloud_models.push_back(std::make_unique<FixedCloudModel>(
    cloud_optical_depth,
    single_scattering_albedo,
    asymmetry_parameter));

}



std::vector<double> TransmissionModel::calcSpectrum(
  const double surface_gravity,
  const double planet_radius,
  const double radius_ratio,
  const std::vector<double>& pressure,
  const std::vector<double>& temperature,
  const std::vector<std::string>& species_symbol,
  const std::vector<std::vector<double>>& mixing_ratios,
  const std::vector<std::vector<double>>& cloud_optical_depth,
  const double use_variable_gravity)
{
  atmosphere.setAtmosphericStructure(
    surface_gravity, 
    planet_radius,
    use_variable_gravity,
    pressure, 
    temperature, 
    species_symbol, 
    mixing_ratios);

  const double bottom_radius = planet_radius;
  const double star_radius = planet_radius/radius_ratio;
  
  setCloudProperties(cloud_optical_depth);

  std::vector<double> spectrum(spectral_grid->nbSpectralPoints(), 0.0);
  
  if (config->use_gpu)
  {
    opacity_calc.calculateGPU(cloud_models, std::vector<double> {});
    
    if (cloud_models.size() > 0)
    { 
      if (cloud_extinction_gpu == nullptr)
        allocateOnDevice(
          cloud_extinction_gpu, 
          nb_grid_points*spectral_grid->nbSpectralPoints());

      initializeOnDevice(
        cloud_extinction_gpu, 
        nb_grid_points*spectral_grid->nbSpectralPoints());

      cloud_models[0]->convertOpticalDepthGPU(
        opacity_calc.cloud_optical_depths_dev,
        atmosphere.altitude_dev,
        nb_grid_points,
        spectral_grid->nbSpectralPoints(),
        cloud_extinction_gpu);
    }

    float* model_spectrum_gpu = nullptr;

    allocateOnDevice(model_spectrum_gpu, spectral_grid->nbSpectralPoints());

    calcTransitDepthGPU(
      model_spectrum_gpu,
      opacity_calc.absorption_coeff_gpu,
      opacity_calc.scattering_coeff_dev,
      cloud_extinction_gpu,
      atmosphere,
      spectral_grid->nbSpectralPoints(),
      bottom_radius,
      star_radius);

    {
      std::vector<float> spectrum_float(spectral_grid->nbSpectralPoints());
      moveToHost(model_spectrum_gpu, spectrum_float);
      deleteFromDevice(model_spectrum_gpu);
      spectrum.assign(spectrum_float.begin(), spectrum_float.end());
    }
    deleteFromDevice(cloud_extinction_gpu);
  }
  else
  {
    opacity_calc.calculate(cloud_models, std::vector<double> {});

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
  }

  cloud_models.clear();

  return spectrum;
}



TransmissionModel::~TransmissionModel()
{
  if (cloud_extinction_gpu != nullptr)
    deleteFromDevice(cloud_extinction_gpu);

  if (spectrum_highres_gpu_ != nullptr)
    deleteFromDevice(spectrum_highres_gpu_);
}


}

