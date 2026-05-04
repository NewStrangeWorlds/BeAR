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

#include "phase_curve.h"

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../observations/observations.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"
#include "../atmosphere/atmosphere.h"
#include "../../radiative_transfer/select_radiative_transfer.h"
#include "../../cloud_model/fixed_cloud_model.h"


namespace bear{


PhaseCurveModel::PhaseCurveModel (
  const PhaseCurveConfig model_config,
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
  radiative_transfer_model_ = model_config.radiative_transfer_model;
  radiative_transfer_parameters_ = model_config.radiative_transfer_parameters;

  std::cout << "Forward model selected: Phase Curve\n\n";

  nb_general_param = 1;

  initModules(model_config);
}



void PhaseCurveModel::extractParameters(
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



bool PhaseCurveModel::calcAtmosphereStructure(const std::vector<double>& parameter)
{
  const double surface_gravity = std::pow(10,parameter[0]);

  bool neglect_model = false;

  neglect_model = atmosphere.calcAtmosphereStructure(
    surface_gravity,
    1.0,
    false,
    temperature_profile.get(),
    temperature_parameters,
    chemistry,
    chemistry_parameters);

  return neglect_model;
}



static size_t moduleParamOffset(
  const std::vector<std::unique_ptr<Module>>& modules, size_t idx)
{
  size_t offset = 0;
  for (size_t k = 0; k < idx; ++k)
    offset += modules[k]->nbParameters();
  return offset;
}


bool PhaseCurveModel::calcModelCPU(
  const std::vector<double>& parameters,
  std::vector<double>& spectrum,
  std::vector<std::vector<double>>& spectrum_obs)
{
  extractParameters(parameters);

  bool neglect = calcAtmosphereStructure(parameters);

  // === Low-res path ===
  if (spectral_grid->nbSpectralPoints() > 0)
  {
    opacity_calc.calculate(cloud_models, cloud_parameters);

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

    if (modules_lowres_idx.empty() && modules_highres_idx.empty())
    {
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
      for (size_t idx : modules_lowres_idx)
      {
        size_t offset = moduleParamOffset(modules, idx);
        std::vector<double> p(
          module_parameters.begin() + offset,
          module_parameters.begin() + offset + modules[idx]->nbParameters());
        modules[idx]->modifySpectrum(p, &atmosphere, spectrum);
      }
    }

    convertSpectrumToObservation(spectrum, true, spectrum_obs);
    applyObservationModifier(spectrum_modifier_parameters, spectrum_obs);
  }

  // === High-res path ===
  if (opacity_calc_highres)
  {
    const size_t nb_hr = spectral_grid_highres->nbSpectralPoints();

    opacity_calc_highres->calculate(cloud_models, std::vector<double>{});

    spectrum_highres_.assign(nb_hr, 0.0);

    radiative_transfer_highres->calcSpectrum(
      atmosphere,
      opacity_calc_highres->absorption_coeff,
      opacity_calc_highres->absorption_coeff,
      opacity_calc_highres->cloud_optical_depths,
      opacity_calc_highres->cloud_single_scattering,
      opacity_calc_highres->cloud_asym_param,
      1.0,
      spectrum_highres_);

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



bool PhaseCurveModel::calcModelGPU(
  const std::vector<double>& parameters,
  float* spectrum,
  std::vector<float*>& spectrum_obs)
{
  extractParameters(parameters);

  bool neglect = calcAtmosphereStructure(parameters);

  neglect = false;

  // === Low-res path ===
  if (spectral_grid->nbSpectralPoints() > 0)
  {
    opacity_calc.calculateGPU(cloud_models, cloud_parameters);

    radiative_transfer->calcSpectrumGPU(
      atmosphere,
      opacity_calc.absorption_coeff_gpu,
      opacity_calc.scattering_coeff_dev,
      opacity_calc.cloud_optical_depths_dev,
      opacity_calc.cloud_single_scattering_dev,
      opacity_calc.cloud_asym_param_dev,
      1.0,
      spectrum);

    if (modules_lowres_idx.empty() && modules_highres_idx.empty())
    {
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

    convertSpectrumToObservationGPU(spectrum, true, spectrum_obs);
    applyObservationModifierGPU(spectrum_modifier_parameters, spectrum_obs);
  }

  // === High-res path ===
  if (opacity_calc_highres)
  {
    opacity_calc_highres->calculateGPU(cloud_models, std::vector<double>{});

    radiative_transfer_highres->calcSpectrumGPU(
      atmosphere,
      opacity_calc_highres->absorption_coeff_gpu,
      opacity_calc_highres->scattering_coeff_dev,
      opacity_calc_highres->cloud_optical_depths_dev,
      opacity_calc_highres->cloud_single_scattering_dev,
      opacity_calc_highres->cloud_asym_param_dev,
      1.0,
      spectrum_highres_gpu_);

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



void PhaseCurveModel::setCloudProperties(
  const std::vector<std::vector<double>>& cloud_optical_depth)
{
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



std::vector<double> PhaseCurveModel::calcSpectrum(
      const double surface_gravity,
      const std::vector<double>& pressure,
      const std::vector<double>& temperature,
      const std::vector<std::string>& species_symbol,
      const std::vector<std::vector<double>>& mixing_ratios,
      const std::vector<std::vector<double>>& cloud_optical_depth)
{
  atmosphere.setAtmosphericStructure(
    surface_gravity,
    1.0,
    false,
    pressure,
    temperature,
    species_symbol,
    mixing_ratios);

  setCloudProperties(cloud_optical_depth);

  std::vector<double> spectrum(spectral_grid->nbSpectralPoints(), 0.0);

  if (config->use_gpu)
  {
    opacity_calc.calculateGPU(cloud_models, std::vector<double> {});

    float* model_spectrum_gpu = nullptr;
    allocateOnDevice(model_spectrum_gpu, spectral_grid->nbSpectralPoints());

    radiative_transfer->calcSpectrumGPU(
      atmosphere,
      opacity_calc.absorption_coeff_gpu,
      opacity_calc.scattering_coeff_dev,
      opacity_calc.cloud_optical_depths_dev,
      opacity_calc.cloud_single_scattering_dev,
      opacity_calc.cloud_asym_param_dev,
      1.0,
      model_spectrum_gpu);

    {
      std::vector<float> spectrum_float(spectral_grid->nbSpectralPoints());
      moveToHost(model_spectrum_gpu, spectrum_float);
      deleteFromDevice(model_spectrum_gpu);
      spectrum.assign(spectrum_float.begin(), spectrum_float.end());
    }
  }
  else
  {
    opacity_calc.calculate(cloud_models, std::vector<double> {});

    radiative_transfer->calcSpectrum(
      atmosphere,
      opacity_calc.absorption_coeff,
      opacity_calc.absorption_coeff,
      opacity_calc.cloud_optical_depths,
      opacity_calc.cloud_single_scattering,
      opacity_calc.cloud_asym_param,
      1.0,
      spectrum);
  }

  cloud_models.clear();

  return spectrum;
}



PhaseCurveModel::~PhaseCurveModel()
{
}



}

