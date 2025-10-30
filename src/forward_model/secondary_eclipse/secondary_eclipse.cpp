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

#include "secondary_eclipse.h"

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../observations/observations.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"
#include "../atmosphere/atmosphere.h"
#include "../../radiative_transfer/select_radiative_transfer.h"

#include "../stellar_spectrum/stellar_spectrum.h"
#include "../stellar_spectrum/star_file_spectrum.h"
#include "../../cloud_model/fixed_cloud_model.h"


namespace bear{


OccultationModel::OccultationModel (
  const OccultationConfig model_config,
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

  std::cout << "Forward model selected: Secondary Eclipse\n\n";

  //this forward model has two free general parameters
  nb_general_param = 2;

  initModules(model_config);
}



OccultationModel::OccultationModel (
  GlobalConfig* config_, 
  SpectralGrid* spectral_grid_,
  const size_t nb_grid_points_,
  const std::vector<double>& stellar_spectrum_wavelengths,
  const std::vector<double>& stellar_spectrum_flux,
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

  stellar_model =  new StarSpectrumFile(
    stellar_spectrum_wavelengths,
    stellar_spectrum_flux,
    spectral_grid);

  radiative_transfer = selectRadiativeTransfer(
    std::string("scm"), 
    std::vector<std::string> {}, 
    nb_grid_points, 
    config, 
    spectral_grid);
}


void OccultationModel::extractParameters(
  const std::vector<double>& parameters)
{
  model_parameters = std::vector<double>(
    parameters.begin(), 
    parameters.begin() + nb_general_param);

  size_t nb_previous_param = nb_general_param;

  stellar_parameters = std::vector<double>(
      parameters.begin() + nb_previous_param,
      parameters.begin() + nb_previous_param + nb_stellar_param);
  
  nb_previous_param += nb_stellar_param;
  
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



//determines the basic atmospheric structure (temperature profile, chemistry...) from the free parameters supplied by MultiNest
bool OccultationModel::calcAtmosphereStructure(const std::vector<double>& parameter)
{
  const double surface_gravity = std::pow(10,parameter[0]);
 
  bool neglect_model = false;
  
  //determine atmosphere structure
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
bool OccultationModel::calcModelCPU(
  const std::vector<double>& parameters, 
  std::vector<double>& spectrum, 
  std::vector<std::vector<double>>& spectrum_obs)
{
  extractParameters(parameters);

  bool neglect = calcAtmosphereStructure(parameters);


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


  //post-process the planet's high-res emission spectrum and bin it to the observational bands
  std::vector<std::vector<double>> planet_spectrum_obs(
    observations.size(), 
    std::vector<double>{});
  
  convertSpectrumToObservation(
    spectrum, 
    true,
    planet_spectrum_obs);
  

  std::vector<double> stellar_spectrum = stellar_model->calcFlux(stellar_parameters);

  //post-process the star's high-res emission spectrum and bin it to the observational bands
  std::vector<std::vector<double>> stellar_spectrum_obs(
    observations.size(), 
    std::vector<double>{});
  
  convertSpectrumToObservation(
    stellar_spectrum, 
    true,
    stellar_spectrum_obs);


  const double radius_ratio = model_parameters[1];

  //model_spectrum_bands.assign(nb_observation_points, 0);

  //convert the spectra into an occulation depth
  for (size_t i=0; i<observations.size(); ++i)
  {
    spectrum_obs[i].assign(observations[i].nbPoints(), 0.0);
    
    for (size_t j=0; j<spectrum_obs[i].size(); ++j)
    {
      spectrum_obs[i][j] = planet_spectrum_obs[i][j]/stellar_spectrum_obs[i][j] 
        * radius_ratio*radius_ratio * 1e6;
    }
  }
  

  applyObservationModifier(spectrum_modifier_parameters, spectrum_obs);


  //convert the high-res planet spectrum to an occulation depth as well
  for (size_t i=0; i<spectral_grid->nbSpectralPoints(); ++i)
    spectrum[i] = spectrum[i]/stellar_spectrum[i] * radius_ratio*radius_ratio * 1e6;

  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool OccultationModel::calcModelGPU(
  const std::vector<double>& parameters, 
  double* spectrum, 
  std::vector<double*>& spectrum_obs)
{
  extractParameters(parameters);

  const double radius_ratio = model_parameters[1];

  bool neglect = calcAtmosphereStructure(parameters);


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


  std::vector<double*> planet_spectrum_obs(observations.size(), nullptr);
  std::vector<double*> stellar_spectrum_obs(observations.size(), nullptr);

  for (size_t i=0; i<observations.size(); ++i)
  {
    allocateOnDevice(planet_spectrum_obs[i], observations[i].nbPoints());
    allocateOnDevice(stellar_spectrum_obs[i], observations[i].nbPoints());
  }

  convertSpectrumToObservationGPU(
    spectrum, 
    true,
    planet_spectrum_obs);


  double* stellar_spectrum = nullptr;
  allocateOnDevice(stellar_spectrum, spectral_grid->nbSpectralPoints());
  
  stellar_model->calcFluxGPU(stellar_parameters, stellar_spectrum);
  

  convertSpectrumToObservationGPU(
    stellar_spectrum, 
    true,
    stellar_spectrum_obs);


  double* albedo_contribution_gpu = nullptr;
  double* albedo_contribution_bands_gpu = nullptr;
  //moveToDevice(albedo_contribution_gpu, albedo_contribution);
  
  for (size_t i=0; i<observations.size(); ++i)
  {
    calcOccultationGPU(
    spectrum_obs[i], 
    planet_spectrum_obs[i], 
    stellar_spectrum_obs[i], 
    observations[i].nbPoints(),
    radius_ratio, 
    albedo_contribution_bands_gpu);
  }
  
  
  deleteFromDevice(albedo_contribution_bands_gpu);

  for (size_t i=0; i<observations.size(); ++i)
  {
    deleteFromDevice(planet_spectrum_obs[i]);
    deleteFromDevice(stellar_spectrum_obs[i]);
  }

 
  applyObservationModifierGPU(spectrum_modifier_parameters, spectrum_obs);


  //convert the original high-res planet spectrum also to a secondary eclipse
  calcOccultationGPU(
    spectrum, 
    spectrum, 
    stellar_spectrum, 
    spectral_grid->nbSpectralPoints(),
    radius_ratio, 
    albedo_contribution_gpu);

  deleteFromDevice(albedo_contribution_gpu);
  deleteFromDevice(stellar_spectrum);

  return neglect;
}



void OccultationModel::setCloudProperties(
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

  FixedCloudModel* model = new FixedCloudModel(
    cloud_optical_depth,
    single_scattering_albedo,
    asymmetry_parameter);

  cloud_models.push_back(model);
}



std::vector<double> OccultationModel::calcSpectrum(
      const double surface_gravity,
      const double radius_ratio,
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

    double* model_spectrum_gpu = nullptr;

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

    double* stellar_spectrum = nullptr;
    allocateOnDevice(stellar_spectrum, spectral_grid->nbSpectralPoints());
    stellar_model->calcFluxGPU(std::vector<double> {}, stellar_spectrum);

    double* albedo_contribution_gpu = nullptr;

    calcOccultationGPU(
      model_spectrum_gpu, 
      model_spectrum_gpu, 
      stellar_spectrum, 
      spectral_grid->nbSpectralPoints(),
      radius_ratio, 
      albedo_contribution_gpu);

    deleteFromDevice(albedo_contribution_gpu);
    deleteFromDevice(stellar_spectrum);

    moveToHost(model_spectrum_gpu, spectrum);
    deleteFromDevice(model_spectrum_gpu);
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
    
    std::vector<double> stellar_spectrum = stellar_model->calcFlux(std::vector<double> {});

    for (size_t i=0; i<spectral_grid->nbSpectralPoints(); ++i)
      spectrum[i] = spectrum[i] / stellar_spectrum[i] * radius_ratio*radius_ratio * 1e6;
  }

  if (cloud_models.size() > 0)
  {
    delete cloud_models[0];
    cloud_models.clear();
  }

  return spectrum;
}



OccultationModel::~OccultationModel()
{
  delete radiative_transfer;
  delete temperature_profile;
  delete stellar_model;
  
  for (auto & i : cloud_models)
    delete i;
  
  for (auto & i : chemistry)
    delete i;
}



}

