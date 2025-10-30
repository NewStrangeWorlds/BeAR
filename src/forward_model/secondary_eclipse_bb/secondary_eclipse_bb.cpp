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

#include "secondary_eclipse_bb.h"

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../observations/observations.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"

#include "../stellar_spectrum/stellar_spectrum.h"

namespace bear{


OccultationBlackBodyModel::OccultationBlackBodyModel (
  const OccultationBlackBodyConfig model_config,
  GlobalConfig* config_,
  SpectralGrid* spectral_grid_,
  std::vector<Observation>& observations_) 
    : ForwardModel(config_, spectral_grid_, observations_)
{
  std::cout << "Forward model selected: Secondary Eclipse Black Body\n\n";

  //this forward model has two free general parameters
  nb_general_param = 2;
  
  initModules(model_config);
}


void OccultationBlackBodyModel::extractParameters(
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

  spectrum_modifier_parameters = std::vector<double>(
    parameters.begin() + nb_previous_param, 
    parameters.begin() + nb_previous_param + nb_spectrum_modifier_param);
}



//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool OccultationBlackBodyModel::calcModelCPU(
  const std::vector<double>& parameters, 
  std::vector<double>& spectrum, 
  std::vector<std::vector<double>>& spectrum_obs)
{
  extractParameters(parameters);
  
  const double planet_temperature = model_parameters[0];

  spectrum.assign(spectral_grid->nbSpectralPoints(), 0.0);

  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<spectral_grid->nbSpectralPoints(); ++i)
    spectrum[i] = aux::planckFunctionWavenumber(planet_temperature, spectral_grid->wavenumber_list[i]) * constants::pi * 1e-3;
  
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


  spectrum_obs.assign(nb_observation_points, std::vector<double>{});

  const double radius_ratio = model_parameters[1];

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


  //convert the high-res spectrum to an occulation depth as well
  for (size_t i=0; i<spectral_grid->nbSpectralPoints(); ++i)
    spectrum[i] = spectrum[i]/stellar_spectrum[i] * radius_ratio*radius_ratio * 1e6;

  return false;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool OccultationBlackBodyModel::calcModelGPU(
  const std::vector<double>& parameters, 
  double* spectrum, 
  std::vector<double*>& spectrum_obs)
{
  extractParameters(parameters);

  const double radius_ratio = model_parameters[1];
  const double planet_temperature = model_parameters[0];

  calcPlanetSpectrumGPU(planet_temperature, spectrum);


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


  //convert the original high-res spectrum also to a secondary eclipse
  calcOccultationGPU(
    spectrum, 
    spectrum, 
    stellar_spectrum, 
    spectral_grid->nbSpectralPoints(),
    radius_ratio, 
    albedo_contribution_gpu);

  deleteFromDevice(albedo_contribution_gpu);
  deleteFromDevice(stellar_spectrum);

  return false;
}



OccultationBlackBodyModel::~OccultationBlackBodyModel()
{
  delete stellar_model;

}


}

