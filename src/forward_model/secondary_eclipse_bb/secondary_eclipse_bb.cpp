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
#include "../../retrieval/priors.h"
#include "../../observations/observations.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"

#include "../stellar_spectrum/stellar_spectrum.h"

namespace bear{


SecondaryEclipseBlackBodyModel::SecondaryEclipseBlackBodyModel (
  const SecondaryEclipseBlackBodyConfig model_config,
  Priors* priors_,
  GlobalConfig* config_,
  SpectralGrid* spectral_grid_,
  std::vector<Observation>& observations_) 
    : config(config_)
    , spectral_grid(spectral_grid_)
    , observations(observations_)
{
  std::cout << "Forward model selected: Secondary Eclipse Black Body\n\n";

  //this forward model has three free general parameters
  nb_general_param = 2;

  for (auto & i : observations)
  {
    nb_observation_points += i.nbPoints();
    nb_spectrum_modifier_param += i.nb_modifier_param;
  }

  initModules(model_config);

  setPriors(priors_);
}


//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool SecondaryEclipseBlackBodyModel::calcModel(
  const std::vector<double>& parameter, 
  std::vector<double>& spectrum, 
  std::vector<double>& model_spectrum_bands)
{
  const double planet_temperature = parameter[0];

  spectrum.assign(spectral_grid->nbSpectralPoints(), 0.0);

  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<spectral_grid->nbSpectralPoints(); ++i)
    spectrum[i] = aux::planckFunctionWavenumber(planet_temperature, spectral_grid->wavenumber_list[i]) * constants::pi * 1e-3;

  //post-process the planet's high-res emission spectrum and bin it to the observational bands
  std::vector<double> planet_spectrum_bands(nb_observation_points, 0);
  postProcessSpectrum(
    spectrum, 
    planet_spectrum_bands);


  std::vector<double> stellar_parameters(
      parameter.begin() + nb_general_param,
      parameter.begin() + nb_general_param + nb_stellar_param);

  std::vector<double> stellar_spectrum = stellar_model->calcFlux(stellar_parameters);

  //post-process the star's high-res emission spectrum and bin it to the observational bands
  std::vector<double> stellar_spectrum_bands(nb_observation_points, 0);
  postProcessSpectrum(
    stellar_spectrum, 
    stellar_spectrum_bands);


  model_spectrum_bands.assign(nb_observation_points, 0);
  
  double radius_distance_ratio = 0;
  std::vector<double> albedo_contribution_bands(nb_observation_points, 0.0);
  std::vector<double> albedo_contribution(spectral_grid->nbSpectralPoints(), 0.0);

  const double radius_ratio = parameter[1];

  //calculate the secondary-eclipse depth with an optional geometric albedo contribution
  for (size_t i=0; i<model_spectrum_bands.size(); ++i)
    model_spectrum_bands[i] = 
      planet_spectrum_bands[i]/stellar_spectrum_bands[i] 
      * radius_ratio*radius_ratio * 1e6 + albedo_contribution_bands[i]*1e6;


  size_t nb_previous_param = nb_general_param + nb_stellar_param;

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


  //convert the high-res spectrum to an eclipse depth as well
  for (size_t i=0; i<spectral_grid->nbSpectralPoints(); ++i)
    spectrum[i] = spectrum[i]/stellar_spectrum[i] * radius_ratio*radius_ratio * 1e6
                  + radius_distance_ratio * albedo_contribution[i] * 1e6;

  return false;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool SecondaryEclipseBlackBodyModel::calcModelGPU(
  const std::vector<double>& parameter, 
  double* model_spectrum_gpu, 
  double* model_spectrum_bands)
{
  const double radius_ratio = parameter[1];
  
  const double planet_temperature = parameter[0];

  calcPlanetSpectrumGPU(planet_temperature, model_spectrum_gpu);

  //post-process the planet's high-res emission spectrum and bin it to the observational bands
  double* planet_spectrum_bands = nullptr;
  allocateOnDevice(planet_spectrum_bands, nb_observation_points);

  postProcessSpectrumGPU(
    model_spectrum_gpu, 
    planet_spectrum_bands);


  std::vector<double> stellar_parameters(
      parameter.begin() + nb_general_param,
      parameter.begin() + nb_general_param + nb_stellar_param);
  
  double* stellar_spectrum = nullptr;
  allocateOnDevice(stellar_spectrum, spectral_grid->nbSpectralPoints());
  stellar_model->calcFluxGPU(stellar_parameters, stellar_spectrum);
  
  //post-process the planet's high-res emission spectrum and bin it to the observational bands
  double* stellar_spectrum_bands = nullptr;
  allocateOnDevice(stellar_spectrum_bands, nb_observation_points);

  postProcessSpectrumGPU(
    stellar_spectrum, 
    stellar_spectrum_bands);


  //std::vector<double> albedo_contribution(retrieval->nb_total_bands, geometric_albedo*radius_distance_ratio*radius_distance_ratio);
  std::vector<double> albedo_contribution(nb_observation_points, 0.0);

  double* albedo_contribution_gpu = nullptr;
  double* albedo_contribution_bands_gpu = nullptr;
  //moveToDevice(albedo_contribution_gpu, albedo_contribution);

  calcSecondaryEclipseGPU(
    model_spectrum_bands, 
    planet_spectrum_bands, 
    stellar_spectrum_bands, 
    nb_observation_points,
    radius_ratio, 
    albedo_contribution_bands_gpu);
  
  deleteFromDevice(albedo_contribution_bands_gpu);
  deleteFromDevice(planet_spectrum_bands);
  deleteFromDevice(stellar_spectrum_bands);


  size_t nb_previous_param = nb_general_param + nb_stellar_param;

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


  //convert the original high-res spectrum also to a secondary eclipse
  calcSecondaryEclipseGPU(
    model_spectrum_gpu, 
    model_spectrum_gpu, 
    stellar_spectrum, 
    spectral_grid->nbSpectralPoints(),
    radius_ratio, 
    albedo_contribution_gpu);

  deleteFromDevice(albedo_contribution_gpu);
  deleteFromDevice(stellar_spectrum);

  return false;
}




//integrate the high-res spectrum to observational bands
//and convolve if necessary 
void SecondaryEclipseBlackBodyModel::postProcessSpectrum(
  std::vector<double>& model_spectrum, 
  std::vector<double>& model_spectrum_bands)
{
  model_spectrum_bands.assign(nb_observation_points, 0.0);

  std::vector<double>::iterator it = model_spectrum_bands.begin();

  for (size_t i=0; i<observations.size(); ++i)
  { 
    const bool is_flux = true;

    std::vector<double> observation_bands = observations[i].processModelSpectrum(
      model_spectrum, 
      is_flux);

    //copy the band-integrated values for this observation into the global
    //vector of all band-integrated points, model_spectrum_bands
    std::copy(observation_bands.begin(), observation_bands.end(), it);
    it += observation_bands.size();
  }
}


//integrate the high-res spectrum to observational bands
//and convolve if necessary 
void SecondaryEclipseBlackBodyModel::postProcessSpectrumGPU(
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


SecondaryEclipseBlackBodyModel::~SecondaryEclipseBlackBodyModel()
{

}


}

