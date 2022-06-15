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


#include "transmission.h"


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>


#include "../../chemistry/chem_species.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/quadrature.h"
#include "../../additional/exceptions.h"
#include "../../retrieval/retrieval.h"
#include "../atmosphere/atmosphere.h"
#include "../../CUDA_kernels/data_management_kernels.h"
#include "../../CUDA_kernels/cross_section_kernels.h"


namespace helios{



TransmissionModel::TransmissionModel (
  Retrieval* retrieval_ptr, 
  const TransmissionModelConfig model_config) 
    : transport_coeff(
        retrieval_ptr->config, 
        &retrieval_ptr->spectral_grid, 
        model_config.opacity_species_symbol, 
        model_config.opacity_species_folder)
    , atmosphere(
        model_config.nb_grid_points, 
        model_config.atmos_boundaries, 
        retrieval_ptr->config->use_gpu)
{
  retrieval = retrieval_ptr;
  nb_grid_points = model_config.nb_grid_points;
  
  std::cout << "Forward model selected: Emission\n\n"; 

  //this forward model has three free general parameters
  nb_general_param = 3;

  //select and set up the modules
  initModules(model_config);

  if (retrieval->config->use_gpu)
    initDeviceMemory();

  setPriors();
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
  

  //determine atmosphere structure
  neglect_model = atmosphere.calcAtmosphereStructure(
    surface_gravity, 
    temperature_profile, 
    temp_parameters, 
    chemistry, 
    chem_parameters);


  return neglect_model;
}



//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool TransmissionModel::calcModel(
  const std::vector<double>& parameter, 
  std::vector<double>& spectrum, 
  std::vector<double>& model_spectrum_bands)
{
  bool neglect = calcAtmosphereStructure(parameter);

  const size_t nb_spectral_points = retrieval->spectral_grid.nbSpectralPoints();

  cloud_optical_depths.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));
  cloud_single_scattering.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));
  cloud_asym_param.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));

  //calculate cloud model if needed
  if (cloud_model != nullptr)
  {
    std::vector<double> cloud_parameters(
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param,
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_cloud_param);
   
    cloud_model->opticalProperties(
      cloud_parameters, 
      atmosphere, 
      &retrieval->spectral_grid, 
      cloud_optical_depths, 
      cloud_single_scattering, 
      cloud_asym_param);
  }


  //calculate gas transport coefficients
  absorption_coeff.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));
  scattering_coeff.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    std::vector<double> absorption_coeff_level(nb_spectral_points, 0.0);
    std::vector<double> scattering_coeff_level(nb_spectral_points, 0.0);


    transport_coeff.calcTransportCoefficients(
      atmosphere.temperature[i],
      atmosphere.pressure[i], 
      atmosphere.number_densities[i], 
      absorption_coeff_level,
      scattering_coeff_level);

    
    for (size_t j=0; j<nb_spectral_points; ++j)
      absorption_coeff[j][i] = absorption_coeff_level[j];
  }


  spectrum.assign(nb_spectral_points, 0.0);

  const double bottom_radius = parameter[1] * constants::radius_jupiter;
  const double star_radius = parameter[2] * constants::radius_sun;

  calcTransmissionSpectrum(bottom_radius, star_radius, spectrum);


  postProcessSpectrum(spectrum, model_spectrum_bands);

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
  

  //calculate cloud model if needed
  if (cloud_model != nullptr)
  { 
    std::vector<double> cloud_parameters(
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param,
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_cloud_param);
   
    cloud_model->opticalPropertiesGPU(
      cloud_parameters, 
      atmosphere, 
      &retrieval->spectral_grid, 
      cloud_optical_depths_dev, 
      cloud_single_scattering_dev, 
      cloud_asym_param_dev);
  }


  initCrossSectionsHost(
    retrieval->spectral_grid.nbSpectralPoints()*nb_grid_points, 
    absorption_coeff_gpu);

  initCrossSectionsHost(
    retrieval->spectral_grid.nbSpectralPoints()*nb_grid_points, 
    scattering_coeff_dev);

  for (size_t i=0; i<nb_grid_points; ++i)
    transport_coeff.calcTransportCoefficientsGPU(
      atmosphere.temperature[i], 
      atmosphere.pressure[i], 
      atmosphere.number_densities[i],
      nb_grid_points, 
      i,
      absorption_coeff_gpu, 
      scattering_coeff_dev);


  const double bottom_radius = parameter[1] * constants::radius_jupiter;
  const double star_radius = parameter[2] * constants::radius_sun;

  calcTransitDepthGPU(
    model_spectrum_gpu, 
    absorption_coeff_gpu, 
    scattering_coeff_dev, 
    atmosphere, 
    retrieval->spectral_grid.nbSpectralPoints(), 
    bottom_radius, 
    star_radius);
  
  postProcessSpectrumGPU(model_spectrum_gpu, model_spectrum_bands);

  return neglect;
}



void TransmissionModel::postProcessSpectrum(
  std::vector<double>& model_spectrum, 
  std::vector<double>& model_spectrum_bands)
{
  model_spectrum_bands.assign(retrieval->nb_observation_points, 0.0);
  
  std::vector<double>::iterator it = model_spectrum_bands.begin();

  for (size_t i=0; i<retrieval->nb_observations; ++i)
  {
    const bool is_flux = false;
    std::vector<double> observation_bands = retrieval->observations[i].processModelSpectrum(model_spectrum, is_flux);

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
  for (size_t i=0; i<retrieval->observations.size(); ++i)
  {
    const bool is_flux = false;
    retrieval->observations[i].processModelSpectrumGPU(
      model_spectrum_gpu, 
      model_spectrum_bands, 
      start_index, 
      is_flux);

    start_index += retrieval->observations[i].spectral_bands.nbBands();
  }
}




TransmissionModel::~TransmissionModel()
{
  if (retrieval->config->use_gpu)
  {
    deleteFromDevice(absorption_coeff_gpu);
    deleteFromDevice(scattering_coeff_dev);

    if (cloud_model != nullptr)
    {
      deleteFromDevice(cloud_optical_depths_dev);
      deleteFromDevice(cloud_single_scattering_dev);
      deleteFromDevice(cloud_asym_param_dev);
    }
  }


  delete temperature_profile;
  delete cloud_model;
  
  for (auto & i : chemistry)
    delete i;
}



}

