/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
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


#include "secondary_eclipse.h"


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>


#include "../../CUDA_kernels/data_management_kernels.h"
#include "../../CUDA_kernels/cross_section_kernels.h"
#include "../../CUDA_kernels/filter_response_kernels.h"
#include "../../CUDA_kernels/band_integration_kernels.h"
#include "../../CUDA_kernels/convolution_kernels.h"


#include "../../chemistry/chem_species.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/quadrature.h"
#include "../../additional/exceptions.h"
#include "../../retrieval/retrieval.h"

#include "../atmosphere/atmosphere.h"



namespace helios{



SecondaryEclipseModel::SecondaryEclipseModel (Retrieval* retrieval_ptr, const SecondaryEclipseConfig model_config) 
 : transport_coeff(retrieval_ptr->config, &retrieval_ptr->spectral_grid, model_config.opacity_species_symbol, model_config.opacity_species_folder),
   atmosphere(model_config.nb_grid_points, model_config.atmos_boundaries, retrieval_ptr->config->use_gpu)
{
  retrieval = retrieval_ptr;
  nb_grid_points = model_config.nb_grid_points;
  use_cloud_layer = model_config.use_cloud_layer;

  std::cout << "Forward model selected: Secondary Eclipse\n\n";

  //this forward model has three free general parameters
  nb_general_param = 3;


  //select and set up the modules
  initStellarSpectrum(model_config);
  initModules(model_config);


  //allocate memory for the absorption coefficients on the GPU if necessary
  if (retrieval->config->use_gpu)
    initDeviceMemory();


  setPriors();
}



//determines the basic atmospheric structure (temperature profile, chemistry...) from the free parameters supplied by MultiNest
//also returns the radius-distance scaling parameter
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
bool SecondaryEclipseModel::calcModel(const std::vector<double>& parameter, std::vector<double>& spectrum)
{
  bool neglect = calcAtmosphereStructure(parameter);


  cloud_optical_depths.assign(retrieval->spectral_grid.nbSpectralPoints(), std::vector<double>(nb_grid_points-1, 0.0));
  cloud_single_scattering.assign(retrieval->spectral_grid.nbSpectralPoints(), std::vector<double>(nb_grid_points-1, 0.0));
  cloud_asym_param.assign(retrieval->spectral_grid.nbSpectralPoints(), std::vector<double>(nb_grid_points-1, 0.0));

  //calculate cloud model if needed
  if (cloud_model != nullptr)
  {
    std::vector<double> cloud_parameters(parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param,
                                      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_cloud_param);
   
    cloud_model->opticalProperties(cloud_parameters, atmosphere, &retrieval->spectral_grid, cloud_optical_depths, cloud_single_scattering, cloud_asym_param);
  }


  absorption_coeff.assign(retrieval->spectral_grid.nbSpectralPoints(), std::vector<double>(nb_grid_points, 0.0));
  scattering_coeff.assign(retrieval->spectral_grid.nbSpectralPoints(), std::vector<double>(nb_grid_points, 0.0));

  
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    std::vector<double> absorption_coeff_level(retrieval->spectral_grid.nbSpectralPoints(), 0.0);
    std::vector<double> scattering_coeff_level(retrieval->spectral_grid.nbSpectralPoints(), 0.0);


    transport_coeff.calcTransportCoefficients(atmosphere.temperature[i], atmosphere.pressure[i], atmosphere.number_densities[i], 
                                              absorption_coeff_level, scattering_coeff_level);

    
    for (size_t j=0; j<retrieval->spectral_grid.nbSpectralPoints(); ++j)
      absorption_coeff[j][i] = absorption_coeff_level[j];
  }


  spectrum.assign(retrieval->spectral_grid.nbSpectralPoints(), 0.0);
 
  radiative_transfer->calcSpectrum(absorption_coeff, scattering_coeff, 
                                   cloud_optical_depths, cloud_single_scattering, cloud_asym_param,
                                   atmosphere.temperature, atmosphere.altitude, 
                                   spectrum);


  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool SecondaryEclipseModel::calcModelGPU(const std::vector<double>& parameter, double* model_spectrum_gpu, double* model_spectrum_bands)
{
  const double radius_ratio = parameter[1];

  bool neglect = calcAtmosphereStructure(parameter);


  //calculate cloud model if needed
  if (cloud_model != nullptr)
  { 
    std::vector<double> cloud_parameters(parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param,
                                      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_cloud_param);
   
    cloud_model->opticalPropertiesGPU(cloud_parameters, atmosphere, &retrieval->spectral_grid, cloud_optical_depths_dev, cloud_single_scattering_dev, cloud_asym_param_dev);
  }


  initCrossSectionsHost(retrieval->spectral_grid.nbSpectralPoints()*nb_grid_points, absorption_coeff_gpu);


  for (size_t i=0; i<nb_grid_points; ++i)
    transport_coeff.calcTransportCoefficientsGPU(atmosphere.temperature[i], atmosphere.pressure[i], atmosphere.number_densities[i],
                                                 nb_grid_points, i,
                                                 absorption_coeff_gpu, nullptr);


  radiative_transfer->calcSpectrumGPU(atmosphere,
                                      model_spectrum_gpu,
                                      absorption_coeff_gpu, 
                                      nullptr,
                                      cloud_optical_depths_dev,
                                      cloud_single_scattering_dev,
                                      cloud_asym_param_dev,
                                      1.0);


  for (size_t i=0; i<retrieval->observations.size(); ++i)
  {
    if (retrieval->observations[i].filter_response.size() != 0) 
      applyFilterResponseGPU(retrieval->spectral_grid.wavenumber_list_gpu,
                             model_spectrum_gpu, 
                             retrieval->observations[i].filter_response_gpu, 
                             retrieval->observations[i].filter_response_weight_gpu, 
                             retrieval->observations[i].filter_response_normalisation,
                             retrieval->spectral_grid.nbSpectralPoints(),
                             retrieval->filter_response_spectra[i]);


    size_t nb_points_observation = retrieval->observations[i].spectral_bands.wavenumbers.size();

    if (retrieval->observations[i].instrument_profile_fwhm.size() != 0) 
      convolveSpectrumGPU(retrieval->filter_response_spectra[i], 
                          retrieval->observation_wavelengths[i], 
                          retrieval->observation_profile_sigma[i], 
                          retrieval->observation_spectral_indices[i],
                          retrieval->convolution_start_index[i], 
                          retrieval->convolution_end_index[i], 
                          nb_points_observation, 
                          retrieval->convolved_spectra[i]);
  }


  double* planet_spectrum_bands = nullptr;
  allocateOnDevice(planet_spectrum_bands, retrieval->nb_total_bands);

 
  //integrate the high-res spectrum to observational bands on the GPU
  bandIntegrationGPU(retrieval->band_spectrum_id,
                     retrieval->band_indices_gpu,
                     retrieval->band_sizes_gpu,
                     retrieval->band_start_index_gpu,
                     retrieval->nb_total_bands,
                     retrieval->spectral_grid.wavenumber_list_gpu,
                     retrieval->spectral_grid.wavelength_list_gpu,
                     planet_spectrum_bands);


  //std::vector<double> albedo_contribution(retrieval->nb_total_bands, geometric_albedo*radius_distance_ratio*radius_distance_ratio);
  std::vector<double> albedo_contribution(retrieval->nb_total_bands, 0.0);


  double* albedo_contribution_gpu = nullptr;

  moveToDevice(albedo_contribution_gpu, albedo_contribution);


  calcSecondaryEclipseGPU(model_spectrum_bands, planet_spectrum_bands, stellar_spectrum_bands_gpu, retrieval->nb_total_bands,
                          radius_ratio, albedo_contribution_gpu);
  
  deleteFromDevice(albedo_contribution_gpu);
  deleteFromDevice(planet_spectrum_bands);


  return neglect;
}




std::vector<double> SecondaryEclipseModel::calcSecondaryEclipse(std::vector<double>& planet_spectrum_bands, const double radius_ratio,
                                                                const double geometric_albedo, const double radius_distance_ratio)
{
  std::vector<double> secondary_eclipse(planet_spectrum_bands.size(), 0.0);

  for (size_t i=0; i<secondary_eclipse.size(); ++i)
    secondary_eclipse[i] = planet_spectrum_bands[i]/stellar_spectrum_bands[i] * radius_ratio*radius_ratio
                           + geometric_albedo * radius_distance_ratio*radius_distance_ratio;
  
  return secondary_eclipse;
}




SecondaryEclipseModel::~SecondaryEclipseModel()
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


  delete radiative_transfer;
  delete temperature_profile;
  delete cloud_model;
  
  for (auto & i : chemistry)
    delete i;
}



}

