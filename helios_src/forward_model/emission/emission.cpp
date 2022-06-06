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


#include "emission.h"


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



EmissionModel::EmissionModel (Retrieval* retrieval_ptr, const EmissionModelConfig model_config) 
 : transport_coeff(retrieval_ptr->config, &retrieval_ptr->spectral_grid, model_config.opacity_species_symbol, model_config.opacity_species_folder),
   atmosphere(model_config.nb_grid_points, model_config.atmos_boundaries, retrieval_ptr->config->use_gpu)
{
  retrieval = retrieval_ptr;
  nb_grid_points = model_config.nb_grid_points;
  
  std::cout << "Forward model selected: Emission\n\n"; 

  //this forward model has three free general parameters
  nb_general_param = 3;

  //select and set up the modules
  initModules(model_config);


  //allocate memory for the absorption coefficients on the GPU if necessary
  if (retrieval->config->use_gpu)
    initDeviceMemory();


  setPriors();
}



//calculates the radius distance scaling factor from the MultiNest parameters
double EmissionModel::radiusDistanceScaling(const std::vector<double>& parameter)
{
  const double scaling_factor = parameter[1];
  const double distance = parameter[2] * constants::parsec;
  const double prior_radius = constants::radius_jupiter;  //we assume a fixed prior radius of 1 Rj

  double scaling = prior_radius/distance;
  scaling = scaling*scaling * scaling_factor;

  return scaling;
}



//determines the basic atmospheric structure (temperature profile, chemistry...) from the free parameters supplied by MultiNest
bool EmissionModel::calcAtmosphereStructure(const std::vector<double>& parameter)
{
  const double surface_gravity = std::pow(10,parameter[0]);
  const double scaling_factor = parameter[1];

  //derived radius in Jupiter radii assuming that the radius prior is 1 Rj
  const double derived_radius = std::sqrt(scaling_factor);  

  //derived mass in Jupiter masses
  const double derived_mass = surface_gravity * std::pow(derived_radius*constants::radius_jupiter, 2) / constants::gravitation_const / constants::mass_jupiter;


  bool neglect_model = false;

  //if derived mass is larger than 80 Jupiter masses, we tell MultiNest to neglect this parameter combination 
  if (derived_mass > 80) neglect_model = true;


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
bool EmissionModel::calcModel(const std::vector<double>& parameter, std::vector<double>& spectrum, std::vector<double>& model_spectrum_bands)
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

  
  //calculate gas absorption coefficients
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
  const double radius_distance_scaling = radiusDistanceScaling(parameter);
 
  radiative_transfer->calcSpectrum(atmosphere,
                                   absorption_coeff, scattering_coeff, 
                                   cloud_optical_depths, cloud_single_scattering, cloud_asym_param,
                                   radius_distance_scaling, 
                                   spectrum);


  postProcessSpectrum(spectrum, model_spectrum_bands);

  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool EmissionModel::calcModelGPU(const std::vector<double>& parameter, double* model_spectrum_gpu, double* model_spectrum_bands)
{ 
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


  const double radius_distance_scaling = radiusDistanceScaling(parameter);

  radiative_transfer->calcSpectrumGPU(atmosphere,
                                      absorption_coeff_gpu, 
                                      nullptr,
                                      cloud_optical_depths_dev,
                                      cloud_single_scattering_dev,
                                      cloud_asym_param_dev,
                                      radius_distance_scaling,
                                      model_spectrum_gpu);


  postProcessSpectrumGPU(model_spectrum_gpu, model_spectrum_bands);

  return neglect;
}



//integrate the high-res spectrum to observational bands
//and convolve if necessary 
void EmissionModel::postProcessSpectrum(std::vector<double>& model_spectrum, std::vector<double>& model_spectrum_bands)
{
  model_spectrum_bands.assign(retrieval->nb_observation_points, 0.0);
  
  std::vector<double>::iterator it = model_spectrum_bands.begin();


  for (size_t i=0; i<retrieval->nb_observations; ++i)
  {
    std::vector<double> observation_bands;

    std::vector<double> spectrum = model_spectrum;
    
    //apply filter function if necessary
    if (retrieval->observations[i].filter_response.size() != 0)
      spectrum = retrieval->observations[i].applyFilterResponseFunction(spectrum);


    if (retrieval->observations[i].instrument_profile_fwhm.size() == 0)
      retrieval->observations[i].spectral_bands.bandIntegrateSpectrum(spectrum, observation_bands);
    else
    { //apply instrument line profile
      std::vector<double> model_spectrum_convolved = retrieval->observations[i].spectral_bands.convolveSpectrum(spectrum);
      retrieval->observations[i].spectral_bands.bandIntegrateSpectrum(model_spectrum_convolved, observation_bands);
    }
  
    //copy the band-integrated values for this observation into the global
    //vector of all band-integrated points, model_spectrum_bands
    std::copy(observation_bands.begin(), observation_bands.end(), it);
    it += observation_bands.size();
  }
}


//integrate the high-res spectrum to observational bands
//and convolve if necessary 
void EmissionModel::postProcessSpectrumGPU(double* model_spectrum_gpu, double* model_spectrum_bands)
{

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


  //integrate the high-res spectrum to observational bands on the GPU
  bandIntegrationGPU(retrieval->band_spectrum_id,
                     retrieval->band_indices_gpu,
                     retrieval->band_sizes_gpu,
                     retrieval->band_start_index_gpu,
                     retrieval->nb_total_bands,
                     retrieval->spectral_grid.wavenumber_list_gpu,
                     retrieval->spectral_grid.wavelength_list_gpu,
                     model_spectrum_bands);
}




EmissionModel::~EmissionModel()
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

