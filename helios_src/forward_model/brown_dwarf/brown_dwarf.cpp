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


#include "brown_dwarf.h"


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



BrownDwarfModel::BrownDwarfModel (Retrieval* retrieval_ptr, const BrownDwarfConfig model_config) 
 : transport_coeff(retrieval_ptr->config, &retrieval_ptr->spectral_grid, model_config.opacity_species_symbol, model_config.opacity_species_folder),
   atmosphere(model_config.nb_grid_points, model_config.atmos_boundaries)
{
  retrieval = retrieval_ptr;
  nb_grid_points = model_config.nb_grid_points;
  use_cloud_layer = model_config.use_cloud_layer;
  
  std::cout << "Forward model selected: Emission\n\n"; 

  //allocate memory for the absorption coefficients on the GPU if necessary
  if (retrieval->config->use_gpu)
    allocateOnDevice(absorption_coeff_gpu, nb_grid_points*retrieval->spectral_grid.nbSpectralPoints());

  //select and set up the modules
  initChemistry(model_config);
  initTemperature(model_config);
  initRadiativeTransfer(model_config);


  setPriors();
}



//calculates the upper and lower grid point of the cloud based on the top and bottom pressure
void BrownDwarfModel::calcCloudPosition(const double top_pressure, const double bottom_pressure, unsigned int& top_index, unsigned int& bottom_index)
{
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    if ((atmosphere.pressure[i] > top_pressure && atmosphere.pressure[i+1] < top_pressure) || atmosphere.pressure[i] == top_pressure )
      top_index = i;

    if ((atmosphere.pressure[i] > bottom_pressure && atmosphere.pressure[i+1] < bottom_pressure) || atmosphere.pressure[i] == bottom_pressure )
      bottom_index = i;
  }


  if (bottom_pressure > atmosphere.pressure[0])
    bottom_index = 0;


  //clouds needs to occupy at least an entire atmospheric layer
  if (top_index == bottom_index)
    bottom_index -= 2;
}




//calculates the radius distance scaling factor from the MultiNest parameters
double BrownDwarfModel::radiusDistanceScaling(const double distance, const double radius, const double scaling_f)
{
  double scaling = radius/distance;
  scaling = scaling*scaling * scaling_f;

  return scaling;
}



//determines the basic atmospheric structure (temperature profile, chemistry...) from the free parameters supplied by MultiNest
//also returns the radius-distance scaling parameter
bool BrownDwarfModel::calcAtmosphereStructure(const std::vector<double>& parameter)
{
  const double surface_gravity = std::pow(10,parameter[0]);
  const double scaling_factor = parameter[1];

  const double distance = parameter[2] * constants::parsec;
  const double prior_radius = constants::radius_jupiter;
  

  radius_distance_scaling = radiusDistanceScaling(distance, prior_radius, scaling_factor);


  //derived radius in Jupiter radii
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


  //optional cloud layer
  if (use_cloud_layer)
  {
    std::vector<double> cloud_parameters(parameter.begin() + nb_general_param + nb_total_chemistry_param + temperature_profile->nbParameters(),
                                         parameter.begin() + nb_general_param + nb_total_chemistry_param + temperature_profile->nbParameters() + nb_cloud_param);
    calcGreyCloudLayer(cloud_parameters);
  }

  return neglect_model;
}



//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool BrownDwarfModel::calcModel(const std::vector<double>& parameter, std::vector<double>& spectrum)
{
  bool neglect = calcAtmosphereStructure(parameter);


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
                                   cloud_optical_depths, atmosphere.temperature, atmosphere.altitude, 
                                   spectrum);


  for (size_t i=0; i<retrieval->spectral_grid.nbSpectralPoints(); ++i)
    spectrum[i] *= radius_distance_scaling;
  

  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool BrownDwarfModel::calcModelGPU(const std::vector<double>& parameter, double* model_spectrum_gpu, double* model_spectrum_bands)
{ 
  bool neglect = calcAtmosphereStructure(parameter);


  initCrossSectionsHost(retrieval->spectral_grid.nbSpectralPoints()*nb_grid_points, absorption_coeff_gpu);


  for (size_t i=0; i<nb_grid_points; ++i)
    transport_coeff.calcTransportCoefficientsGPU(atmosphere.temperature[i], atmosphere.pressure[i], atmosphere.number_densities[i],
                                                 nb_grid_points, i,
                                                 absorption_coeff_gpu, nullptr);

 
  radiative_transfer->calcSpectrumGPU(model_spectrum_gpu,
                                      absorption_coeff_gpu, 
                                      nullptr,
                                      retrieval->spectral_grid.wavenumber_list_gpu,
                                      cloud_optical_depths,
                                      atmosphere.temperature, atmosphere.altitude,
                                      radius_distance_scaling);


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

  //convolveHSTSpectrumGPU(model_spectrum_bands, retrieval->nb_total_bands);

  return neglect;
}



//calculates the vertical distribution of the grey layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
void BrownDwarfModel::calcGreyCloudLayer(const std::vector<double>& cloud_parameters)
{
  double cloud_top_pressure = cloud_parameters[0];
  double cloud_bottom_pressure = cloud_top_pressure * cloud_parameters[1];
  double cloud_optical_depth = cloud_parameters[2];

  
  unsigned int cloud_top_index = 0;
  unsigned int cloud_bottom_index = 0;

  calcCloudPosition(cloud_top_pressure, cloud_bottom_pressure, cloud_top_index, cloud_bottom_index);

  double cloud_optical_depth_layer = cloud_optical_depth/static_cast<double>(cloud_top_index - cloud_bottom_index); 

  cloud_optical_depths.assign(nb_grid_points-1, 0);

  for (size_t i=cloud_bottom_index; i<cloud_top_index; ++i)
    cloud_optical_depths[i] = cloud_optical_depth_layer;
}




BrownDwarfModel::~BrownDwarfModel()
{
  if (retrieval->config->use_gpu)
    deleteFromDevice(absorption_coeff_gpu);

  delete radiative_transfer;
  
  for (auto & i : chemistry)
    delete i;

  delete temperature_profile;
}



}

