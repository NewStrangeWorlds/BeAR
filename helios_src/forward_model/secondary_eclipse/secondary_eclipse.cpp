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


#include "../../chemistry/chem_species.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/quadrature.h"
#include "../../additional/exceptions.h"
#include "../../retrieval/retrieval.h"




namespace helios{



SecondaryEclipseModel::SecondaryEclipseModel (Retrieval* retrieval_ptr, const SecondaryEclipseConfig model_config) 
 : transport_coeff(retrieval_ptr->config, &retrieval_ptr->spectral_grid, model_config.opacity_species_symbol, model_config.opacity_species_folder)
{
  retrieval = retrieval_ptr;
  nb_grid_points = model_config.nb_grid_points;
  use_cloud_layer = model_config.use_cloud_layer;


  //allocate memory for the absorption coefficients on the GPU if necessary
  if (retrieval->config->use_gpu)
    allocateOnDevice(absorption_coeff_gpu, nb_grid_points*retrieval->spectral_grid.nbSpectralPoints());

  
  createPressureGrid(model_config.atmos_boundaries);

  
  //initialise temperatures, altitudes, and number_densities to 0
  temperature.assign(nb_grid_points, 0.0);
  z_grid.assign(nb_grid_points, 0.0);
  number_densities.assign(nb_grid_points, std::vector<double>(constants::species_data.size(), 0.0));

  
  //select and set up the modules
  initStellarSpectrum(model_config);
  initChemistry(model_config);
  initTemperature(model_config);
  initRadiativeTransfer(model_config);

  
  setPriors();
}


//Creates a pressure grid between the two boundary points
//nb_grid_point pressures vales are equidistantly spread in the log(p) space
void SecondaryEclipseModel::createPressureGrid(const double atmos_boundaries [2])
{
  pressure.assign(nb_grid_points, 0.0);

  const double min_pressure = atmos_boundaries[1];
  const double max_pressure = atmos_boundaries[0];

  
  pressure.front() = std::log10(max_pressure);

  const double log_step = (std::log10(max_pressure) - std::log10(min_pressure)) / (nb_grid_points - 1.0);


  for (size_t i=1; i<nb_grid_points-1; ++i)
    pressure[i] = pressure[i-1] - log_step;


  for (size_t i=0; i<nb_grid_points-1; ++i)
    pressure[i] = pow(10.0, pressure[i]);

  pressure.back() = min_pressure;
}



//calculates the upper and lower grid point of the cloud based on the top and bottom pressure
void SecondaryEclipseModel::calcCloudPosition(const double top_pressure, const double bottom_pressure, unsigned int& top_index, unsigned int& bottom_index)
{
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    if ((pressure[i] > top_pressure && pressure[i+1] < top_pressure) || pressure[i] == top_pressure )
      top_index = i;

    if ((pressure[i] > bottom_pressure && pressure[i+1] < bottom_pressure) || pressure[i] == bottom_pressure )
      bottom_index = i;
  }


  if (bottom_pressure > pressure[0])
    bottom_index = 0;


  //clouds needs to occupy at least an entire atmospheric layer
  if (top_index == bottom_index)
    bottom_index -= 2;
}



//determine the vertical grid via hydrostatic equilibrium
void SecondaryEclipseModel::calcAltitude(const double g, const std::vector<double>& mean_molecular_weights)
{
  z_grid.assign(nb_grid_points, 0.0);
  std::vector<double> mass_density(nb_grid_points, 0.0);

  for (size_t i=0; i<nb_grid_points; ++i)
    mass_density[i] = mean_molecular_weights[i] * pressure[i]*1e6 / (constants::gas_constant  * temperature[i]);

  for (size_t i=1; i<nb_grid_points; i++)
  {
    double delta_z = (1.0/(mass_density[i]*g) + 1.0/(mass_density[i-1]*g)) * 0.5 * (pressure[i-1]*1e6 - pressure[i]*1e6);

    z_grid[i] = z_grid[i-1] + delta_z;
  }
  
}


//determines the basic atmospheric structure (temperature profile, chemistry...) from the free parameters supplied by MultiNest
//also returns the radius-distance scaling parameter
bool SecondaryEclipseModel::calcAtmosphereStructure(const std::vector<double>& parameter)
{
  const double surface_gravity = std::pow(10,parameter[0]);
 
  bool neglect_model = false;


  //temperature profile
  std::vector<double> temp_parameters(parameter.begin() + nb_general_param + nb_total_chemistry_param, 
                                      parameter.begin() + nb_general_param + nb_total_chemistry_param + temperature_profile->nbParameters());

  bool neglect_temperature = temperature_profile->calcProfile(temp_parameters, pressure, temperature);


  if (neglect_temperature) neglect_model = true;

  //chemical composition
  std::vector<double> mean_molecular_weights(nb_grid_points, 0.0);
  number_densities.assign(nb_grid_points, std::vector<double>(constants::species_data.size(), 0.0));
  size_t nb_chem_param = 0;

  for (auto & i : chemistry)
  {
    std::vector<double> chem_parameters(parameter.begin() + nb_general_param + nb_chem_param, 
                                        parameter.begin() + nb_general_param + nb_chem_param + i->nbParameters());
    nb_chem_param += i->nbParameters();
    
    bool neglect = i->calcChemicalComposition(chem_parameters, temperature, pressure, number_densities, mean_molecular_weights);

    if (neglect) neglect_model = true;
  } 
 
 
  //optional cloud layer
  if (use_cloud_layer)
  {
    std::vector<double> cloud_parameters(parameter.begin() + nb_general_param + nb_total_chemistry_param + temperature_profile->nbParameters(),
                                         parameter.begin() + nb_general_param + nb_total_chemistry_param + temperature_profile->nbParameters() + nb_cloud_param);
    calcGreyCloudLayer(cloud_parameters);
  }
  

  calcAltitude(surface_gravity, mean_molecular_weights);


  return neglect_model;
}



//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool SecondaryEclipseModel::calcModel(const std::vector<double>& parameter, std::vector<double>& spectrum)
{
  bool neglect = calcAtmosphereStructure(parameter);


  absorption_coeff.assign(retrieval->spectral_grid.nbSpectralPoints(), std::vector<double>(nb_grid_points, 0.0));
  scattering_coeff.assign(retrieval->spectral_grid.nbSpectralPoints(), std::vector<double>(nb_grid_points, 0.0));

  
  for (size_t i=0; i<nb_grid_points; ++i)
  {
    std::vector<double> absorption_coeff_level(retrieval->spectral_grid.nbSpectralPoints(), 0.0);
    std::vector<double> scattering_coeff_level(retrieval->spectral_grid.nbSpectralPoints(), 0.0);


    transport_coeff.calcTransportCoefficients(temperature[i], pressure[i], number_densities[i], absorption_coeff_level, scattering_coeff_level);

    
    for (size_t j=0; j<retrieval->spectral_grid.nbSpectralPoints(); ++j)
      absorption_coeff[j][i] = absorption_coeff_level[j];
  }


  spectrum.assign(retrieval->spectral_grid.nbSpectralPoints(), 0.0);
 
  radiative_transfer->calcSpectrum(absorption_coeff, scattering_coeff, cloud_optical_depths, temperature, z_grid, spectrum);


  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool SecondaryEclipseModel::calcModelGPU(const std::vector<double>& parameter, double* model_spectrum_gpu, double* model_spectrum_bands)
{
  const double radius_ratio = parameter[1];
  const double radius_distance_ratio = parameter[1]/parameter[2];
  const double geometric_albedo = parameter[3];


  bool neglect = calcAtmosphereStructure(parameter);


  initCrossSectionsHost(retrieval->spectral_grid.nbSpectralPoints()*nb_grid_points, absorption_coeff_gpu);


  for (size_t i=0; i<nb_grid_points; ++i)
    transport_coeff.calcTransportCoefficientsGPU(temperature[i], pressure[i], number_densities[i],
                                                 nb_grid_points, i,
                                                 absorption_coeff_gpu, nullptr);


  radiative_transfer->calcSpectrumGPU(model_spectrum_gpu,
                                      absorption_coeff_gpu, 
                                      nullptr,
                                      retrieval->spectral_grid.wavenumber_list_gpu,
                                      cloud_optical_depths,
                                      temperature, z_grid,
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
                             retrieval->convolved_spectra[i]);
    
    /*std::cout << i << "\t" << retrieval->convolved_spectra[i] << "\t" << model_spectrum_gpu << "\n";
    std::vector<double> test_spectrum(retrieval->spectral_grid.nbSpectralPoints(), 0.0);
    moveToHost(retrieval->convolved_spectra[i], test_spectrum);

    for (size_t j=0; j<retrieval->spectral_grid.nbSpectralPoints(); ++j)
      std::cout << j << "\t" << retrieval->spectral_grid.wavelength_list[j] << "\t" << test_spectrum[j] << "\n";
    std::cout << "\n\n";*/
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
  //albedo_contribution[retrieval->nb_total_bands-1] = 0.0;
  //albedo_contribution[retrieval->nb_total_bands-2] = 0.0;

  albedo_contribution[0] = geometric_albedo*radius_distance_ratio*radius_distance_ratio;
  //albedo_contribution[1] = geometric_albedo*radius_distance_ratio*radius_distance_ratio;


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



//calculates the vertical distribution of the grey layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
void SecondaryEclipseModel::calcGreyCloudLayer(const std::vector<double>& cloud_parameters)
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




SecondaryEclipseModel::~SecondaryEclipseModel()
{
  if (retrieval->config->use_gpu)
    deleteFromDevice(absorption_coeff_gpu);


  delete radiative_transfer;
  
  for (auto & i : chemistry)
    delete i;

  delete temperature_profile;
}



}

