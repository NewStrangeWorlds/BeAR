/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2019 Daniel Kitzmann
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


#include "../CUDA_kernels/data_management_kernels.h"
#include "../CUDA_kernels/cross_section_kernels.h"


#include "../transport_coeff/chem_species.h"
#include "../retrieval/retrieval.h"
#include "../additional/aux_functions.h"
#include "../additional/physical_const.h"
#include "../additional/quadrature.h"
#include "piecewise_poly.h"




namespace helios{


BrownDwarfModel::BrownDwarfModel (Retrieval* retrieval_ptr, const size_t& nb_points, const double atmos_boundaries [2])
 : temperature_pol(retrieval_ptr->config->nb_temperature_elements, retrieval_ptr->config->temperature_poly_degree, atmos_boundaries)
{

  retrieval = retrieval_ptr;
  nb_grid_points = nb_points;


  transport_coeff.init(retrieval_ptr->config, &retrieval_ptr->spectral_grid);

  
  //allocate memory for the absorption coefficients on the GPU if necessary
  if (retrieval->config->useGPU())
    allocateOnDevice(absorption_coeff_gpu, nb_grid_points*retrieval->spectral_grid.nbSpectralPoints());


  createPressureGrid(atmos_boundaries);

  
  //initialise temperatures, altitudes, and number_densities to 0
  temperature.assign(nb_grid_points, 0.0);
  z_grid.assign(nb_grid_points, 0.0);
  number_densities.assign(nb_grid_points, std::vector<double>(chemical_symbols.size(), 0.0));

  
  //initialise radiative transfer
  ShortCharacteristics* scm = new ShortCharacteristics(retrieval_ptr);
  radiative_transfer = scm;

  //DiscreteOrdinates* disort = new DiscreteOrdinates(retrieval_ptr, 4, nb_points); 
  //radiative_transfer = disort;
 
}



void BrownDwarfModel::createPressureGrid(const double atmos_boundaries [2])
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
void BrownDwarfModel::calcCloudPosition(const double top_pressure, const double bottom_pressure, unsigned int& top_index, unsigned int& bottom_index)
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
void BrownDwarfModel::calcAltitude(const double g, const std::vector<double>& mean_molecular_weights)
{

  z_grid.assign(nb_grid_points, 0.0);
  std::vector<double> mass_density(nb_grid_points, 0.0);

  for (size_t i=0; i<nb_grid_points; ++i)
    mass_density[i] = mean_molecular_weights[i] * pressure[i]*1e6 / (helios::CONST_R  * temperature[i]);

  for (size_t i=1; i<nb_grid_points; i++)
  {

    double delta_z = (1.0/(mass_density[i]*g) + 1.0/(mass_density[i-1]*g)) * 0.5 * (pressure[i-1]*1e6 - pressure[i]*1e6);

    z_grid[i] = z_grid[i-1] + delta_z;

  }
  
  

}




std::vector<double> BrownDwarfModel::temperatureProfile(const std::vector<double>& parameter, std::vector<double>& pressure_profile)
{

  calcAtmosphereStructure(parameter);

  pressure_profile = pressure;

  return temperature;

}



//calculates the radius distance scaling factor from the MultiNest parameters
double BrownDwarfModel::radiusDistanceScaling(const std::vector<double>& parameter)
{

  const double distance = parameter[2] * CONST_PARSEC;
  const double radius =  CONST_JUPITER_RADIUS;

  const double scaling_f = parameter[1];

  double scaling = radius/distance;
  scaling = scaling*scaling * scaling_f;

  return scaling;
}



//determines the basic atmospheric structure (temperature profile, chemistry...) from the free parameters supplied by MultiNest
//also returns the radius-distance scaling parameter
bool BrownDwarfModel::calcAtmosphereStructure(const std::vector<double>& parameter)
{
  
  //double metallicity_factor = parameter[0];
  //double co_ratio = parameter[1];
  const double surface_gravity = std::pow(10,parameter[0]);
  const double scaling_factor = parameter[1];
  
  radius_distance_scaling = radiusDistanceScaling(parameter);

  

  //derived radius in Jupiter radii
  const double derived_radius = std::sqrt(scaling_factor);  

  //derived mass in Jupiter masses
  const double derived_mass = surface_gravity * std::pow(derived_radius*CONST_JUPITER_RADIUS, 2) / CONST_G / CONST_JUPITER_MASS;


  bool neglect = false;

  //if derived mass is larger than 80 Jupiter masses, we tell MultiNest to neglect this parameter combination 
  if (derived_mass > 80) neglect = true;


  std::vector<double> temp_parameters(parameter.begin()+10, 
                                      parameter.begin()+10+retrieval->config->nb_temperature_elements*retrieval->config->temperature_poly_degree + 1);

  calcTemperature(temp_parameters);

  //neglect models with too low temperatures
  for (auto & i : temperature)
    if (i < 50) {i = 50; neglect = true;}


  std::vector<double> chem_parameters(parameter.begin()+3, parameter.begin()+3+9);
  std::vector<double> mean_molecular_weights(nb_grid_points, 0.0);
  
  calcFreeChemistry(chem_parameters, mean_molecular_weights);


  
  //std::vector<double> cloud_parameters(parameter.begin()+17, parameter.begin()+17+3);

  //calcGreyCloudLayer(cloud_parameters);



  calcAltitude(surface_gravity, mean_molecular_weights);


  return neglect;
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


    transport_coeff.calcTransportCoefficients(temperature[i], pressure[i], number_densities[i], absorption_coeff_level, scattering_coeff_level);


    for (size_t j=0; j<retrieval->spectral_grid.nbSpectralPoints(); ++j)
      absorption_coeff[j][i] = absorption_coeff_level[j];
  }


  spectrum.assign(retrieval->spectral_grid.nbSpectralPoints(), 0.0);
 
  radiative_transfer->calcSpectrum(absorption_coeff, scattering_coeff, temperature, z_grid, spectrum);

  

  for (size_t i=0; i<retrieval->spectral_grid.nbSpectralPoints(); ++i)
    spectrum[i] *= radius_distance_scaling;


  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool BrownDwarfModel::calcModelGPU(const std::vector<double>& parameter, double* model_spectrum_gpu)
{
 
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
                                      radius_distance_scaling);


  return neglect;
}



//calculate the temperature vie piecewise polynomials
void BrownDwarfModel::calcTemperature(const std::vector<double>& temp_parameter)
{


  size_t nb_dof = retrieval->config->nb_temperature_elements*retrieval->config->temperature_poly_degree + 1;
  
  std::vector<double> temperature_dof(nb_dof, 0.0);

  temperature_dof[0] = temp_parameter[0];

  for (size_t i=1; i<nb_dof; ++i)
    temperature_dof[i] = temperature_dof[i-1] * temp_parameter[i];

  temperature_pol.setDOFvalues(temperature_dof);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    double p_log = std::log10(pressure[i]);

    temperature[i] = temperature_pol.getValue(p_log);
  }


}



void BrownDwarfModel::calcFreeChemistry(const std::vector<double>& chem_parameters, std::vector<double>& mean_molecular_weights)
{

  const double mr_h2o = chem_parameters[0];
  const double mr_ch4 = chem_parameters[1];
  const double mr_nh3 = chem_parameters[2];
  const double mr_co2 = chem_parameters[3];
  const double mr_co = chem_parameters[4];
  const double mr_h2s = chem_parameters[5];
  const double mr_k = chem_parameters[6];


  const double solar_na_k = 16.2181;
  const double solar_h2 = 0.5;
  const double solar_he = 0.085114;
  const double solar_h2_he = solar_h2 + solar_he;
  
  const double mr_na = mr_k * solar_na_k;


  double mr_rest = 1.0 -  mr_h2o - mr_co - mr_co2 - mr_ch4 - mr_nh3 - mr_h2s - mr_na - mr_k;
  
  double mr_h2 = mr_rest * solar_h2 / solar_h2_he;
  double mr_he = mr_rest * solar_he / solar_h2_he;


  for (size_t i=0; i<nb_grid_points; ++i)
  {
    number_densities[i][_TOTAL] = pressure[i] * 1.e6 / helios::CONST_K / temperature[i];

    number_densities[i][_H2O] = number_densities[i][_TOTAL] * mr_h2o;
    number_densities[i][_CO]  = number_densities[i][_TOTAL] * mr_co;
    number_densities[i][_CO2] = number_densities[i][_TOTAL] * mr_co2;
    number_densities[i][_CH4] = number_densities[i][_TOTAL] * mr_ch4;
    number_densities[i][_NH3] = number_densities[i][_TOTAL] * mr_nh3;
    number_densities[i][_H2S] = number_densities[i][_TOTAL] * mr_h2s;

    number_densities[i][_Na] = number_densities[i][_TOTAL] * mr_na;
    number_densities[i][_K] = number_densities[i][_TOTAL] * mr_k;

    number_densities[i][_H2] = number_densities[i][_TOTAL] * mr_h2;
    number_densities[i][_He] = number_densities[i][_TOTAL] * mr_he;
  }


  //the mean molecular weight
  double mu = 18.*mr_h2o
            + 28.0101*mr_co 
            + 44.*mr_co2 
            + 16.04246*mr_ch4
            + 17.03052*mr_nh3
            + 34.08088*mr_h2s
            + 22.98977*mr_na
            + 39.0983*mr_k
            + 2.01588*mr_h2
            + 4.002602*mr_he;


  mean_molecular_weights.assign(nb_grid_points, mu);
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

  if (retrieval->config->useGPU())
    deleteFromDevice(absorption_coeff_gpu);


  delete radiative_transfer;

}



}

