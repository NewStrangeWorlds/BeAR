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


#include "short_characteristics.h"


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "../retrieval/retrieval.h"

#include "../additional/aux_functions.h"
#include "../additional/physical_const.h"
#include "../additional/quadrature.h"

#include "../CUDA_kernels/data_management_kernels.h"
#include "../CUDA_kernels/short_characteristics_kernels.h"




namespace helios{


void ShortCharacteristics::calcSpectrum(const std::vector< std::vector<double> >& absorption_coeff, 
                                        const std::vector< std::vector<double> >& scattering_coeff,
                                        const std::vector<double>& temperature, const std::vector<double>& vertical_grid,
                                        std::vector<double>& spectrum)
{

  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<spectrum.size(); ++i)
    spectrum[i] = calcSpectrum(absorption_coeff[i], temperature, vertical_grid, i);

}



double ShortCharacteristics::calcSpectrum(const std::vector<double>& absorption_coeff, const std::vector<double>& temperature, const std::vector<double>& vertical_grid,
                                          const size_t nu_index)
{
  size_t nb_grid_points = absorption_coeff.size();


  std::vector<double> optical_depth_layer(nb_grid_points-1, 0.0);

  for (size_t i=0; i<nb_grid_points-1; ++i)
    optical_depth_layer[i] = (vertical_grid[i+1] - vertical_grid[i]) * (absorption_coeff[i+1] + absorption_coeff[i])/2.;

  
  std::vector<double> planck_function(nb_grid_points, 0.0);

  for (size_t i=0; i<nb_grid_points; ++i)
    planck_function[i] = aux::planckFunctionWavenumber(temperature[i], retrieval->spectral_grid.wavenumber_list[nu_index]);


  //inner boundary condition
  double intensity_mu1 = planck_function[0];
  double intensity_mu2 = planck_function[0];


  for (size_t i=0; i<nb_grid_points-1; ++i)
  {
    const double delta1 = optical_depth_layer[i]/gauss_nodes[0];
    const double attenuation_factor1 = std::exp(-delta1);

    intensity_mu1 = intensity_mu1 * attenuation_factor1;

    const double beta1 = 1.0 + (attenuation_factor1 - 1.0)/delta1;
    const double gamma1 = -attenuation_factor1 - (attenuation_factor1 - 1.0)/delta1;

    intensity_mu1 += beta1 * planck_function[i+1] + gamma1 * planck_function[i];



    const double delta2 = optical_depth_layer[i]/gauss_nodes[1];
    const double attenuation_factor2 = std::exp(-delta2);

    intensity_mu2 = intensity_mu2 * attenuation_factor2;

    const double beta2 = 1.0 + (attenuation_factor2 - 1.0)/delta2;
    const double gamma2 = -attenuation_factor2 - (attenuation_factor2 - 1.0)/delta2;

    intensity_mu2 += beta2 * planck_function[i+1] + gamma2 * planck_function[i];
  }

  
  const double flux = 2.0 * helios::CONST_PI * (intensity_mu1 * gauss_nodes[0] * gauss_weights[0] + intensity_mu2 * gauss_nodes[1] * gauss_weights[1]) * 1e-3; //in W m-2 cm-1


  return flux;
}



void ShortCharacteristics::calcSpectrumGPU(double* model_spectrum_dev,
                                           double* absorption_coeff_dev, 
                                           double* scattering_coeff_dev, 
                                           double* wavenumber_list_dev,
                                           const std::vector<double>& cloud_optical_depth,
                                           const std::vector<double>& temperature, 
                                           const std::vector<double>& vertical_grid,
                                           const double radius_distance_scaling)
{

  if (cloud_optical_depth.size() != 0)
    helios::shortCharacteristicsGPU(model_spectrum_dev,
                                    absorption_coeff_dev, 
                                    retrieval->spectral_grid.wavenumber_list_gpu,
                                    cloud_optical_depth,
                                    temperature, vertical_grid,
                                    radius_distance_scaling,
                                    retrieval->spectral_grid.nbSpectralPoints());
  else
    helios::shortCharacteristicsGPU(model_spectrum_dev,
                                    absorption_coeff_dev, 
                                    retrieval->spectral_grid.wavenumber_list_gpu,
                                    temperature, vertical_grid,
                                    radius_distance_scaling,
                                    retrieval->spectral_grid.nbSpectralPoints());

}



}
