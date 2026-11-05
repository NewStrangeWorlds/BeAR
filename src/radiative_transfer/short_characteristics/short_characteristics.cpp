// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "short_characteristics.h"

#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/quadrature.h"


namespace bear{


void ShortCharacteristics::calcSpectrum(
  const Atmosphere& atmosphere,
  const std::vector< std::vector<double> >& absorption_coeff,
  const std::vector< std::vector<double> >& scattering_coeff,
  const std::vector< std::vector<double> >& cloud_optical_depth,
  const std::vector< std::vector<double> >& cloud_single_scattering,
  const std::vector< std::vector<double> >& cloud_asym_param,
  const double spectrum_scaling,
  std::vector<double>& spectrum)
{

  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<spectrum.size(); ++i)
    spectrum[i] = calcSpectrum(absorption_coeff[i], cloud_optical_depth[i], atmosphere.temperature, atmosphere.altitude, i);


  for (auto & i : spectrum) i *= spectrum_scaling;
}



double ShortCharacteristics::calcSpectrum(
  const std::vector<double>& absorption_coeff,
  const std::vector<double>& cloud_optical_depth,
  const std::vector<double>& temperature,
  const std::vector<double>& vertical_grid,
  const size_t nu_index)
{
  size_t nb_grid_points = absorption_coeff.size();


  std::vector<double> optical_depth_layer(nb_grid_points-1, 0.0);

  for (size_t i=0; i<nb_grid_points-1; ++i)
    optical_depth_layer[i] = (vertical_grid[i+1] - vertical_grid[i]) * (absorption_coeff[i+1] + absorption_coeff[i])/2.;

  if (cloud_optical_depth.size() != 0)
    for (size_t i=0; i<nb_grid_points-1; ++i) optical_depth_layer[i] += cloud_optical_depth[i];
  
  
  std::vector<double> planck_function(nb_grid_points, 0.0);

  for (size_t i=0; i<nb_grid_points; ++i)
    planck_function[i] = aux::planckFunctionWavenumber(temperature[i], spectral_grid->wavenumber_list[nu_index]);


  //inner boundary condition
  double intensity_mu1 = planck_function[0];
  double intensity_mu2 = planck_function[0];


  for (size_t i=0; i<nb_grid_points-1; ++i)
  {
    if (optical_depth_layer[i] == 0)
      continue;

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

  const double flux = 
    2.0 * constants::pi 
    * (intensity_mu1 * gauss_nodes[0] * gauss_weights[0] 
     + intensity_mu2 * gauss_nodes[1] * gauss_weights[1]) * 1e-3; //in W m-2 cm-1


  return flux;
}


}
