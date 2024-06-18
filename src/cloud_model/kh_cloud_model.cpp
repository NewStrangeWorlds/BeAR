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

#include "kh_cloud_model.h"
#include "grey_cloud_model.h"

#include "../forward_model/atmosphere/atmosphere.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/physical_const.h"
#include "../additional/exceptions.h"


namespace helios{


KHCloudModel::KHCloudModel(const std::vector<std::string>& parameters)
 : GreyCloudModel(parameters)
{
  //general case
  nb_parameters = 6;

  if (parameters.size() < 1)
  {
    std::string error_message = 
      "Expected at least one parameter for the KH non-grey cloud model, but only found " 
      + std::to_string(parameters.size()) + "\n";
    throw InvalidInput(std::string ("forward_model.config"), error_message);
  }

  reference_wavelength = std::stod(parameters[0]);

  //fixed bottom case
  if (parameters.size() == 2 && parameters[1] == "fb") 
  {
    fixed_bottom = true;
    nb_parameters = 5;
  }

}


//calculates the vertical distribution of the cloud layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
void KHCloudModel::opticalProperties(
  const std::vector<double>& parameters, 
  const Atmosphere& atmosphere,
  SpectralGrid* spectral_grid,
  std::vector<std::vector<double>>& optical_depth, 
  std::vector<std::vector<double>>& single_scattering, 
  std::vector<std::vector<double>>& asym_param)
{ 
  double cloud_optical_depth = parameters[0];
  double q0 = parameters[1];
  double a0 = parameters[2];
  double particle_size = parameters[3];
  double cloud_top_pressure = parameters[4];

  double reference_size_parameter = 2*constants::pi * particle_size / reference_wavelength;


  double cloud_bottom_pressure = 0;

  if (fixed_bottom == false)
    cloud_bottom_pressure = cloud_top_pressure * parameters[5];


  if (cloud_optical_depth < 0) cloud_optical_depth = 0;


  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  size_t nb_grid_points = atmosphere.nb_grid_points;


  //pre-allocate the data
  optical_depth.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));
  single_scattering.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));
  asym_param.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));


  unsigned int cloud_top_index = 0;
  unsigned int cloud_bottom_index = 0;

  if (fixed_bottom == true)
    cloudPosition(atmosphere, cloud_top_pressure, cloud_top_index, cloud_bottom_index);
  else
    cloudPosition(atmosphere, cloud_top_pressure, cloud_bottom_pressure, cloud_top_index, cloud_bottom_index);


  double optical_depth_layer_reference = cloud_optical_depth/static_cast<double>(cloud_top_index - cloud_bottom_index); 


  std::vector<double> size_parameter (nb_spectral_points, 0);

  for (size_t i=0; i<nb_spectral_points; ++i)
    size_parameter[i] = 2*constants::pi * particle_size / spectral_grid->wavelength_list[i];

  double reference_q = q0 * std::pow(reference_size_parameter, -a0) + std::pow(reference_size_parameter, 0.2);
  
  //the cloud here only has absorption, no scattering
  for (size_t j=0; j<nb_spectral_points; ++j)
    for (size_t i=cloud_top_index; i>cloud_bottom_index; --i)
    {
      optical_depth[j][i] = 
        optical_depth_layer_reference * reference_q / (q0 * std::pow(size_parameter[j], -a0) 
        + std::pow(size_parameter[j], 0.2) );
    }
}


}

