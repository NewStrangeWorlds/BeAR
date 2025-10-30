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


#include "fixed_cloud_model.h"


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>


#include "../forward_model/atmosphere/atmosphere.h"
#include "../CUDA_kernels/data_management_kernels.h"


namespace bear{


FixedCloudModel::FixedCloudModel(
  const std::vector<std::vector<double>>& optical_depth,
  const std::vector<std::vector<double>>& single_scattering_albedo,
  const std::vector<std::vector<double>>& asymmetry_parameter)
{
  nb_parameters = 0;

  size_t nb_spectral_points = optical_depth[0].size();
  size_t nb_layers = optical_depth.size();
  
  fixed_optical_depth.assign(nb_spectral_points, std::vector<double>(nb_layers, 0.0));
  fixed_single_scattering_albedo.assign(nb_spectral_points, std::vector<double>(nb_layers, 0.0));
  fixed_asymmetry_parameter.assign(nb_spectral_points, std::vector<double>(nb_layers, 0.0));
  
  for (size_t i = 0; i < nb_spectral_points; i++)
    for (size_t j = 0; j < nb_layers; j++)
    {
      fixed_optical_depth[i][j] = optical_depth[j][i];
      fixed_single_scattering_albedo[i][j] = single_scattering_albedo[j][i];
      fixed_asymmetry_parameter[i][j] = asymmetry_parameter[j][i];
    }
}


//calculates the vertical distribution of the grey layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
void FixedCloudModel::opticalProperties(
  const std::vector<double>& parameters, 
  const Atmosphere& atmosphere,
  SpectralGrid* spectral_grid,
  std::vector<std::vector<double>>& optical_depth, 
  std::vector<std::vector<double>>& single_scattering, 
  std::vector<std::vector<double>>& asym_param)
{
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  size_t nb_layers = atmosphere.nb_grid_points - 1;

  //pre-allocate the data
  optical_depth.assign(nb_spectral_points, std::vector<double>(nb_layers, 0.0));
  single_scattering.assign(nb_spectral_points, std::vector<double>(nb_layers, 0.0));
  asym_param.assign(nb_spectral_points, std::vector<double>(nb_layers, 0.0));


  optical_depth = fixed_optical_depth;
  single_scattering = fixed_single_scattering_albedo;
  asym_param = fixed_asymmetry_parameter;
}


void FixedCloudModel::opticalPropertiesGPU(
  const std::vector<double>& parameters, 
  const Atmosphere& atmosphere,
  SpectralGrid* spectral_grid,
  double* optical_depth_dev, 
  double* single_scattering_dev, 
  double* asym_param)
{
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  size_t nb_layers = atmosphere.nb_grid_points - 1;

  std::vector<double> optical_depth_flat(nb_spectral_points*nb_layers, 0.0);
  std::vector<double> single_scattering_flat(nb_spectral_points*nb_layers, 0.0);
  std::vector<double> asymmetry_param_flat(nb_spectral_points*nb_layers, 0.0);
  
  for (size_t j = 0; j<nb_spectral_points; ++j)
  {
    for (size_t i=0; i<nb_layers; ++i)
    {
      const int idx = i*nb_spectral_points + j;

      optical_depth_flat[idx] = fixed_optical_depth[j][i];
      single_scattering_flat[idx] = fixed_optical_depth[j][i];
      asymmetry_param_flat[idx] = fixed_asymmetry_parameter[j][i];
    }
  }

  moveToDevice(optical_depth_dev, optical_depth_flat, false);
  moveToDevice(single_scattering_dev, single_scattering_flat, false);
  moveToDevice(asym_param, asymmetry_param_flat, false);
}


}
