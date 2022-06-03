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


#include "grey_cloud_model.h"


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>


#include "../forward_model/atmosphere/atmosphere.h"
#include "../CUDA_kernels/data_management_kernels.h"


namespace helios{



//calculates the vertical distribution of the grey layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
void GreyCloudModel::opticalProperties(const std::vector<double>& parameters, const Atmosphere& atmosphere,
                                       SpectralGrid* spectral_grid,
                                       std::vector<std::vector<double>>& optical_depth, 
                                       std::vector<std::vector<double>>& single_scattering, 
                                       std::vector<std::vector<double>>& asym_param)
{
  double cloud_top_pressure = parameters[0];
  double cloud_bottom_pressure = cloud_top_pressure * parameters[1];
  double cloud_optical_depth = parameters[2];


  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  size_t nb_grid_points = atmosphere.nb_grid_points;


  //pre-allocate the data
  optical_depth.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));
  single_scattering.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));
  asym_param.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));


  unsigned int cloud_top_index = 0;
  unsigned int cloud_bottom_index = 0;

  cloudPosition(atmosphere, cloud_top_pressure, cloud_bottom_pressure, cloud_top_index, cloud_bottom_index);


  double optical_depth_layer = cloud_optical_depth/static_cast<double>(cloud_top_index - cloud_bottom_index); 

  //the grey cloud here only has absorption, no scattering
  for (size_t j=0; j<nb_spectral_points; ++j)
    for (size_t i=cloud_bottom_index; i<cloud_top_index; ++i)
      optical_depth[j][i] = optical_depth_layer;
}




//calculates the vertical distribution of the grey layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
void GreyCloudModel::opticalPropertiesGPU(const std::vector<double>& parameters, const Atmosphere& atmosphere,
                                          SpectralGrid* spectral_grid,
                                          double* optical_depth_dev, 
                                          double* single_scattering_dev, 
                                          double* asym_param)
{
  double cloud_top_pressure = parameters[0];
  double cloud_bottom_pressure = cloud_top_pressure * parameters[1];
  double cloud_optical_depth = parameters[2];


  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  size_t nb_grid_points = atmosphere.nb_grid_points;

  unsigned int cloud_top_index = 0;
  unsigned int cloud_bottom_index = 0;

  cloudPosition(atmosphere, cloud_top_pressure, cloud_bottom_pressure, cloud_top_index, cloud_bottom_index);


  double optical_depth_layer = cloud_optical_depth/static_cast<double>(cloud_top_index - cloud_bottom_index);

  std::vector<double> optical_depth(nb_spectral_points*(nb_grid_points-1), 0.0); 

  //the grey cloud here only has absorption, no scattering
  for (size_t j=0; j<nb_spectral_points; ++j)
    for (size_t i=cloud_bottom_index; i<cloud_top_index; ++i)
    {
      unsigned int index = i*nb_spectral_points + j;

      optical_depth[index] = optical_depth_layer;
      std::cout << j << "\t" << i << "\t" << index << "\t" << optical_depth[index] << "\n";
    }

  moveToDevice(optical_depth_dev, optical_depth, false);
}





//calculates the upper and lower grid point of the cloud based on the top and bottom pressure
void GreyCloudModel::cloudPosition(const Atmosphere& atmosphere, const double top_pressure, const double bottom_pressure, 
                                   unsigned int& top_index, unsigned int& bottom_index)
{
  for (size_t i=0; i<atmosphere.nb_grid_points; ++i)
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



}

