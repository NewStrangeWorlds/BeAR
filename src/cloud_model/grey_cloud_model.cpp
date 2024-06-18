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


GreyCloudModel::GreyCloudModel(const std::vector<std::string>& parameters)
{
  //general case
  nb_parameters = 3;

  //fixed bottom case
  if (parameters.size() == 1 && parameters[0] == "fb") 
  {
    fixed_bottom = true;
    nb_parameters = 2;
  }

}



//calculates the vertical distribution of the grey layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
void GreyCloudModel::opticalProperties(
  const std::vector<double>& parameters, 
  const Atmosphere& atmosphere,
  SpectralGrid* spectral_grid,
  std::vector<std::vector<double>>& optical_depth, 
  std::vector<std::vector<double>>& single_scattering, 
  std::vector<std::vector<double>>& asym_param)
{
  double cloud_optical_depth = parameters[0];
  double cloud_top_pressure = parameters[1];
  double cloud_bottom_pressure = 0;

  if (fixed_bottom == false)
    cloud_bottom_pressure = cloud_top_pressure * parameters[2];


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
    cloudPosition(
      atmosphere, cloud_top_pressure, cloud_top_index, cloud_bottom_index);
  else
    cloudPosition(
      atmosphere, cloud_top_pressure, cloud_bottom_pressure, cloud_top_index, cloud_bottom_index);


  double optical_depth_layer = 
    cloud_optical_depth/static_cast<double>(cloud_top_index - cloud_bottom_index); 

  //the grey cloud here only has absorption, no scattering
  for (size_t j=0; j<nb_spectral_points; ++j)
    for (size_t i=cloud_top_index; i>cloud_bottom_index; --i)
      optical_depth[j][i] = optical_depth_layer;
}




//calculates the upper and lower grid point of the cloud based on the top and bottom pressure
void GreyCloudModel::cloudPosition(
  const Atmosphere& atmosphere, 
  const double top_pressure, 
  const double bottom_pressure, 
  unsigned int& top_index, 
  unsigned int& bottom_index)
{
  for (size_t i=0; i<atmosphere.nb_grid_points-1; ++i)
  {
    if ((atmosphere.pressure[i] > top_pressure && atmosphere.pressure[i+1] < top_pressure) 
      || atmosphere.pressure[i] == top_pressure)
      top_index = i;

    if ((atmosphere.pressure[i] > bottom_pressure && atmosphere.pressure[i+1] < bottom_pressure) 
      || atmosphere.pressure[i] == bottom_pressure)
      bottom_index = i;
  }


  if (top_pressure < atmosphere.pressure.back())
    top_index = atmosphere.nb_grid_points - 1;

  if (bottom_pressure > atmosphere.pressure[0])
    bottom_index = 0;


  //clouds needs to occupy at least an entire atmospheric layer
  if (top_index == bottom_index)
    bottom_index -= 2;

  if (top_index < 2)
  {
    top_index = 2;
    bottom_index = 0;
  }
}


//calculates the upper and lower grid point of the cloud based on the top and bottom pressure
//this version assumes that the cloud extends one atmospheric scale height
void GreyCloudModel::cloudPosition(
  const Atmosphere& atmosphere,
  const double top_pressure, 
  unsigned int& top_index,
  unsigned int& bottom_index)
{ 
  //find the top pressure index first
  for (size_t i=0; i<atmosphere.nb_grid_points-1; ++i)
  {
    if ((atmosphere.pressure[i] > top_pressure && atmosphere.pressure[i+1] < top_pressure) 
      || atmosphere.pressure[i] == top_pressure)
      top_index = i;
  }


  if (top_pressure < atmosphere.pressure.back())
    top_index = atmosphere.nb_grid_points - 1;


  double bottom_altitude = atmosphere.altitude[top_index] - atmosphere.scale_height[top_index];
  if (bottom_altitude < 0) bottom_altitude = 0;


  for (size_t i=0; i<atmosphere.nb_grid_points-1; ++i)
  {
    if ((atmosphere.altitude[i] < bottom_altitude && atmosphere.altitude[i+1] > bottom_altitude) 
      || atmosphere.altitude[i] == bottom_altitude)
      bottom_index = i;
  }


  //clouds needs to occupy at least an entire atmospheric layer
  if (top_index == bottom_index)
    bottom_index -= 2;

  if (top_index == bottom_index)
    bottom_index -= 2;

  if (top_index < 2)
  {
    top_index = 2;
    bottom_index = 0;
  }
}



}

