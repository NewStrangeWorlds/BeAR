// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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
#include "../additional/exceptions.h"


namespace bear{


GreyCloudModel::GreyCloudModel(const std::vector<std::string>& parameters)
{
  //optional model options:
  //  "fb"          fixed cloud bottom, one scale height below the cloud top
  //  "scattering"  the cloud scatters with a constant single scattering albedo
  //                and Henyey-Greenstein asymmetry parameter (two extra free parameters)
  for (auto & p : parameters)
  {
    if (p == "fb")
      fixed_bottom = true;
    else if (p == "scattering")
      scattering = true;
    else
    {
      std::string error_message =
        "Unknown option \"" + p + "\" for the grey cloud model (allowed: fb, scattering)\n";
      throw InvalidInput(std::string ("GreyCloudModel::GreyCloudModel"), error_message);
    }
  }

  //populate the parameter names in the same order readParameters() reads them
  parameter_names = {
    "cloud_optical_depth",
    "cloud_top_pressure"};

  if (fixed_bottom == false)
    parameter_names.push_back("cloud_bottom_fraction");

  if (scattering == true)
  {
    parameter_names.push_back("cloud_single_scattering_albedo");
    parameter_names.push_back("cloud_asym_param");
  }
}



//extracts the free parameters; the single scattering albedo and the asymmetry
//parameter are zero unless the scattering option is used
void GreyCloudModel::readParameters(
  const std::vector<double>& parameters,
  double& cloud_optical_depth,
  double& cloud_top_pressure,
  double& cloud_bottom_pressure,
  double& single_scattering_albedo,
  double& asym_param) const
{
  cloud_optical_depth = parameters[0];
  cloud_top_pressure = parameters[1];
  cloud_bottom_pressure = 0;

  size_t idx = 2;

  if (fixed_bottom == false)
    cloud_bottom_pressure = cloud_top_pressure * parameters[idx++];

  if (cloud_optical_depth < 0) cloud_optical_depth = 0;

  single_scattering_albedo = 0;
  asym_param = 0;

  if (scattering == true)
  {
    single_scattering_albedo = std::min(std::max(parameters[idx], 0.0), 1.0);
    asym_param = std::min(std::max(parameters[idx+1], -1.0), 1.0);
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
  double cloud_optical_depth = 0;
  double cloud_top_pressure = 0;
  double cloud_bottom_pressure = 0;
  double cloud_single_scattering_albedo = 0;
  double cloud_asym_param = 0;

  readParameters(
    parameters,
    cloud_optical_depth,
    cloud_top_pressure,
    cloud_bottom_pressure,
    cloud_single_scattering_albedo,
    cloud_asym_param);


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

  //constant optical properties in all cloud layers
  //without the scattering option, the cloud is a pure absorber
  for (size_t j=0; j<nb_spectral_points; ++j)
    for (size_t i=cloud_top_index; i>cloud_bottom_index; --i)
    {
      optical_depth[j][i] = optical_depth_layer;
      single_scattering[j][i] = cloud_single_scattering_albedo;
      asym_param[j][i] = cloud_asym_param;
    }
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


  if (top_pressure <= atmosphere.pressure.back())
    top_index = atmosphere.nb_grid_points - 1;

  if (bottom_pressure >= atmosphere.pressure[0])
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

