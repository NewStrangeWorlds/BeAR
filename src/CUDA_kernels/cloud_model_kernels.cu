// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "../cloud_model/grey_cloud_model.h"
#include "../cloud_model/kh_cloud_model.h"
#include "../cloud_model/power_law_cloud_model.h"
#include "../spectral_grid/spectral_grid.h"

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "error_check.h"
#include "reduce_kernels.h"
#include "data_management_kernels.h"
#include "../additional/physical_const.h"


namespace bear{


  __global__ void  convertOpticalDepthDev(
  float* optical_depth,
  float* altitude,
  const int nb_grid_points,
  const int nb_spectral_points,
  float* extinction_coeff)
{
  //the thread index tid is the wavelength
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    for (int i=1; i<nb_grid_points; ++i)
    { 
      const int j = i*nb_spectral_points + tid;
      const int jm = (i-1)*nb_spectral_points + tid;
      const double delta_z = altitude[i] - altitude[i-1];
      
      //convert optical depth to extinction coefficient
      //uses a finite difference to approximate the derivative
      extinction_coeff[j] = optical_depth[jm] / delta_z;
    }
  }

}


//mixes the optical properties of a new cloud layer into the already accumulated
//cloud properties at array index j
//the single scattering albedo is weighted by the optical depth, the asymmetry
//parameter by the scattering optical depth (same scheme as OpacityCalculation::calculate)
__device__ __forceinline__ void mixCloudLayer(
  float* optical_depth,
  float* single_scattering_albedo,
  float* asym_param,
  const int j,
  const double layer_optical_depth,
  const double layer_single_scattering,
  const double layer_asym_param)
{
  const double tau_scattering1 = optical_depth[j] * single_scattering_albedo[j];
  const double tau_scattering2 = layer_optical_depth * layer_single_scattering;
  const double tau_scattering_mix = tau_scattering1 + tau_scattering2;

  const double optical_depth_mix = optical_depth[j] + layer_optical_depth;

  double single_scattering_mix = 0;

  if (optical_depth_mix > 0)
    single_scattering_mix = tau_scattering_mix / optical_depth_mix;

  double asym_param_mix = 0;

  if (tau_scattering_mix > 0)
    asym_param_mix = (tau_scattering1 * asym_param[j]
                    + tau_scattering2 * layer_asym_param) / tau_scattering_mix;

  optical_depth[j] = optical_depth_mix;
  single_scattering_albedo[j] = single_scattering_mix;
  asym_param[j] = asym_param_mix;
}


__global__ void greyCloudModel(
  float* optical_depth, 
  float* single_scattering_albedo, 
  float* asym_param,
  const double layer_optical_depth,
  const int cloud_top,
  const int cloud_bottom,
  const int nb_spectral_points)
{
  //the thread index tid is the wavelength
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    for (int i=cloud_top; i>cloud_bottom; --i)
    { 
      const int j = i*nb_spectral_points + tid;

      //the grey cloud here only has absorption, no scattering
      mixCloudLayer(
        optical_depth, single_scattering_albedo, asym_param,
        j, layer_optical_depth, 0.0, 0.0);
    }

  }

}



__global__ void KHnongreyCloudModel(
  float* optical_depth, 
  float* single_scattering_albedo, 
  float* asym_param,
  double* wavelengths,
  const double optical_depth_layer_ref,
  const double reference_q,
  const double particle_size,
  const double a0,
  const double q0,
  const int cloud_top,
  const int cloud_bottom,
  const int nb_spectral_points)
{
  //the thread index tid is the wavelength
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    const double size_parameter = 2*constants::pi * particle_size / wavelengths[tid];

    const double optical_depth_layer = optical_depth_layer_ref * reference_q 
      / (q0 * pow(size_parameter, -a0) + pow(size_parameter, 0.2) );

    for (int i=cloud_top; i>cloud_bottom; --i)
    {
      const int j = i*nb_spectral_points + tid;

      //the cloud here only has absorption, no scattering
      mixCloudLayer(
        optical_depth, single_scattering_albedo, asym_param,
        j, optical_depth_layer, 0.0, 0.0);
    }
  }

}



__global__ void powerLawCloudModel(
  float* optical_depth, 
  float* single_scattering_albedo, 
  float* asym_param,
  double* wavelengths,
  const double optical_depth_layer_ref,
  const double reference_value,
  const double exponent,
  const int cloud_top,
  const int cloud_bottom,
  const int nb_spectral_points)
{
  //the thread index tid is the wavelength
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    const double optical_depth_layer = optical_depth_layer_ref * std::pow(wavelengths[tid], exponent) 
                                      / reference_value;

    for (int i=cloud_top; i>cloud_bottom; --i)
    {
      const int j = i*nb_spectral_points + tid;

      //the cloud here only has absorption, no scattering
      mixCloudLayer(
        optical_depth, single_scattering_albedo, asym_param,
        j, optical_depth_layer, 0.0, 0.0);
    }
  }

}



//converts layer optical depths to level-based extinction coefficients
__host__ void CloudModel::convertOpticalDepthGPU(
  float* optical_depth_dev,
  float* altitude,
  const size_t nb_grid_points,
  const size_t nb_spectral_points,
  float* extinction_coeff_dev)
{

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  convertOpticalDepthDev<<<blocks,threads>>>(
    optical_depth_dev,
    altitude,
    nb_grid_points,
    nb_spectral_points,
    extinction_coeff_dev);

  CUDA_CHECK_AFTER_KERNEL();
}



//calculates the vertical distribution of the grey layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
__host__ void GreyCloudModel::opticalPropertiesGPU(const std::vector<double>& parameters, const Atmosphere& atmosphere,
  SpectralGrid* spectral_grid,
  float* optical_depth_dev, 
  float* single_scattering_dev, 
  float* asym_param)
{
  double cloud_optical_depth = parameters[0];
  double cloud_top_pressure = parameters[1];
  double cloud_bottom_pressure = 0;

  if (fixed_bottom == false)
    cloud_bottom_pressure = cloud_top_pressure * parameters[2];


  if (cloud_optical_depth < 0) cloud_optical_depth = 0;


  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  unsigned int cloud_top_index = 0;
  unsigned int cloud_bottom_index = 0;

  if (fixed_bottom == true)
    cloudPosition(atmosphere, cloud_top_pressure, cloud_top_index, cloud_bottom_index);
  else
    cloudPosition(atmosphere, cloud_top_pressure, cloud_bottom_pressure, cloud_top_index, cloud_bottom_index);


  const double optical_depth_layer = cloud_optical_depth/static_cast<double>(cloud_top_index - cloud_bottom_index);

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  greyCloudModel<<<blocks,threads>>>(
    optical_depth_dev,
    single_scattering_dev,
    asym_param,
    optical_depth_layer,
    cloud_top_index,
    cloud_bottom_index,
    nb_spectral_points);

  CUDA_CHECK_AFTER_KERNEL(); 
}



//calculates the vertical distribution of the KH non-grey layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
__host__ void KHCloudModel::opticalPropertiesGPU(
  const std::vector<double>& parameters, 
  const Atmosphere& atmosphere,
  SpectralGrid* spectral_grid,
  float* optical_depth_dev, 
  float* single_scattering_dev, 
  float* asym_param)
{
  double cloud_optical_depth = parameters[0];
  double q0 = parameters[1];
  double a0 = parameters[2];
  double particle_size = parameters[3];
  double cloud_top_pressure = parameters[4];


  double cloud_bottom_pressure = 0;

  if (fixed_bottom == false)
    cloud_bottom_pressure = cloud_top_pressure * parameters[5];


  if (cloud_optical_depth < 0) cloud_optical_depth = 0;


  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  unsigned int cloud_top_index = 0;
  unsigned int cloud_bottom_index = 0;

  if (fixed_bottom == true)
    cloudPosition(atmosphere, cloud_top_pressure, cloud_top_index, cloud_bottom_index);
  else
    cloudPosition(atmosphere, cloud_top_pressure, cloud_bottom_pressure, cloud_top_index, cloud_bottom_index);


  double optical_depth_layer_reference = cloud_optical_depth/static_cast<double>(cloud_top_index - cloud_bottom_index);

  double reference_size_parameter = 2*constants::pi * particle_size / reference_wavelength;
  double reference_q = q0 * std::pow(reference_size_parameter, -a0) + std::pow(reference_size_parameter, 0.2);

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  KHnongreyCloudModel<<<blocks,threads>>>(
    optical_depth_dev,
    single_scattering_dev,
    asym_param,
    spectral_grid->wavelength_list_gpu,
    optical_depth_layer_reference,
    reference_q,
    particle_size,
    a0,
    q0,
    cloud_top_index,
    cloud_bottom_index,
    nb_spectral_points);

  CUDA_CHECK_AFTER_KERNEL();
}


//calculates the vertical distribution of the power-law cloud layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
__host__ void PowerLawCloudModel::opticalPropertiesGPU(
  const std::vector<double>& parameters, 
  const Atmosphere& atmosphere,
  SpectralGrid* spectral_grid,
  float* optical_depth_dev, 
  float* single_scattering_dev, 
  float* asym_param)
{
  double cloud_optical_depth = parameters[0];
  double exponent = parameters[1];
  double cloud_top_pressure = parameters[2];

  double cloud_bottom_pressure = 0;

  if (fixed_bottom == false)
    cloud_bottom_pressure = cloud_top_pressure * parameters[3];


  if (cloud_optical_depth < 0) cloud_optical_depth = 0;


  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  unsigned int cloud_top_index = 0;
  unsigned int cloud_bottom_index = 0;

  if (fixed_bottom == true)
    cloudPosition(atmosphere, cloud_top_pressure, cloud_top_index, cloud_bottom_index);
  else
    cloudPosition(atmosphere, cloud_top_pressure, cloud_bottom_pressure, cloud_top_index, cloud_bottom_index);


  const double optical_depth_layer_reference = cloud_optical_depth
    /static_cast<double>(cloud_top_index - cloud_bottom_index); 
  const double reference_value =  std::pow(reference_wavelength, exponent);
  
  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  powerLawCloudModel<<<blocks,threads>>>(
    optical_depth_dev,
    single_scattering_dev,
    asym_param,
    spectral_grid->wavelength_list_gpu,
    optical_depth_layer_reference,
    reference_value,
    exponent,
    cloud_top_index,
    cloud_bottom_index,
    nb_spectral_points);

  CUDA_CHECK_AFTER_KERNEL();
}


}
