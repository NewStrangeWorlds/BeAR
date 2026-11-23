// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef OPACITY_CALC_H
#define OPACITY_CALC_H


#include <vector>
#include <memory>

#include "transport_coeff.h"
#include "../config/global_config.h"
#include "../spectral_grid/spectral_grid.h"
#include "../forward_model/atmosphere/atmosphere.h"
#include "../cloud_model/cloud_model.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../CUDA_kernels/cross_section_kernels.h"
#include "../CUDA_kernels/batched_opacity_kernels.h"


namespace bear{


class OpacityCalculation {
  public:
    OpacityCalculation(
      GlobalConfig* config_, 
      SpectralGrid* spectral_grid_,
      Atmosphere* atmosphere_,
      const std::vector<std::string> opacity_species_symbol,
      const std::vector<std::string> opacity_species_folder,
      const bool use_gpu_,
      const bool use_cloud_,
      const bool extinction_only_ = false)
      : spectral_grid(spectral_grid_)
      , atmosphere(atmosphere_)
      , transport_coeff(
          config_, 
          spectral_grid_, 
          opacity_species_symbol, 
          opacity_species_folder)
      , use_gpu(use_gpu_)
      , use_cloud(use_cloud_)
      , extinction_only(extinction_only_)
      {
        initDeviceMemory();
      }
    ~OpacityCalculation();

    void calculate(
      std::vector<std::unique_ptr<CloudModel>>& cloud_models,
      const std::vector<double>& cloud_parameter);
    void calculateGPU(
      std::vector<std::unique_ptr<CloudModel>>& cloud_models,
      const std::vector<double>& cloud_parameter);

    std::vector< std::vector<double> > absorption_coeff;
    std::vector< std::vector<double> > scattering_coeff;

    std::vector< std::vector<double> > cloud_optical_depths;
    std::vector< std::vector<double> > cloud_single_scattering;
    std::vector< std::vector<double> > cloud_asym_param;

    //coefficients on the GPU
    //in the extinction-only mode, absorption_coeff_gpu contains the total gas 
    //extinction (absorption and scattering) and scattering_coeff_dev as well as the
    //cloud single scattering albedos and asymmetry parameters are not used (nullptr)
    float* absorption_coeff_gpu = nullptr;
    float* scattering_coeff_dev = nullptr;

    float* cloud_optical_depths_dev = nullptr;
    float* cloud_single_scattering_dev = nullptr;
    float* cloud_asym_param_dev = nullptr;
    
  private:
    SpectralGrid* spectral_grid;
    Atmosphere* atmosphere;
    TransportCoefficients transport_coeff;

    bool use_gpu = false;
    bool use_cloud = false;
    //GPU only: for models that need no separate scattering properties (transmission spectra)
    bool extinction_only = false;

    BatchedDeviceBuffers batch_buffers;

    void initDeviceMemory();
};



inline void OpacityCalculation::calculateGPU(
  std::vector<std::unique_ptr<CloudModel>>& cloud_models,
  const std::vector<double>& cloud_parameter)
{
  const size_t nb_grid_points = atmosphere->nb_grid_points;
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  // Collect all cross-section and Rayleigh metadata on the CPU
  std::vector<float*> cs1_ptrs, cs2_ptrs, cs3_ptrs, cs4_ptrs;
  std::vector<float> temp_factors, pres_factors;
  std::vector<float> cs_log_number_densities;
  std::vector<int> cs_grid_points;
  std::vector<float*> ray_ptrs;
  std::vector<double> ray_number_densities;
  std::vector<int> ray_grid_points;

  transport_coeff.prepareBatchedGPU(
    *atmosphere,
    cs1_ptrs, cs2_ptrs, cs3_ptrs, cs4_ptrs,
    temp_factors, pres_factors,
    cs_log_number_densities, cs_grid_points,
    ray_ptrs, ray_number_densities, ray_grid_points);

  //gas absorption and Rayleigh scattering in a single kernel launch
  //the kernel writes all elements of the coefficient arrays
  if (!cs_grid_points.empty() || !ray_grid_points.empty())
    launchBatchedOpacities(
      cs1_ptrs, cs2_ptrs, cs3_ptrs, cs4_ptrs,
      temp_factors, pres_factors,
      cs_log_number_densities, cs_grid_points,
      ray_ptrs, ray_number_densities, ray_grid_points,
      batch_buffers, nb_spectral_points, nb_grid_points,
      absorption_coeff_gpu,
      scattering_coeff_dev);
  else
  {
    initializeOnDevice(absorption_coeff_gpu, nb_spectral_points*nb_grid_points);
    initializeOnDevice(scattering_coeff_dev, nb_spectral_points*nb_grid_points);
  }

  //H- continuum, all grid points in a single batched launch
  transport_coeff.calculateContinuumGPU(
    atmosphere->temperature,
    atmosphere->number_densities,
    absorption_coeff_gpu);

  if (use_cloud)
  {
    const size_t nb_layers = nb_grid_points - 1;

    initializeOnDevice(cloud_optical_depths_dev, nb_layers*nb_spectral_points);
    initializeOnDevice(cloud_single_scattering_dev, nb_layers*nb_spectral_points);
    initializeOnDevice(cloud_asym_param_dev, nb_layers*nb_spectral_points);

    size_t nb_param = 0;

    for (auto & cm : cloud_models)
    {
      std::vector<double> parameter(
        cloud_parameter.begin() + nb_param,
        cloud_parameter.begin() + nb_param + cm->nbParameters());

      nb_param += cm->nbParameters();

      cm->opticalPropertiesGPU(
        parameter,
        *atmosphere,
        spectral_grid,
        cloud_optical_depths_dev,
        cloud_single_scattering_dev,
        cloud_asym_param_dev);
    }
  }
}


inline void OpacityCalculation::calculate(
  std::vector<std::unique_ptr<CloudModel>>& cloud_models,
  const std::vector<double>& cloud_parameter)
{ 
  const size_t nb_grid_points = atmosphere->nb_grid_points;
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  //calculate gas transport coefficients
  absorption_coeff.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));
  scattering_coeff.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    std::vector<double> absorption_coeff_level(nb_spectral_points, 0.0);
    std::vector<double> scattering_coeff_level(nb_spectral_points, 0.0);

    transport_coeff.calculate(
      atmosphere->temperature[i],
      atmosphere->pressure[i], 
      atmosphere->number_densities[i], 
      absorption_coeff_level,
      scattering_coeff_level);

    for (size_t j=0; j<nb_spectral_points; ++j)
    {
      absorption_coeff[j][i] = absorption_coeff_level[j];
      scattering_coeff[j][i] = scattering_coeff_level[j];
    }
  }

  cloud_optical_depths.assign(
    nb_spectral_points, 
    std::vector<double>(nb_grid_points-1, 0.0));
  cloud_single_scattering.assign(
    nb_spectral_points, 
    std::vector<double>(nb_grid_points-1, 0.0));
  cloud_asym_param.assign(
    nb_spectral_points, 
    std::vector<double>(nb_grid_points-1, 0.0));


  if (use_cloud)
  { 
    size_t nb_param = 0;

    for (auto & i : cloud_models)
    {
      std::vector<double> parameter(
        cloud_parameter.begin() + nb_param, 
        cloud_parameter.begin() + nb_param + i->nbParameters());

      nb_param += i->nbParameters();

      std::vector<std::vector<double>> optical_depths;
      std::vector<std::vector<double>> single_scattering;
      std::vector<std::vector<double>> asym_param;
      i->opticalProperties(
        parameter, 
        *atmosphere, 
        spectral_grid, 
        optical_depths, 
        single_scattering, 
        asym_param);

      for (size_t j=0; j<nb_spectral_points; ++j)
        for (size_t k=0; k<nb_grid_points-1; ++k)
        {
          //scattering optical depths of the already accumulated clouds and the new model
          const double tau_scattering1 = cloud_optical_depths[j][k] * cloud_single_scattering[j][k];
          const double tau_scattering2 = optical_depths[j][k] * single_scattering[j][k];
          const double tau_scattering_mix = tau_scattering1 + tau_scattering2;
          
          const double optical_depth_mix = cloud_optical_depths[j][k] + optical_depths[j][k];
          
          double single_scattering_mix = 0;
          
          if (optical_depth_mix > 0)
            single_scattering_mix = tau_scattering_mix / optical_depth_mix;

          //asymmetry parameter is weighted by the scattering optical depth;
          //without any scattering it is left at zero
          double mix_asym_param = 0;

          if (tau_scattering_mix > 0)
            mix_asym_param = (tau_scattering1 * cloud_asym_param[j][k]
                            + tau_scattering2 * asym_param[j][k]) / tau_scattering_mix;

          cloud_optical_depths[j][k] = optical_depth_mix;
          cloud_single_scattering[j][k] = single_scattering_mix;
          cloud_asym_param[j][k] = mix_asym_param;
        }
    }
  }


}



inline void OpacityCalculation::initDeviceMemory()
{ 
  if (!use_gpu) return;
  
  const size_t nb_grid_points = atmosphere->nb_grid_points;
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  allocateOnDevice(absorption_coeff_gpu, nb_grid_points*nb_spectral_points);

  if (!extinction_only)
    allocateOnDevice(scattering_coeff_dev, nb_grid_points*nb_spectral_points);

  if (use_cloud)
  { 
    const size_t nb_layers = nb_grid_points - 1;

    allocateOnDevice(cloud_optical_depths_dev, nb_layers*nb_spectral_points);

    if (!extinction_only)
    {
      allocateOnDevice(cloud_single_scattering_dev, nb_layers*nb_spectral_points);
      allocateOnDevice(cloud_asym_param_dev, nb_layers*nb_spectral_points);
    }

    initializeOnDevice(cloud_optical_depths_dev, nb_layers*nb_spectral_points);
    initializeOnDevice(cloud_single_scattering_dev, nb_layers*nb_spectral_points);
    initializeOnDevice(cloud_asym_param_dev, nb_layers*nb_spectral_points);
  }
}



inline OpacityCalculation::~OpacityCalculation()
{
  if (use_gpu)
  {
    deleteFromDevice(absorption_coeff_gpu);

    if (scattering_coeff_dev != nullptr)
      deleteFromDevice(scattering_coeff_dev);

    freeBatchBuffers(batch_buffers);

    if (use_cloud)
    {
      deleteFromDevice(cloud_optical_depths_dev);

      if (cloud_single_scattering_dev != nullptr)
        deleteFromDevice(cloud_single_scattering_dev);

      if (cloud_asym_param_dev != nullptr)
        deleteFromDevice(cloud_asym_param_dev);
    }
  }
}


}

#endif
