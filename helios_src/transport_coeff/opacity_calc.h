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


#ifndef OPACITY_CALC_H
#define OPACITY_CALC_H


#include <vector>

#include "transport_coeff.h"
#include "../config/global_config.h"
#include "../spectral_grid/spectral_grid.h"
#include "../forward_model/atmosphere/atmosphere.h"
#include "../cloud_model/cloud_model.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../CUDA_kernels/cross_section_kernels.h"


namespace helios{


class OpacityCalculation {
  public:
    OpacityCalculation(
      GlobalConfig* config_, 
      SpectralGrid* spectral_grid_,
      Atmosphere* atmosphere_,
      const std::vector<std::string> opacity_species_symbol,
      const std::vector<std::string> opacity_species_folder,
      const bool use_gpu_,
      const bool use_cloud_)
      : spectral_grid(spectral_grid_)
      , atmosphere(atmosphere_)
      , transport_coeff(
          config_, 
          spectral_grid_, 
          opacity_species_symbol, 
          opacity_species_folder)
      , use_gpu(use_gpu_)
      , use_cloud(use_cloud_)
      {
        initDeviceMemory();
      }
    ~OpacityCalculation();

    void calculate(
      std::vector<CloudModel*>& cloud_models,
      const std::vector<double>& cloud_parameter);
    void calculateGPU(
      std::vector<CloudModel*>& cloud_models,
      const std::vector<double>& cloud_parameter);

    std::vector< std::vector<double> > absorption_coeff;
    std::vector< std::vector<double> > scattering_coeff;

    std::vector< std::vector<double> > cloud_optical_depths;
    std::vector< std::vector<double> > cloud_single_scattering;
    std::vector< std::vector<double> > cloud_asym_param;

    //pointer to the array that holds the pointers to the coefficients on the GPU
    double* absorption_coeff_gpu = nullptr;
    double* scattering_coeff_dev = nullptr;

    double* cloud_optical_depths_dev = nullptr;
    double* cloud_single_scattering_dev = nullptr;
    double* cloud_asym_param_dev = nullptr;
    
  private:
    SpectralGrid* spectral_grid;
    Atmosphere* atmosphere;
    TransportCoefficients transport_coeff;

    bool use_gpu = false;
    bool use_cloud = false;

    void initDeviceMemory();
};



inline void OpacityCalculation::calculateGPU(
  std::vector<CloudModel*>& cloud_models,
  const std::vector<double>& cloud_parameter)
{
  const size_t nb_grid_points = atmosphere->nb_grid_points;
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  initCrossSectionsHost(
    nb_spectral_points*nb_grid_points, 
    absorption_coeff_gpu);

  initCrossSectionsHost(
    nb_spectral_points*nb_grid_points, 
    scattering_coeff_dev);

  for (size_t i=0; i<nb_grid_points; ++i)
    transport_coeff.calculateGPU(
      atmosphere->temperature[i], 
      atmosphere->pressure[i], 
      atmosphere->number_densities[i],
      nb_grid_points, 
      i,
      absorption_coeff_gpu, 
      scattering_coeff_dev);


  if (use_cloud)
  { 
    const size_t nb_layers = nb_grid_points - 1;

    intializeOnDevice(cloud_optical_depths_dev, nb_layers*nb_spectral_points);
    intializeOnDevice(cloud_single_scattering_dev, nb_layers*nb_spectral_points);
    intializeOnDevice(cloud_asym_param_dev, nb_layers*nb_spectral_points);

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
  std::vector<CloudModel*>& cloud_models,
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

  cloud_optical_depths.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));
  cloud_single_scattering.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));
  cloud_asym_param.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));


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
          const double tau_scattering1 = cloud_optical_depths[j][k] * cloud_single_scattering[j][k];
          const double tau_scattering2 = optical_depths[j][k] * cloud_single_scattering[j][k];
          const double tau_scattering_mix = tau_scattering1 + tau_scattering2;
          
          const double optical_depth_mix = cloud_optical_depths[j][k] + optical_depths[j][k];
          
          double single_scattering_mix = 0;
          
          if (optical_depth_mix > 0)
            single_scattering_mix = tau_scattering_mix / optical_depth_mix;

          const double mix_asym_param = tau_scattering1/tau_scattering_mix * cloud_asym_param[j][k]
                                      + tau_scattering2/tau_scattering_mix * asym_param[j][k];

          cloud_optical_depths[j][k] = optical_depth_mix;
          cloud_single_scattering[j][k] = single_scattering_mix;
          cloud_asym_param[j][k] = mix_asym_param;
        }
    }
  }


}



inline void OpacityCalculation::initDeviceMemory()
{ 
  const size_t nb_grid_points = atmosphere->nb_grid_points;
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  allocateOnDevice(absorption_coeff_gpu, nb_grid_points*nb_spectral_points);
  allocateOnDevice(scattering_coeff_dev, nb_grid_points*nb_spectral_points);

  if (use_cloud)
  { 
    const size_t nb_layers = nb_grid_points - 1;

    allocateOnDevice(cloud_optical_depths_dev, nb_layers*nb_spectral_points);
    allocateOnDevice(cloud_single_scattering_dev, nb_layers*nb_spectral_points);
    allocateOnDevice(cloud_asym_param_dev, nb_layers*nb_spectral_points);

    intializeOnDevice(cloud_optical_depths_dev, nb_layers*nb_spectral_points);
    intializeOnDevice(cloud_single_scattering_dev, nb_layers*nb_spectral_points);
    intializeOnDevice(cloud_asym_param_dev, nb_layers*nb_spectral_points);
  }
}



inline OpacityCalculation::~OpacityCalculation()
{
  if (use_gpu)
  {
    deleteFromDevice(absorption_coeff_gpu);
    deleteFromDevice(scattering_coeff_dev);

    if (use_cloud)
    {
      deleteFromDevice(cloud_optical_depths_dev);
      deleteFromDevice(cloud_single_scattering_dev);
      deleteFromDevice(cloud_asym_param_dev);
    }
  }     
}


}

#endif
