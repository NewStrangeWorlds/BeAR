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

#include "adding_doubling.h"

#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/quadrature.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../additional/exceptions.h"
#include "../../CUDA_kernels/data_management_kernels.h"

#include <adding_doubling.hpp>


namespace bear{


AddingDoubling::AddingDoubling(
  SpectralGrid* spectral_grid_ptr,
  const size_t nb_quadrature_,
  const size_t nb_grid_points_,
  const bool use_gpu)
   : RadiativeTransfer(spectral_grid_ptr)
   , nb_quadrature(nb_quadrature_)
   , nb_grid_points(nb_grid_points_)
   , gpu_enabled(use_gpu)
{
  if (nb_quadrature < 2)
  {
    std::string error_message = "Adding-doubling RT requires at least 2 quadrature points\n";
    throw InvalidInput(std::string ("AddingDoubling::AddingDoubling"), error_message);
  }

  const size_t nb_layers = nb_grid_points - 1;

  if (!use_gpu)
  {
    const int nb_cores = omp_get_max_threads();

    adrt::ADConfig config_template(nb_layers, nb_quadrature);

    config_template.index_from_bottom = true;
    config_template.surface_albedo = 0;
    config_template.solar_flux = 0;
    config_template.solar_mu = 0.5;
    config_template.top_emission = 0;
    config_template.surface_emission = 0;

    // Allocate with thermal emission enabled so temperature array is always present
    config_template.use_thermal_emission = true;
    config_template.allocate();

    // Set constant optical properties
    for (size_t j = 0; j < nb_layers; ++j)
    {
      config_template.single_scat_albedo[j] = 0.0;
      config_template.setIsotropic(j);
    }

    configs.assign(nb_cores, config_template);
    workspaces.resize(nb_cores);
  }
}



// Destructor is defined in adding_doubling_kernels.cu
// because it needs to clean up CUDA resources



void AddingDoubling::calcSpectrum(
  const Atmosphere& atmosphere,
  const std::vector< std::vector<double> >& absorption_coeff,
  const std::vector< std::vector<double> >& scattering_coeff,
  const std::vector< std::vector<double> >& cloud_optical_depth,
  const std::vector< std::vector<double> >& cloud_single_scattering,
  const std::vector< std::vector<double> >& cloud_asym_param,
  const double spectrum_scaling,
  std::vector<double>& spectrum)
{
  const size_t nb_layers = nb_grid_points - 1;

  // Set per-call constants: temperature structure
  for (size_t t = 0; t < configs.size(); ++t)
  {
    for (size_t i = 0; i < atmosphere.temperature.size(); ++i)
      configs[t].temperature[i] = atmosphere.temperature[i];
  }


  std::string error_message = "";

  #pragma omp parallel for
  for (size_t i = 0; i < spectrum.size(); ++i)
  {
    if (!error_message.empty()) continue;

    const int thread = omp_get_thread_num();
    auto& config = configs[thread];

    // Set per-wavelength quantities: optical depth and single scattering albedo
    for (size_t j = 0; j < nb_layers; ++j)
    {
      const double dz = atmosphere.altitude[j+1] - atmosphere.altitude[j];

      const double abs_depth = dz
        * (absorption_coeff[i][j+1] + absorption_coeff[i][j]) / 2.;
      const double scat_depth = dz
        * (scattering_coeff[i][j+1] + scattering_coeff[i][j]) / 2.;

      double optical_depth = abs_depth + scat_depth;

      if (cloud_optical_depth[i].size() != 0)
        optical_depth += cloud_optical_depth[i][j];

      if (optical_depth < 0.0) optical_depth = 0.0;

      config.delta_tau[j] = optical_depth;
      config.single_scat_albedo[j] = (optical_depth > 0.0)
        ? scat_depth / optical_depth : 0.0;
    }

    // Set wavenumber for Planck function
    const double wavenumber = spectral_grid->wavenumber_list[i];
    config.wavenumber_low = wavenumber;
    config.wavenumber_high = wavenumber;

    // Enable thermal emission if Planck function is non-negligible
    config.use_thermal_emission =
      adrt::planckFunction(wavenumber, wavenumber, atmosphere.temperature[0]) > 1.e-35;

    try
    {
      // Solve and extract upward flux at TOA
      // With index_from_bottom=true, index 0 = BOA, last = TOA
      auto result = adrt::solve(config, workspaces[thread]);
      spectrum[i] = result.flux_up.back() * spectrum_scaling;
    }
    catch (const std::exception& e)
    {
      #pragma omp critical
      {
        if (error_message.empty())
          error_message = std::string("Adding-doubling solver error at spectral point ")
            + std::to_string(i) + ": " + e.what();
      }
    }
  }

  if (!error_message.empty())
    throw std::runtime_error(error_message);
}


}
