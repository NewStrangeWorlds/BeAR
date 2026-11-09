// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "adding_doubling.h"

#include "../layer_optical_properties.h"
#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/quadrature.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../additional/exceptions.h"
#include "../../CUDA_kernels/data_management_kernels.h"

#include <adding_doubling.hpp>


namespace bear{


// Legendre moments of the mixed gas + cloud phase function of a layer, weighted by
// the scattering optical depths. Gas: Rayleigh (chi_0 = 1, chi_2 = 0.1),
// cloud: Henyey-Greenstein (chi_l = g^l). The (2l+1) factor is applied by the solver.
static void setPhaseMoments(
  std::vector<double>& moments,
  const double cloud_weight,
  const double asym_param)
{
  const double gas_weight = 1.0 - cloud_weight;
  double gk = 1.0;

  for (size_t m = 0; m < moments.size(); ++m)
  {
    double rayleigh = 0.0;
    if (m == 0) rayleigh = 1.0;
    if (m == 2) rayleigh = 0.1;

    moments[m] = gas_weight * rayleigh + cloud_weight * gk;
    gk *= asym_param;
  }
}



AddingDoubling::AddingDoubling(
  SpectralGrid* spectral_grid_ptr,
  const size_t nb_quadrature_,
  const size_t nb_grid_points_,
  const bool use_gpu,
  const bool use_delta_m_)
   : RadiativeTransfer(spectral_grid_ptr)
   , nb_quadrature(nb_quadrature_)
   , nb_grid_points(nb_grid_points_)
   , gpu_enabled(use_gpu)
   , use_delta_m(use_delta_m_)
{
  if (nb_quadrature < 2)
  {
    std::string error_message = "Adding-doubling RT requires at least 2 quadrature points\n";
    throw InvalidInput(std::string ("AddingDoubling::AddingDoubling"), error_message);
  }

  const size_t nb_layers = nb_grid_points - 1;

  // delta-M scaling needs the 2N-th moment in addition to the 2N used by the solver
  nb_phase_moments = use_delta_m ? 2 * nb_quadrature + 1 : 2 * nb_quadrature;

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
    config_template.use_delta_m = use_delta_m;
    config_template.allocate();

    // Default optical properties: non-scattering layers with the gas (Rayleigh) phase function
    for (size_t j = 0; j < nb_layers; ++j)
    {
      config_template.single_scat_albedo[j] = 0.0;
      config_template.setRayleigh(j);
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
      const LayerOpticalProperties p = layerOpticalProperties(
        atmosphere.altitude,
        absorption_coeff[i],
        scattering_coeff[i],
        cloud_optical_depth[i],
        cloud_single_scattering[i],
        cloud_asym_param[i],
        j);

      config.delta_tau[j] = p.optical_depth;
      config.single_scat_albedo[j] = p.single_scattering;

      // Phase function: gas Rayleigh + cloud Henyey-Greenstein, weighted by
      // their scattering optical depths
      setPhaseMoments(config.phase_function_moments[j], p.cloud_scattering_fraction, p.cloud_asym_param);
    }

    // Set wavenumber for Planck function
    const double wavenumber = spectral_grid->wavenumber_list[i];
    config.wavenumber_low = wavenumber;
    config.wavenumber_high = wavenumber;

    // Enable thermal emission if Planck function is non-negligible
    config.use_thermal_emission =
      adrt::planckFunction(wavenumber, wavenumber, atmosphere.temperature[0]) > 1.e-35;

    // Cold space at the top: no thermal downwelling at TOA.
    // The library only accepts top_temperature >= 0 together with thermal emission
    config.top_temperature = config.use_thermal_emission ? 0.0 : -1.0;

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
