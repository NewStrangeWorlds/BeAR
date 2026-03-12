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

#include "discrete_ordinate.h"

#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/quadrature.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../additional/exceptions.h"

#include <DisortFluxConfig.hpp>
#include <FluxSolver.hpp>
#include <Planck.hpp>


namespace bear{


template<int NStr>
struct FluxSolverWrapper : FluxSolverBase {
  disortpp::DisortFluxSolver<NStr> solver;
  disortpp::FluxResult solve(disortpp::DisortFluxConfig& config) override {
    return solver.solve(config);
  }
};


namespace {

std::unique_ptr<FluxSolverBase> createSolver(size_t nb_streams)
{
  switch (nb_streams)
  {
    case 4:  return std::make_unique<FluxSolverWrapper<4>>();
    case 8:  return std::make_unique<FluxSolverWrapper<8>>();
    case 16: return std::make_unique<FluxSolverWrapper<16>>();
    case 32: return std::make_unique<FluxSolverWrapper<32>>();
    default:
    {
      std::string error_message =
        "DisORT flux solver only supports 4, 8, 16, or 32 streams. Got: "
        + std::to_string(nb_streams) + "\n";
      throw InvalidInput(std::string ("DiscreteOrdinates"), error_message);
    }
  }
}

} // anonymous namespace


DiscreteOrdinates::DiscreteOrdinates(
  SpectralGrid* spectral_grid_ptr,
  const size_t nb_streams_,
  const size_t nb_grid_points_,
  const bool use_gpu)
   : RadiativeTransfer(spectral_grid_ptr)
   , nb_streams(nb_streams_)
   , nb_grid_points(nb_grid_points_)
{
  if (use_gpu)
  {
    std::string error_message = "Radiative transfer model DisORT cannot run on the GPU\n";
    throw InvalidInput(std::string ("DiscreteOrdinates::DiscreteOrdinates"), error_message);
  }

  if (nb_streams < 4 || nb_streams % 2 != 0)
  {
    std::string error_message = "DisORT requires an even number of streams >= 4\n";
    throw InvalidInput(std::string ("DiscreteOrdinates::DiscreteOrdinates"), error_message);
  }

  const int nb_cores = omp_get_max_threads();

  solvers.resize(nb_cores);
  for (int i = 0; i < nb_cores; ++i)
    solvers[i] = createSolver(nb_streams);

  const size_t nb_layers = nb_grid_points - 1;

  disortpp::DisortFluxConfig config_template(nb_layers, nb_streams, nb_streams);

  config_template.index_from_bottom = true;

  config_template.direct_beam_flux = 0;
  config_template.direct_beam_mu = 0.5;
  config_template.surface_albedo = 0;
  config_template.isotropic_flux_top = 0;
  config_template.isotropic_flux_bottom = 0;
  config_template.temperature_top = 0;
  config_template.emissivity_top = 0;

  // Allocate with thermal emission enabled so the temperature array
  // is always present; the flag is toggled per wavelength
  config_template.use_thermal_emission = true;
  config_template.allocate();

  // Set single scattering albedo and phase function to isotropic (constant)
  for (size_t j = 0; j < nb_layers; ++j)
  {
    config_template.single_scat_albedo[j] = 0.0;
    config_template.setIsotropic(j);
  }

  configs.assign(nb_cores, config_template);
}



void DiscreteOrdinates::calcSpectrum(
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

  // Set per-call constants: temperature structure and bottom boundary temperature
  for (size_t t = 0; t < configs.size(); ++t)
  {
    configs[t].temperature_bottom = atmosphere.temperature[0];

    for (size_t i = 0; i < atmosphere.temperature.size(); ++i)
      configs[t].temperature[i] = atmosphere.temperature[i];
  }


  #pragma omp parallel for
  for (size_t i = 0; i < spectrum.size(); ++i)
  {
    const int thread = omp_get_thread_num();
    auto& config = configs[thread];

    // Set per-wavelength quantities: optical depth
    for (size_t j = 0; j < nb_layers; ++j)
    {
      double optical_depth = (atmosphere.altitude[j+1] - atmosphere.altitude[j])
                             * (absorption_coeff[i][j+1] + absorption_coeff[i][j]) / 2.;

      if (cloud_optical_depth[i].size() != 0)
        optical_depth += cloud_optical_depth[i][j];

      config.delta_tau[j] = optical_depth;
    }

    // Set wavenumber for Planck function
    const double wavenumber = spectral_grid->wavenumber_list[i];
    config.wavenumber_low = wavenumber;
    config.wavenumber_high = wavenumber;

    // Enable thermal emission if Planck function is non-negligible
    config.use_thermal_emission =
      disortpp::planckFunction2(wavenumber, wavenumber, config.temperature_bottom) > 1.e-35;

    // Solve and extract upward flux at TOA
    // With index_from_bottom=true, index 0 = BOA, last = TOA
    auto result = solvers[thread]->solve(config);
    spectrum[i] = result.flux_up.back() * spectrum_scaling;
  }
}


}
