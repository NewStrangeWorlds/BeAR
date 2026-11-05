// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef TRANSPORT_COEFF_H
#define TRANSPORT_COEFF_H

#include <vector>
#include <memory>

#include "opacity_species.h"
#include "../spectral_grid/spectral_grid.h"
#include "../config/global_config.h"
#include "../forward_model/atmosphere/atmosphere.h"


namespace bear{


class TransportCoefficients {
  public:
    TransportCoefficients(
      GlobalConfig* config_ptr,
      SpectralGrid* grid_ptr, 
      const std::vector<std::string>& opacity_species_symbol,
      const std::vector<std::string>& opacity_species_folder);
    ~TransportCoefficients();

    void calculate(
      const double temperature,
      const double pressure,
      const std::vector<double>& number_densities,
      std::vector<double>& absorption_coeff,
      std::vector<double>& scattering_coeff);

    void prepareBatchedGPU(
      const Atmosphere& atmosphere,
      std::vector<float*>& cs1_ptrs, std::vector<float*>& cs2_ptrs,
      std::vector<float*>& cs3_ptrs, std::vector<float*>& cs4_ptrs,
      std::vector<float>& temp_factors, std::vector<float>& pres_factors,
      std::vector<float>& cs_log_number_densities, std::vector<int>& cs_grid_points,
      std::vector<float*>& ray_ptrs,
      std::vector<double>& ray_number_densities, std::vector<int>& ray_grid_points);

    void calculateContinuumGPU(
      const std::vector<double>& temperature,
      const std::vector<std::vector<double>>& number_densities,
      float* absorption_coeff_device);

  private:
    GlobalConfig* config = nullptr;
    SpectralGrid* spectral_grid = nullptr;

    std::vector<std::unique_ptr<OpacitySpecies>> gas_species;

    bool addOpacitySpecies(
      const std::string& species_symbol, const std::string& species_folder);
};


}

#endif
