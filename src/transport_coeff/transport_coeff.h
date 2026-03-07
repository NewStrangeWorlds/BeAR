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

    void calculateGPU(
      const double temperature,
      const double pressure,
      const std::vector<double>& number_densities,
      const size_t nb_grid_points,
      const size_t grid_point,
      float* absorption_coeff_device,
      float* scattering_coeff_device);

    void prepareBatchedGPU(
      const Atmosphere& atmosphere,
      std::vector<float*>& cs1_ptrs, std::vector<float*>& cs2_ptrs,
      std::vector<float*>& cs3_ptrs, std::vector<float*>& cs4_ptrs,
      std::vector<float>& temp_factors, std::vector<float>& pres_factors,
      std::vector<float>& cs_log_number_densities, std::vector<int>& cs_grid_points,
      std::vector<float*>& ray_ptrs,
      std::vector<double>& ray_number_densities, std::vector<int>& ray_grid_points);

    void calculateContinuumGPU(
      const double temperature,
      const double pressure,
      const std::vector<double>& number_densities,
      const size_t nb_grid_points,
      const size_t grid_point,
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
