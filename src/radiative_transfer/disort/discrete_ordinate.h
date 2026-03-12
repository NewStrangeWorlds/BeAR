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


#ifndef _discrete_ordinate_h
#define _discrete_ordinate_h


#include <vector>
#include <iostream>
#include <cmath>
#include <memory>

#include "../radiative_transfer.h"
#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../spectral_grid/spectral_grid.h"

#include <DisortFluxConfig.hpp>
#include <FluxResult.hpp>


namespace bear {


struct FluxSolverBase {
  virtual ~FluxSolverBase() = default;
  virtual disortpp::FluxResult solve(disortpp::DisortFluxConfig& config) = 0;
};

template<int NStr>
struct FluxSolverWrapper;


class DiscreteOrdinates : public RadiativeTransfer{
  public:
    DiscreteOrdinates(
      SpectralGrid* spectral_grid_ptr,
      const size_t nb_streams,
      const size_t nb_grid_points,
      const bool use_gpu);
    virtual ~DiscreteOrdinates() {}

    virtual void calcSpectrum(
      const Atmosphere& atmosphere,
      const std::vector< std::vector<double> >& absorption_coeff,
      const std::vector< std::vector<double> >& scattering_coeff,
      const std::vector< std::vector<double> >& cloud_optical_depth,
      const std::vector< std::vector<double> >& cloud_single_scattering,
      const std::vector< std::vector<double> >& cloud_asym_param,
      const double spectrum_scaling,
      std::vector<double>& spectrum);
    virtual void calcSpectrumGPU(
      const Atmosphere& atmosphere,
      float* absorption_coeff_dev,
      float* scattering_coeff_dev,
      float* cloud_optical_depth,
      float* cloud_single_scattering,
      float* cloud_asym_param,
      const double spectrum_scaling,
      float* model_spectrum_dev)
      {
        std::cout << "Sorry, DisORT has no GPU option :(\n";
      }
  private:
    size_t nb_streams;
    size_t nb_grid_points;

    std::vector<std::unique_ptr<FluxSolverBase>> solvers;
    std::vector<disortpp::DisortFluxConfig> configs;
};


}
#endif


