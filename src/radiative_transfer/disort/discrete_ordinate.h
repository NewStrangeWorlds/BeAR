// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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


