// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _adding_doubling_h
#define _adding_doubling_h


#include <vector>
#include <iostream>
#include <cmath>
#include <memory>

#include "../radiative_transfer.h"
#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../spectral_grid/spectral_grid.h"

#include <adding_doubling.hpp>


namespace bear {


class AddingDoubling : public RadiativeTransfer{
  public:
    AddingDoubling(
      SpectralGrid* spectral_grid_ptr,
      const size_t nb_quadrature,
      const size_t nb_grid_points,
      const bool use_gpu,
      const bool use_delta_m = true);
    virtual ~AddingDoubling();

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
      float* model_spectrum_dev);

  private:
    size_t nb_quadrature;
    size_t nb_grid_points;
    bool gpu_enabled = false;
    bool use_delta_m = true;

    std::vector<adrt::ADConfig> configs;
    std::vector<adrt::SolverWorkspace> workspaces;

    // CUDA batched solver data
    float* phase_moments_dev = nullptr;
    int nb_phase_moments = 0;
    void* gpu_workspace_ptr = nullptr;
};


}
#endif


