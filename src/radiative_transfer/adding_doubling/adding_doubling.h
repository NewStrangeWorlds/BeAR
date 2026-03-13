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
#include <cuda_solver.cuh>


namespace bear {


class AddingDoubling : public RadiativeTransfer{
  public:
    AddingDoubling(
      SpectralGrid* spectral_grid_ptr,
      const size_t nb_quadrature,
      const size_t nb_grid_points,
      const bool use_gpu);
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

    std::vector<adrt::ADConfig> configs;
    std::vector<adrt::SolverWorkspace> workspaces;

    // CUDA batched solver data
    float* phase_moments_dev = nullptr;
    int nb_phase_moments = 0;
    adrt::cuda::SolverWorkspaceGPU gpu_workspace;
};


}
#endif


