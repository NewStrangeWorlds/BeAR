// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _discrete_ordinate_adding_h
#define _discrete_ordinate_adding_h


#include <vector>
#include <cstddef>

#include "../radiative_transfer.h"
#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../spectral_grid/spectral_grid.h"


namespace bear {


//Thermal emission with multiple scattering using the four-stream discrete ordinate
//adding method of Zhang et al. (2013, JQSRT 122, 64), as implemented in Exo-FMS (lw_DDA).
//Each layer is solved analytically for its reflection and transmission matrices and
//thermal sources, the layers are then combined with the adding method from the bottom up.
//Uses double-Gauss quadrature (two Gauss-Legendre angles per hemisphere), the Rayleigh
//phase function for the gas, a Henyey-Greenstein phase function for the clouds and
//optional delta-M scaling. Without scattering the result is identical to the
//short characteristics solver.
class DiscreteOrdinateAdding : public RadiativeTransfer{
  public:
    DiscreteOrdinateAdding(
      SpectralGrid* spectral_grid_ptr,
      const size_t nb_grid_points,
      const bool use_gpu,
      const bool use_delta_m);
    virtual ~DiscreteOrdinateAdding() {}

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
    size_t nb_grid_points = 0;
    bool gpu_enabled = false;
    bool use_delta_m = true;

    double calcSpectrum(
      const Atmosphere& atmosphere,
      const std::vector<double>& absorption_coeff,
      const std::vector<double>& scattering_coeff,
      const std::vector<double>& cloud_optical_depth,
      const std::vector<double>& cloud_single_scattering,
      const std::vector<double>& cloud_asym_param,
      const size_t nu_index) const;
};


}
#endif
