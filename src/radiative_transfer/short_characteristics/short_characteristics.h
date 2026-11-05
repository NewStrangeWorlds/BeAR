// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _short_char_h
#define _short_char_h


#include <vector>
#include <iostream>
#include <cmath>

#include "../radiative_transfer.h"

#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../spectral_grid/spectral_grid.h"


namespace bear {


class ShortCharacteristics : public RadiativeTransfer{
  public:
    ShortCharacteristics(SpectralGrid* spectral_grid_ptr)
    : RadiativeTransfer(spectral_grid_ptr) {};
    virtual ~ShortCharacteristics() {}
    
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
    const std::vector<double> gauss_nodes{0.211324865405187, 0.788675134594813};
    const std::vector<double> gauss_weights{0.5, 0.5};
    const size_t nb_angles = gauss_nodes.size();

    double calcSpectrum(
      const std::vector<double>& absorption_coeff,
      const std::vector<double>& cloud_optical_depth,
      const std::vector<double>& temperature,
      const std::vector<double>& vertical_grid,
      const size_t nu_index);
};


}
#endif


