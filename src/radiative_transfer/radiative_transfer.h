// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _radiative_transfer_h
#define _radiative_transfer_h

#include <vector>

#include "../forward_model/atmosphere/atmosphere.h"
#include "../spectral_grid/spectral_grid.h"


namespace bear {


class RadiativeTransfer{
  public:
    RadiativeTransfer(SpectralGrid* spectral_grid_) {
      spectral_grid = spectral_grid_;};
    virtual ~RadiativeTransfer() {}
    virtual void calcSpectrum(
      const Atmosphere& atmosphere,
      const std::vector< std::vector<double> >& absorption_coeff, 
      const std::vector< std::vector<double> >& scattering_coeff,
      const std::vector< std::vector<double> >& cloud_optical_depth,
      const std::vector< std::vector<double> >& cloud_single_scattering,
      const std::vector< std::vector<double> >& cloud_asym_param,
      const double spectrum_scaling,
      std::vector<double>& spectrum) = 0;
    virtual void calcSpectrumGPU(
      const Atmosphere& atmosphere,
      float* absorption_coeff_dev,
      float* scattering_coeff_dev,
      float* cloud_optical_depth,
      float* cloud_single_scattering,
      float* cloud_asym_param,
      const double spectrum_scaling,
      float* model_spectrum_dev) = 0;
    
    //change units of high-res spectrum from cm to micron^-1
    void changeSpectrumUnits(std::vector<double>& spectrum) {
      for (size_t i=0; i<spectrum.size(); ++i)
        spectrum[i] = spectrum[i]/spectral_grid->wavelength_list[i]/spectral_grid->wavelength_list[i]*10000.0;};
    void changeSpectrumUnitsGPU(float* spectrum);
  
  protected:
    SpectralGrid* spectral_grid = nullptr;
};


}
#endif

