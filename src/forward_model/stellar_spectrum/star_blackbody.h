// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _star_blackbody_h
#define _star_blackbody_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "stellar_spectrum.h"
#include "../../spectral_grid/spectral_grid.h"


namespace bear {


class StarBlackBody : public StellarSpectrumModel{
  public:
    StarBlackBody (SpectralGrid* spectral_grid_) 
      : spectral_grid(spectral_grid_)
      {parameter_names = {"star_t_eff"};}
    virtual ~StarBlackBody() {}
    
    virtual std::vector<double> calcFlux(
      const std::vector<double>& parameter);
    
    virtual void calcFluxGPU(
      const std::vector<double>& parameter,
      float* spectrum_gpu);
  protected:
    SpectralGrid* spectral_grid;
};


}

#endif