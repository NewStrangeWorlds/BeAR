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


#ifndef _velocity_broadening_h
#define _velocity_broadening_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <memory>

#include "../module.h"

#include "../../../spectral_grid/spectral_grid.h"
#include "../../../CUDA_kernels/data_management_kernels.h"

namespace bear {


class VelocityBroadening : public Module{
  public:
    VelocityBroadening (
      const std::vector<std::string>& velocity_broadening_parameters,
      SpectralGrid* spectral_grid_);
    virtual ~VelocityBroadening() {}
    
    virtual void modifySpectrum(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      std::vector<double>& spectrum);
    
    virtual void modifySpectrumGPU(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      float* spectrum_gpu);
  protected:
    SpectralGrid* spectral_grid;
};


}

#endif