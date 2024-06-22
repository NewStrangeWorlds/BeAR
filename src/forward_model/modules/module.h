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


#ifndef _module_h
#define _module_h

#include <vector>
#include <iostream>
#include "../../CUDA_kernels/data_management_kernels.h"

#include "../../spectral_grid/spectral_grid.h"
#include "../atmosphere/atmosphere.h"

namespace bear {


class Module{
  public:
    virtual ~Module() {}
    virtual void modifySpectrum(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      std::vector<double>& spectrum) = 0;
    virtual void modifySpectrumGPU(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      double* spectrum_gpu) = 0;
    size_t nbParameters() {return nb_parameters;}
  protected:
    size_t nb_parameters {};
};

}

#endif