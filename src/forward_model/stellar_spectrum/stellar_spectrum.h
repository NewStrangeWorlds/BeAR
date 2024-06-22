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


#ifndef _stellar_spectrum_h
#define _stellar_spectrum_h

#include <vector>
#include <iostream>
#include "../../CUDA_kernels/data_management_kernels.h"

namespace bear {


class StellarSpectrumModel{
  public:
    virtual ~StellarSpectrumModel() {}
    virtual std::vector<double> calcFlux(
      const std::vector<double>& parameter) = 0;
    virtual void calcFluxGPU(
      const std::vector<double>& parameter,
      double* spectrum_gpu) = 0;
    size_t nbParameters() {return nb_parameters;}
  protected:
    size_t nb_parameters {};
};

}

#endif