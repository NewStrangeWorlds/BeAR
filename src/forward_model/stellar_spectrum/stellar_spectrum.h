// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _stellar_spectrum_h
#define _stellar_spectrum_h

#include <vector>
#include <string>
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
      float* spectrum_gpu) = 0;
    size_t nbParameters() {return parameter_names.size();}
    const std::vector<std::string>& parameterNames() const {return parameter_names;}
  protected:
    std::vector<std::string> parameter_names;
};

}

#endif