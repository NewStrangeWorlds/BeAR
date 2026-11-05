// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _module_h
#define _module_h

#include <vector>
#include <string>
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
      float* spectrum_gpu) = 0;
    size_t nbParameters() {return parameter_names.size();}
    const std::vector<std::string>& parameterNames() const {return parameter_names;}
  protected:
    std::vector<std::string> parameter_names;
};

}

#endif