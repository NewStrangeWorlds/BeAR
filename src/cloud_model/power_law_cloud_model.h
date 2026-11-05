// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _power_law_cloud_model_h
#define _power_law_cloud_model_h


#include <vector>

#include "cloud_model.h"
#include "grey_cloud_model.h"
#include "cloud_model.h"
#include "../forward_model/atmosphere/atmosphere.h"


namespace bear {


class PowerLawCloudModel: public GreyCloudModel{
  public:
    PowerLawCloudModel(const std::vector<std::string>& parameters);
    virtual ~PowerLawCloudModel() {}
    virtual void opticalProperties(
      const std::vector<double>& parameters, 
      const Atmosphere& atmosphere,
      SpectralGrid* spectral_grid,
      std::vector<std::vector<double>>& optical_depth, 
      std::vector<std::vector<double>>& single_scattering, 
      std::vector<std::vector<double>>& asym_param);
    virtual void opticalPropertiesGPU(
      const std::vector<double>& parameters, 
      const Atmosphere& atmosphere,
      SpectralGrid* spectral_grid,
      float* optical_depth_dev, 
      float* single_scattering_dev, 
      float* asym_param);

  protected:
    double reference_wavelength = 0;
};


}
#endif 

