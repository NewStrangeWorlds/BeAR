// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _fixed_cloud_model_h
#define _fixed_cloud_model_h


#include <vector>

#include "cloud_model.h"
#include "../forward_model/atmosphere/atmosphere.h"


namespace bear {


class FixedCloudModel: public CloudModel{
  public:
    FixedCloudModel(
      const std::vector<std::vector<double>>& optical_depth,
      const std::vector<std::vector<double>>& single_scattering_albedo,
      const std::vector<std::vector<double>>& asymmetry_parameter);
    virtual ~FixedCloudModel() {}
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
      std::vector<std::vector<double>> fixed_optical_depth;
      std::vector<std::vector<double>> fixed_single_scattering_albedo;
      std::vector<std::vector<double>> fixed_asymmetry_parameter;
};


}
#endif 

