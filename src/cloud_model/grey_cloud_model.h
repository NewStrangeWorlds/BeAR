// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _grey_cloud_model_h
#define _grey_cloud_model_h


#include <vector>

#include "cloud_model.h"
#include "../forward_model/atmosphere/atmosphere.h"


namespace bear {


class GreyCloudModel: public CloudModel{
  public:
    GreyCloudModel(const std::vector<std::string>& parameters);
    virtual ~GreyCloudModel() {}
    virtual void opticalProperties(
      const std::vector<double>& parameters, 
      const Atmosphere& atmosphere,
      SpectralGrid* spectral_grid,
      std::vector<std::vector<double>>& optical_depth, 
      std::vector<std::vector<double>>& single_scattering, 
      std::vector<std::vector<double>>& asym_param);
    virtual void opticalPropertiesGPU(
      const std::vector<double>& parameters, const Atmosphere& atmosphere,
      SpectralGrid* spectral_grid,
      float* optical_depth_dev, 
      float* single_scattering_dev, 
      float* asym_param);

  protected:
    bool fixed_bottom = false;

    void cloudPosition(
      const Atmosphere& atmosphere, 
      const double top_pressure, 
      const double bottom_pressure, 
      unsigned int& top_index, 
      unsigned int& bottom_index);
    void cloudPosition(
      const Atmosphere& atmosphere, 
      const double top_pressure, 
      unsigned int& top_index, 
      unsigned int& bottom_index);
};


}
#endif 

