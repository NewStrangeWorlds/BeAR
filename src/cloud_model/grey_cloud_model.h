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
      double* optical_depth_dev, 
      double* single_scattering_dev, 
      double* asym_param);

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

