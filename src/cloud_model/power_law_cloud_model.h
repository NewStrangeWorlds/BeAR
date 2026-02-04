/*
* This file is part of the BeAr code (https://github.com/newstrangeworlds/bear).
* Copyright (C) 2024 Daniel Kitzmann
*
* BeAr is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* BeAr is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* BeAr directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


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

