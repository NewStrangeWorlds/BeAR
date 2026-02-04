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


#ifndef _cloud_model_h
#define _cloud_model_h

#include <vector>

#include "../forward_model/atmosphere/atmosphere.h"
#include "../spectral_grid/spectral_grid.h"


namespace bear {


class CloudModel{
  public:
    virtual ~CloudModel() {}
    virtual void opticalProperties(
      const std::vector<double>& parameters, 
      const Atmosphere& atmosphere,
      SpectralGrid* spectral_grid,
      std::vector<std::vector<double>>& optical_depth, 
      std::vector<std::vector<double>>& single_scattering, 
      std::vector<std::vector<double>>& asym_param) = 0;
    virtual void opticalPropertiesGPU(
      const std::vector<double>& parameters, 
      const Atmosphere& atmosphere,
      SpectralGrid* spectral_grid,
      float* optical_depth_dev, 
      float* single_scattering_dev, 
      float* asym_param) = 0;
    void convertOpticalDepth(
      std::vector<std::vector<double>>& optical_depth,
      std::vector<std::vector<double>>& extinction_coeff,
      std::vector<double>& altitude);
    void convertOpticalDepthGPU(
      float* optical_depth_dev,
      double* altitude,
      const size_t nb_grid_points,
      const size_t nb_spectral_points,
      float* extinction_coeff_dev);
    size_t nbParameters() {return nb_parameters;}
  protected:
    size_t nb_parameters {};
};


//converts the layer optical depth to a level-based extinction coefficient
inline void CloudModel::convertOpticalDepth(
  std::vector<std::vector<double>>& optical_depth,
  std::vector<std::vector<double>>& extinction_coeff,
  std::vector<double>& altitude)
{ 
  size_t nb_spectral_points = optical_depth.size();
  size_t nb_grid_points = optical_depth[0].size() + 1;

  extinction_coeff.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));

  for (size_t i=0; i<nb_spectral_points; ++i)
    for (size_t j=1; j<nb_grid_points; ++j)
    {
      double delta_z = altitude[j] - altitude[j-1];

      //convert optical depth to extinction coefficient
      //uses a finite difference to approximate the derivative
      extinction_coeff[i][j] = optical_depth[i][j-1] / delta_z;
    }


}


}
#endif 

