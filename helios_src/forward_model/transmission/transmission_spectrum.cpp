/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2022 Daniel Kitzmann
*
* Helios-r2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Helios-r2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* Helios-r2 directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#include "transmission.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>


namespace helios{


void TransmissionModel::calcTransmissionSpectrum(
  const double bottom_radius, 
  const double star_radius, 
  std::vector<double>& spectrum)
{

  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<spectrum.size(); ++i)
  {
    const double planet_radius = calcTransitRadius(i, bottom_radius);

    spectrum[i] = planet_radius*planet_radius/star_radius/star_radius *1e6;
  }

}



double TransmissionModel::calcTransitRadius(const unsigned int wavelength, const double bottom_radius)
{
  const double effective_tangent_height = integrateEffectiveTangentHeight(wavelength, bottom_radius);
  
  return effective_tangent_height + bottom_radius;
}



double TransmissionModel::integrateEffectiveTangentHeight(
  const unsigned int wavelength, 
  const double bottom_radius)
{
  std::vector<double> path_transmission = tangentPathsTransmission(wavelength, bottom_radius);

  double effective_tangent_height = 0;

  for (size_t i=0; i<nb_grid_points-1; i++)
  {
    effective_tangent_height += 2. * ( (bottom_radius + atmosphere.altitude[i]) * (1. - path_transmission[i])
                                     + (bottom_radius + atmosphere.altitude[i+1]) * (1. - path_transmission[i+1]) )
                                   * (atmosphere.altitude[i+1] - atmosphere.altitude[i]) * 0.5;
  }
    
  effective_tangent_height = std::sqrt(effective_tangent_height + bottom_radius*bottom_radius) - bottom_radius;

  return effective_tangent_height;
}




std::vector<double> TransmissionModel::tangentPathsTransmission(
  const unsigned int wavelength, 
  const double bottom_radius)
{
  std::vector<double> path_transmission(nb_grid_points, 0);
  
  //the top limb is transparent (with only a single point we can't perform an integration)
  path_transmission.back() = 1;
  
  //we start the calculation at the radius just below the top and proceed downwards
  for (int i=nb_grid_points-2; i>-1; --i)
  {
    const double optical_depth = tangentOpticalDepth(i, wavelength, bottom_radius);
    
    path_transmission[i] = std::exp(-optical_depth);
    
    //if the transmission becomes too small, we skip the rest
    //if (path_transmission[i] < transmission_cutoff) break;
  }

  return path_transmission;
}




double TransmissionModel::tangentOpticalDepth(
  const unsigned int tangent_radius, 
  const unsigned int wavelength, 
  const double bottom_radius)
{
  double tangent_optical_depth = 0;
  
  //trapezoidal rule to calculate the optical depth along the path
  //since the problem is symmetric around the tangent point, 
  //we multiply the optical depth by 2 to account for both hemispheres
  for (size_t i=tangent_radius; i<nb_grid_points-1; ++i)
  {
    const double path_length = distanceToTangentCenter(tangent_radius, i+1, bottom_radius) 
                            -  distanceToTangentCenter(tangent_radius, i, bottom_radius);

    const double extinction_coeff1 = opacity_calc.absorption_coeff[wavelength][i] + opacity_calc.scattering_coeff[wavelength][i];
    const double extinction_coeff2 = opacity_calc.absorption_coeff[wavelength][i+1] + opacity_calc.scattering_coeff[wavelength][i+1];
    
    //account for both hemispheres by multiplying the results by 2
    //note: this cancels the factor of 0.5 from the trapezoidal rule
    tangent_optical_depth += path_length * (extinction_coeff1 + extinction_coeff2);

    if (tangent_optical_depth > transmission_optical_depth_cutoff)
      break;
  }

  return tangent_optical_depth;
}



double TransmissionModel::distanceToTangentCenter(
  const unsigned int tangent_altitude, 
  const unsigned int altitude, 
  const double bottom_radius)
{
  const double a = bottom_radius + atmosphere.altitude[altitude];
  const double b = bottom_radius + atmosphere.altitude[tangent_altitude];

  return std::sqrt(a*a - b*b);
}


}