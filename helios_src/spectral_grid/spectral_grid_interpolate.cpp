/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
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


#include "spectral_band.h"


#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>


#include "spectral_grid.h"
#include "../additional/quadrature.h"
#include "../additional/aux_functions.h"


namespace helios{



std::vector<double> SpectralGrid::interpolateToWavenumberGrid(const std::vector<double>& data_x, const std::vector<double>& data_y, const bool log_interpolation)
{
  std::vector<double> x = data_x;
  std::vector<double> y = data_y;

  if (x[0] > x[1])
  {
    std::reverse(x.begin(), x.end());
    std::reverse(y.begin(), y.end());
  }

  x = convertWavenumbersToWavelengths(x);

  return interpolateToWavelengthGrid(x, y, log_interpolation);
}



std::vector<double> SpectralGrid::interpolateToWavelengthGrid(const std::vector<double>& data_x, const std::vector<double>& data_y, const bool log_interpolation)
{
  std::vector<double> x = data_x;
  std::vector<double> y = data_y;

  if (x[0] < x[1])
  {
    std::reverse(x.begin(), x.end());
    std::reverse(y.begin(), y.end());
  }

  if (log_interpolation == true)
    for (auto & i : y) i = std::log10(i);

  std::vector<double> interpolated_data(nb_spectral_points, 0.0);


  auto linearInterpolation = [] (const double x1, const double x2, const double y1, const double y2, const double x){return y1 + (y2 - y1) * (x - x1)/(x2 - x1);};


  size_t x_start = 0;

  for (size_t i=0; i<nb_spectral_points; ++i)
  {
    if (wavelength_list[i] > x.front()) continue;
    if (wavelength_list[i] < x.back()) break;

    auto it = std::find_if(x.cbegin()+x_start, x.cend(), [this, i](double val){return (val <= wavelength_list[i]); } );

    std::size_t index = std::distance(x.cbegin(), it);

    if (*it == wavelength_list[i])
      interpolated_data[i] = y[index];
    else
      interpolated_data[i] = linearInterpolation(*(it-1), *it, y[index-1], y[index], wavelength_list[i]);
    
    if (x_start > 0)
      x_start = index-1;
  }


  if (log_interpolation == true)
    for (auto & i : interpolated_data) if (i != 0.0) i = std::pow(10, i);


  return interpolated_data;
}






}
