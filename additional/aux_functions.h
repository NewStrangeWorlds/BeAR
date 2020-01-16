/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2019 Daniel Kitzmann
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


#ifndef _aux_functions_h
#define _aux_functions_h

#include <vector>


namespace helios{ namespace aux{


double planckFunctionWavenumber(const double temperature, const double wavenumber);
double linearInterpolation(const double x1, const double x2, const double y1, const double y2, const double x);
double voigtProfile(const double x, const double gaussian_width, const double lorentz_width);

double normalDistribution(const double mu, const double sigma, const double x);
double normalDistribution(const double sigma, const double x);

std::vector<double> interpolateToWavenumberGrid(const std::vector<double>& wavenumber_data, const std::vector<double>& data,
                                                const std::vector<double>& wavenumber_grid,
                                                const double outside_range_value,
                                                const bool interpolate_log);

}}



#endif
