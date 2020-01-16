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


#ifndef _cross_section_kernels_h
#define _cross_section_kernels_h

#include <vector>

namespace helios{


void calcCrossSectionsHost(double* cross_sections1, double* cross_sections2, double* cross_sections3, double* cross_sections4,
                           const double temperature1, const double temperature2,
                           const double pressure1, const double pressure2,
                           const double temperature, const double pressure, const double number_density,
                           const size_t nb_spectral_points, const size_t nb_grid_points, const size_t grid_point,
                           double* absorption_coeff_device, double* scattering_coeff_device);


void calcCIACoefficientsHost(double* cross_sections1, double* cross_sections2,
                             const double temperature1, const double temperature2,
                             const double temperature, const double number_densities,
                             const size_t nb_spectral_points, const size_t nb_grid_points, const size_t grid_point,
                             double* absorption_coeff_device);


void initCrossSectionsHost(const size_t nb_points,  double* absorption_coeff_device);

void checkCrossSectionsHost(const size_t nb_spectral_points, const size_t nb_grid_points,
                                     double* absorption_coeff_device);

}

#endif
