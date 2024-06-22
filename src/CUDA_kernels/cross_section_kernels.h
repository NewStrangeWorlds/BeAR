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


#ifndef _cross_section_kernels_h
#define _cross_section_kernels_h

#include <vector>

namespace bear{


void calcCrossSectionsHost(double* cross_sections1, double* cross_sections2, double* cross_sections3, double* cross_sections4,
                           const double temperature1, const double temperature2,
                           const double log_pressure1, const double log_pressure2,
                           const double temperature, const double log_pressure, const double number_density,
                           const size_t nb_spectral_points, const size_t nb_grid_points, const size_t grid_point,
                           double* absorption_coeff_device, double* scattering_coeff_device);


void calcCIACoefficientsHost(double* cross_sections1, double* cross_sections2,
                             const double temperature1, const double temperature2,
                             const double temperature, const double number_densities,
                             const size_t nb_spectral_points, const size_t nb_grid_points, const size_t grid_point,
                             double* absorption_coeff_device);


void calcHmContinuumHost(const double hm_number_density,
                         const double h_number_density,
                         const double e_pressure,
                         const double temperature,
                         const int nb_spectral_points, const int grid_point,
                         double* wavelengths_device,
                         double* absorption_coeff_device);


void initCrossSectionsHost(const size_t nb_points,  double* absorption_coeff_device);

void checkCrossSectionsHost(const size_t nb_spectral_points, const size_t nb_grid_points,
                                     double* absorption_coeff_device);

}

#endif
