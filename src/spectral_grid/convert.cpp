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


#include "spectral_grid.h"


#include <iostream>
#include <string>
#include <fstream>
#include <vector>


#include "../config/global_config.h"



namespace helios{


//convert wavenumbers in cm-1 to wavelengths in microns
std::vector<double> SpectralGrid::wavenumberToWavelength(
  const std::vector<double>& wavenumbers)
{
  const size_t nb_wavenumbers = wavenumbers.size();

  if (nb_wavenumbers == 0) return std::vector<double>(0,0);

  std::vector<double> wavelengths(nb_wavenumbers, 0);

  for (size_t i=0; i<nb_wavenumbers; ++i)
    wavelengths[i] = wavenumberToWavelength(wavenumbers[i]);

  return wavelengths;
}



std::vector<double> SpectralGrid::wavelengthToWavenumber(
  const std::vector<double>& wavelengths)
{
  size_t nb_wavelengths = wavelengths.size();

   if (nb_wavelengths == 0) return std::vector<double>(0,0);

  std::vector<double> wavenumbers(nb_wavelengths, 0);

  for (size_t i=0; i<nb_wavelengths; ++i)
    wavenumbers[i] = wavelengthToWavenumber(wavelengths[i]);

  return wavenumbers;
}




}

