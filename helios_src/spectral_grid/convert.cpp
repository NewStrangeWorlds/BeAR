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


#include "spectral_grid.h"


#include <iostream>
#include <string>
#include <fstream>
#include <vector>


#include "../config/global_config.h"



namespace helios{



//convert wavenumbers in cm-1 to wavelengths in microns
void SpectralGrid::convertWavenumbersToWavelengths(const std::vector<double>& wavenumbers, std::vector<double>& wavelengths)
{
  size_t nb_wavenumbers = wavenumbers.size();

  if (nb_wavenumbers == 0) return;


  wavelengths.assign(nb_wavenumbers, 0);


  for (size_t i=0; i<nb_wavenumbers; ++i)
    wavelengths[i] = 1.0/wavenumbers[i] * 1e4;
}


//convert wavenumbers in cm-1 to wavelengths in microns
std::vector<double> SpectralGrid::convertWavenumbersToWavelengths(const std::vector<double>& wavenumbers)
{
  size_t nb_wavenumbers = wavenumbers.size();

  std::vector<double> wavelengths(nb_wavenumbers, 0);

  if (nb_wavenumbers == 0) return wavelengths;


  for (size_t i=0; i<nb_wavenumbers; ++i)
    wavelengths[i] = 1.0/wavenumbers[i] * 1e4;

  return wavelengths;
}



//convert wavelengths in microns to wavenumbers in cm-1
void SpectralGrid::convertWavelengthsToWavenumbers(const std::vector<double>& wavelengths, std::vector<double>& wavenumbers)
{
  size_t nb_wavelengths = wavelengths.size();

  if (nb_wavelengths == 0) return;


  wavenumbers.assign(nb_wavelengths, 0);


  for (size_t i=0; i<nb_wavelengths; ++i)
    wavenumbers[i] = 1.0/wavelengths[i] * 1e4;
}



std::vector<double> SpectralGrid::convertWavelengthsToWavenumbers(const std::vector<double>& wavelengths)
{
  size_t nb_wavelengths = wavelengths.size();

  std::vector<double> wavenumbers(nb_wavelengths, 0);

  if (nb_wavelengths == 0) return wavenumbers;


  for (size_t i=0; i<nb_wavelengths; ++i)
    wavenumbers[i] = 1.0/wavelengths[i] * 1e4;

  return wavenumbers;
}




}

