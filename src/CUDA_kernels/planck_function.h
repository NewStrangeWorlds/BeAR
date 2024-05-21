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


#ifndef _planck_function_kernel_h
#define _planck_function_kernel_h

#include "../additional/physical_const.h"

namespace helios{


__forceinline__ __device__ double planckFunction(const double temperature, const double wavenumber)
{

  return 2. * constants::planck_h * constants::light_c * constants::light_c * wavenumber*wavenumber*wavenumber 
         / ( exp(constants::planck_h * wavenumber * constants::light_c / constants::boltzmann_k / temperature) - 1.0);

}


}

#endif