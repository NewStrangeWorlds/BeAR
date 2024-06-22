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


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>

#include "star_blackbody.h"

#include "../../spectral_grid/spectral_grid.h"
#include "../../additional/physical_const.h"
#include "../../additional/aux_functions.h"
#include "../../additional/exceptions.h"


namespace bear{


std::vector<double> StarBlackBody::calcFlux(
  const std::vector<double>& parameter)
{
  std::vector<double> flux(spectral_grid->nbSpectralPoints(), 0);

  const double effective_temperature = parameter[0];

  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<flux.size(); ++i)
    flux[i] = aux::planckFunctionWavenumber(
      effective_temperature, 
      spectral_grid->wavenumber_list[i]) * constants::pi * 1e-3;

  return flux;
}



}