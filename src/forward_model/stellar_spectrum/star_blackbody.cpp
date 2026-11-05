// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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