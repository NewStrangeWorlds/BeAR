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

#include "stellar_contamination.h"

#include "../../spectral_grid/spectral_grid.h"
#include "../../additional/physical_const.h"
#include "../../additional/aux_functions.h"
#include "../../additional/exceptions.h"
#include "../stellar_spectrum/select_stellar_model.h"


namespace bear{


StellarContamination::StellarContamination (
  const std::vector<std::string>& stellar_model_parameters,
  SpectralGrid* spectral_grid_)
  : spectral_grid(spectral_grid_)
{
  std::vector<std::string> parameters(stellar_model_parameters.begin()+1, stellar_model_parameters.end());

  stellar_model = selectStellarModel(
    stellar_model_parameters[0],
    parameters,
    spectral_grid);
  
  nb_stellar_model_param = stellar_model->nbParameters();
  nb_parameters = nb_stellar_model_param + 4;
}


void StellarContamination::modifySpectrum(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      std::vector<double>& spectrum)
{
  std::vector<double> stellar_param(
    parameter.begin(), 
    parameter.begin()+nb_stellar_model_param);
  
  const double temperature_phot = parameter[0];

  const double temperature_fac = temperature_phot + parameter[nb_stellar_model_param];
  const double temperature_spot = temperature_phot - parameter[nb_stellar_model_param+1];
  const double fraction_fac = parameter[nb_stellar_model_param+2];
  const double fraction_spot = parameter[nb_stellar_model_param+3];

  const std::vector<double> spectrum_phot = stellar_model->calcFlux(stellar_param);


  std::vector<double> spectrum_fac(spectrum.size(), 0.0);

  if (fraction_fac > 0)
  {
    stellar_param[0] = temperature_fac;
    spectrum_fac = stellar_model->calcFlux(stellar_param);
  }

  std::vector<double> spectrum_spot(spectrum.size(), 0.0);

  if (fraction_spot > 0)
  {
    stellar_param[0] = temperature_spot;
    spectrum_spot = stellar_model->calcFlux(stellar_param);
  }

  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<spectrum.size(); ++i)
  {
    const double stellar_activity = 1.0 
      - fraction_spot * (1.0 - spectrum_spot[i]/spectrum_phot[i])
      - fraction_fac * (1.0 - spectrum_fac[i]/spectrum_phot[i]);

    spectrum[i] /= stellar_activity;
  }

}



}