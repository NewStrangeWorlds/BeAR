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


#include "brown_dwarf.h"


#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>


#include "../../observations/observations.h"
#include "../../additional/physical_const.h"
#include "../../retrieval/retrieval.h"



namespace helios{


//Calculate the effective temperatures from the spectra
//it uses the "observation" at the back of the observation vector
//assuming that this is the postprocess_spectrum_data
//this could in principle be a multi-band "observation" here, the method will sum up all fluxes
double BrownDwarfModel::postProcessEffectiveTemperature(const std::vector<double>& model_spectrum_bands)
{
  //the required "observation" is the last one
  size_t nb_bands = retrieval->observations.back().spectral_bands.nbBands();

  //the start index for this observation in the vector of band-integrated fluxes
  size_t band_start_index = 0;

  for (size_t i=0; i< retrieval->observations.size()-1; ++i)
    band_start_index += retrieval->observations[i].spectral_bands.nbBands();



  std::vector<double> fluxes(nb_bands, 0);
  double total_flux = 0;


  //sum up the total flux
  for (size_t i=0; i<nb_bands; ++i)
  {
    const double wavelength_edge_left =  1.0/retrieval->observations.back().spectral_bands.band_wavenumbers[i].front() * 1e4;
    const double wavelength_edge_right =  1.0/retrieval->observations.back().spectral_bands.band_wavenumbers[i].back() * 1e4;


    //convert back from flux per wavelength to total flux
    fluxes[i] = model_spectrum_bands[band_start_index + i] * (wavelength_edge_left - wavelength_edge_right);
      
    //"undo" the radius-distance scaling from the forward model
    //we need the flux at the top of the planet's atmosphere
    total_flux += fluxes[i]/(radius_distance_scaling);
  }


  
  
  //use the Stefan-Boltzmann law to convert the total flux to effective temperatures
  //total_flux has units of W m-2 and needs to be converted to cgs
  double effective_temperature = std::pow(total_flux*1000.0/constants::stefan_boltzmann, 0.25);

  
  return effective_temperature;
}


}
