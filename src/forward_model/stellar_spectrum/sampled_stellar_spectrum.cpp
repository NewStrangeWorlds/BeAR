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
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>
#include <algorithm> 
#include <assert.h>

#include "stellar_spectrum_grid.h"

#include "../../spectral_grid/spectral_grid.h"
#include "../../additional/exceptions.h"


namespace bear{


void SampledStellarSpectrum::deleteSampledData()
{
  std::vector<double>().swap (spectrum);

  is_sampled = false;
}


void SampledStellarSpectrum::sampleSpectrum(
  SpectralGrid* spectral_grid,
  const std::vector<double>& grid_wavelengths,
  const bool use_gpu)
{
  if (is_sampled == true) return;

  if (!data_file.is_loaded) data_file.loadFile();

  if (data_file.spectrum.size() != grid_wavelengths.size()) 
    std::cout << "spectrum file " << data_file.file_path << " has an incorrect number of spectral points!\n";

  assert (data_file.spectrum.size() != grid_wavelengths.size());

  //convert from W m-2 mu-1 to W m-2 cm
  for (size_t i=0; i<data_file.spectrum.size(); ++i)
    data_file.spectrum[i] = data_file.spectrum[i]*grid_wavelengths[i]*grid_wavelengths[i]/10000.;


  spectrum = spectral_grid->interpolateToWavelengthGrid(grid_wavelengths, data_file.spectrum, false);


  //remove zeros from the spectrum
  for (auto & s : spectrum)
    if (s < 1e-40) s = 1e-40;

  if (use_gpu)
    moveToDevice(spectrum_gpu, spectrum, true);

  is_sampled = true;

  data_file.unloadData();
}



SampledStellarSpectrum::~SampledStellarSpectrum()
{
  deleteFromDevice(spectrum_gpu);
}

} 
