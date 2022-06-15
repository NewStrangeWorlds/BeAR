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


#include "observations.h"


#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "../spectral_grid/spectral_band_type.h"
#include "../spectral_grid/spectral_band.h"
#include "../retrieval/retrieval.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"


namespace helios{



void Observation::init (Retrieval* retrieval_ptr, const std::string& file_name)
{
  retrieval = retrieval_ptr;


  loadFile(file_name);
}



Observation::~Observation()
{
  if (retrieval->config->use_gpu)
  {
    deleteFromDevice(filter_response_gpu);
    deleteFromDevice(filter_response_weight_gpu);
    deleteFromDevice(spectrum_filter_dev);
    deleteFromDevice(spectrum_convolved_dev);
  }
}



void Observation::initDeviceMemory()
{
  if (retrieval->config->use_gpu)
  {
    if (filter_response.size() != 0)
    allocateOnDevice(spectrum_filter_dev, retrieval->spectral_grid.nbSpectralPoints());

    if (instrument_profile_fwhm.size() != 0)
      allocateOnDevice(spectrum_convolved_dev, retrieval->spectral_grid.nbSpectralPoints());
  }
  
}



std::vector<double> Observation::processModelSpectrum(
  const std::vector<double> spectrum, 
  const bool is_flux)
{
  std::vector<double> spectrum_filter = applyFilterResponseFunction(spectrum);

  std::vector<double> convolved_spectrum = spectral_bands.convolveSpectrum(spectrum_filter);

  std::vector<double> spectrum_bands = spectral_bands.bandIntegrateSpectrum(convolved_spectrum, is_flux);

  return spectrum_bands;
}



void Observation::processModelSpectrumGPU(
  double* spectrum,
  double* spectrum_bands,
  const unsigned int start_index,
  const bool is_flux)
{
  if (filter_response.size() != 0)
    applyFilterResponseGPU(spectrum);
  else
    spectrum_filter_dev = spectrum;

  if (instrument_profile_fwhm.size() != 0)
    spectral_bands.convolveSpectrumGPU(
      spectrum_filter_dev, 
      spectrum_convolved_dev);
  else
    spectrum_convolved_dev = spectrum_filter_dev;

  spectral_bands.bandIntegrateSpectrumGPU(
    spectrum_convolved_dev, 
    spectrum_bands, 
    start_index, 
    is_flux);
}




}
