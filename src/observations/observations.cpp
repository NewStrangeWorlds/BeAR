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
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"


namespace helios{


void Observation::init(const std::string& file_name, const std::string spectrum_modifier_id)
{
  setSpectrumModifier(spectrum_modifier_id);
  loadFile(file_name);
}



Observation::~Observation()
{
  if (config->use_gpu)
  {
    deleteFromDevice(filter_response_gpu);
    deleteFromDevice(filter_response_weight_gpu);
    deleteFromDevice(spectrum_filter_dev);
    deleteFromDevice(spectrum_convolved_dev);
  }
}



void Observation::initDeviceMemory()
{
  if (config->use_gpu)
  {
    if (filter_response.size() != 0)
    allocateOnDevice(spectrum_filter_dev, spectral_grid->nbSpectralPoints());

    if (instrument_profile_fwhm.size() != 0)
      allocateOnDevice(spectrum_convolved_dev, spectral_grid->nbSpectralPoints());
  }
  
}



std::vector<double> Observation::processModelSpectrum(
  const std::vector<double> spectrum, 
  const bool is_flux)
{
  std::vector<double> spectrum_filter = applyFilterResponseFunction(spectrum);

  std::vector<double> convolved_spectrum = spectral_bands.convolveSpectrum(spectrum_filter);

  std::vector<double> spectrum_bands = spectral_bands.bandIntegrateSpectrum(
    convolved_spectrum, 
    is_flux,
    filter_response.size() != 0);

  return spectrum_bands;
}



void Observation::processModelSpectrumGPU(
  double* spectrum,
  double* spectrum_bands,
  const unsigned int start_index,
  const bool is_flux)
{ 
  bool use_filter_response = filter_response.size() != 0;

  if (use_filter_response)
  {
    applyFilterResponseGPU(spectrum);

    if (instrument_profile_fwhm.size() != 0)
    {
      spectral_bands.convolveSpectrumGPU(
        spectrum_filter_dev,
        spectrum_convolved_dev);

      spectral_bands.bandIntegrateSpectrumGPU(
        spectrum_convolved_dev,
        spectrum_bands,
        start_index,
        is_flux,
        use_filter_response);
    }
    else
    {
      spectral_bands.bandIntegrateSpectrumGPU(
        spectrum_filter_dev,
        spectrum_bands,
        start_index,
        is_flux,
        use_filter_response);
    }
  }
  else if (instrument_profile_fwhm.size() != 0)
  {
    spectral_bands.convolveSpectrumGPU(
      spectrum,
      spectrum_convolved_dev);

    spectral_bands.bandIntegrateSpectrumGPU(
      spectrum_convolved_dev,
      spectrum_bands,
      start_index,
      is_flux,
      use_filter_response);
  }
  else
    spectral_bands.bandIntegrateSpectrumGPU(
      spectrum,
      spectrum_bands,
      start_index,
      is_flux,
      use_filter_response);
}



void Observation::addShiftToSpectrum(
  std::vector<double>& spectrum_bands,
  const double spectrum_shift)
{

  for (auto & s : spectrum_bands)
    s += spectrum_shift;

}



void Observation::setSpectrumModifier(const std::string modifier_id)
{
  if (modifier_id == "none" or modifier_id == "")
  {
    spectrum_modifier = observation_modifiers::id::none;
    nb_modifier_param = 0;
  }
  else if (modifier_id == "shift_const")
  {
    spectrum_modifier = observation_modifiers::id::shift_const;
    nb_modifier_param = 1;
  }
  else
  {
    std::string error_message = "Observation modifier " + modifier_id + " in observations.list unknown.\n";
    throw InvalidInput(std::string ("Observation::setSpectrumModifier"), error_message);
  }

}




}
