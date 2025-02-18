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


namespace bear{



void Observation::init(const ObservationInput& input)
{
  setSpectrumModifier(input.spectrum_modifier_id);

  observation_name = input.name;

  band_type::id band_type = selectBandType(input.type, observation_name);

  instrument_profile_fwhm = input.instrument_profile_fwhm;
  data = input.data;
  data_error = input.data_error;

  if (input.likelihood_weight.size() == 0)
    likelihood_weight = std::vector<double>(data.size(), 1.0);
  else
    likelihood_weight = input.likelihood_weight;

  filter_response_file = input.filter_response;
  
  if (filter_response_file.size() > 0)
  {
    filter_detector_type = input.filter_detector_type;

    if (input.filter_detector_type == "Energy")
      filter_detector_type = "energy";

    if (filter_detector_type == "Photon")
      filter_detector_type = "photon";

    if (filter_detector_type != "energy" && filter_detector_type != "photon")
    {
      std::string error_message = "Unsupported detector type *" + filter_detector_type + "\n";
      throw InvalidInput(std::string ("Observation::init"), error_message);
    }

    if (filter_response_file[0].front() < filter_response_file[0].back())
    {
      std::reverse(filter_response_file[0].begin(), filter_response_file[0].end());
      std::reverse(filter_response_file[1].begin(), filter_response_file[1].end());
    }
  }
  
  std::vector<std::vector<double>> bin_edges;
  std::vector<double> bin_centers;

  if (band_type == band_type::id::spectroscopy)
  {
    bin_centers = input.wavelengths;
    
    ascending_wavelengths = areWavelengthsAscending(bin_centers);

    bin_edges = calcBinEdges(bin_centers);
  }
  else
  { 
    bin_edges = input.bin_wavelength_edges;
    
    ascending_wavelengths = areWavelengthsAscending(bin_edges);

    bin_centers = calcBinCenters(bin_edges);
  }
 
  setObservationEdges(bin_edges);

  spectral_bands.init(
    wavelength_edges, 
    bin_edges, 
    bin_centers, 
    band_type);
}



void Observation::init(
  const std::string& file_name, 
  const std::string spectrum_modifier_id)
{
  setSpectrumModifier(spectrum_modifier_id);
  loadFile(file_name);
}



//wavelengths should be ordered in descending order because the 
//opacity data is organized in ascending wavenumbers
bool Observation::areWavelengthsAscending(
  std::vector<double>& wavelengths)
{
  bool ascending = true;

  if (wavelengths.size() > 1)
  {
    if (wavelengths[0] < wavelengths[1])
    {
      ascending = true;

      std::reverse(wavelengths.begin(), wavelengths.end());
      std::reverse(data.begin(), data.end());
      std::reverse(data_error.begin(), data_error.end());
      std::reverse(instrument_profile_fwhm.begin(), instrument_profile_fwhm.end());
      std::reverse(likelihood_weight.begin(), likelihood_weight.end());
    }
    else
      ascending = false;
  }
  else
    return false;

  return ascending;
}


//wavelengths should be ordered in descending order because the 
//opacity data is organized in ascending wavenumbers
bool Observation::areWavelengthsAscending(
  std::vector<std::vector<double>>& bin_edges)
{
  bool ascending = true;

  if (bin_edges.size() > 1)
  {
    if (bin_edges.front()[0] < bin_edges.back()[0])
    {
      ascending = true;

      std::reverse(bin_edges.begin(), bin_edges.end());
      std::reverse(data.begin(), data.end());
      std::reverse(data_error.begin(), data_error.end());
      std::reverse(instrument_profile_fwhm.begin(), instrument_profile_fwhm.end());
      std::reverse(likelihood_weight.begin(), likelihood_weight.end());
    }
    else
      ascending = false;
  }
  else
    ascending = false;

  for (size_t i=0; i<bin_edges.size(); ++i)
    if (bin_edges[i][0] < bin_edges[i][1]) 
      std::reverse(bin_edges[i].begin(), bin_edges[i].end());

  return ascending;
}



std::vector<std::vector<double>> Observation::calcBinEdges(
  const std::vector<double>& wavelengths)
{
  std::vector< std::vector<double> > bin_edges(
    wavelengths.size(), 
    std::vector<double>(2, 0));


  //set up the bin edges
  for (size_t i=0; i<wavelengths.size(); ++i)
  {
    
    if (i == 0)
    {
      //first bin
      bin_edges.front()[1] = wavelengths[0] - (wavelengths[i] - wavelengths[i+1]) * 0.5;
      bin_edges.front()[0] = wavelengths[0] + (wavelengths[i] - wavelengths[i+1]) * 0.5;
    }
    else if (i == wavelengths.size()-1)
    {
      //last bin
      bin_edges.back()[1] = wavelengths[i] - (wavelengths[i-1] - wavelengths[i]) * 0.5;
      bin_edges.back()[0] = wavelengths[i] + (wavelengths[i-1] - wavelengths[i]) * 0.5;
    }
    else
    {
      bin_edges[i][0] = wavelengths[i-1] - (wavelengths[i-1] - wavelengths[i]) * 0.5;
      bin_edges[i][1] = wavelengths[i] - (wavelengths[i] - wavelengths[i+1]) * 0.5;
    }

  }

  return bin_edges;
}


std::vector<double> Observation::calcBinCenters(
  const std::vector<std::vector<double>>& bin_edges)
{
  std::vector<double> bin_centers(bin_edges.size(), 0.0);

  for (size_t i=0; i<bin_centers.size(); ++i)
    bin_centers[i] = bin_edges[i][0] - (bin_edges[i][0] - bin_edges[i][1]) * 0.5; 
  
  return bin_centers;
}



void Observation::setObservationEdges(
  const std::vector<std::vector<double>>& bin_edges)
{
  //if we have an instrument profile, we extend the wavelength edges by 5 sigma
  if (instrument_profile_fwhm.size() != 0)
  {
    wavelength_edges[0] = bin_edges.front()[0] + 5.0 * instrument_profile_fwhm.front()/ 2.355;
    wavelength_edges[1] = bin_edges.back()[1] - 5.0 * instrument_profile_fwhm.back()/ 2.355;
  }
  else
  {
    wavelength_edges[0] = bin_edges.front()[0];
    wavelength_edges[1] = bin_edges.back()[1];
  }
}



Observation::~Observation()
{ 
  if (config->use_gpu)
  {
    deleteFromDevice(data_gpu);
    deleteFromDevice(data_error_gpu);
    deleteFromDevice(likelihood_weight_gpu);
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
    moveToDevice(data_gpu, data, true);
    moveToDevice(data_error_gpu, data_error, true);
    moveToDevice(likelihood_weight_gpu, likelihood_weight, true);

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
  double* spectrum_obs,
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
        spectrum_obs,
        is_flux,
        use_filter_response);
    }
    else
    {
      spectral_bands.bandIntegrateSpectrumGPU(
        spectrum_filter_dev,
        spectrum_obs,
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
      spectrum_obs,
      is_flux,
      use_filter_response);
  }
  else
    spectral_bands.bandIntegrateSpectrumGPU(
      spectrum,
      spectrum_obs,
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
