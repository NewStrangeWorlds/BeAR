// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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


  spectrum = spectral_grid->interpolateToWavelengthGrid(
    grid_wavelengths, 
    data_file.spectrum, 
    false);


  //remove zeros from the spectrum
  for (auto & s : spectrum)
    if (s < 1e-37) s = 1e-37;

  if (use_gpu)
    moveToDevice(spectrum_gpu, spectrum);

  is_sampled = true;

  data_file.unloadData();
}



SampledStellarSpectrum::~SampledStellarSpectrum()
{
  if (spectrum_gpu != nullptr)
    deleteFromDevice(spectrum_gpu);
}

} 
