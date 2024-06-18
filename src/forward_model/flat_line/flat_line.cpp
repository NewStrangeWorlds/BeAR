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

#include "flat_line.h"

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../retrieval/priors.h"
#include "../../observations/observations.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"
#include "../../CUDA_kernels/data_management_kernels.h"


namespace helios{


FlatLine::FlatLine (
  Priors* priors_,
  GlobalConfig* config_,
  SpectralGrid* spectral_grid_,
  std::vector<Observation>& observations_)
    : config(config_)
    , spectral_grid(spectral_grid_)
    , observations(observations_)
{
  std::cout << "Forward model selected: Flat line\n\n"; 

  //this forward model has three free general parameters
  nb_general_param = 1;

  setPriors(priors_);

  for (auto & i : observations)
    nb_observation_points += i.nbPoints();
}


//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool FlatLine::calcModel(
  const std::vector<double>& parameter, 
  std::vector<double>& spectrum, 
  std::vector<double>& model_spectrum_bands)
{
  bool neglect = false;

  const double spectrum_value = parameter[0];

  spectrum.assign(spectral_grid->nbSpectralPoints(), spectrum_value);

  postProcessSpectrum(spectrum, model_spectrum_bands);

  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool FlatLine::calcModelGPU(
  const std::vector<double>& parameter, 
  double* model_spectrum_gpu, 
  double* model_spectrum_bands)
{ 
  bool neglect = false;

  const double spectrum_value = parameter[0];
  
  std::vector<double> spectrum(spectral_grid->nbSpectralPoints(), spectrum_value);
  
  moveToDevice(model_spectrum_gpu, spectrum, false);

  postProcessSpectrumGPU(model_spectrum_gpu, model_spectrum_bands);

  return neglect;
}



//integrate the high-res spectrum to observational bands
//and convolve if necessary 
void FlatLine::postProcessSpectrum(
  std::vector<double>& model_spectrum, 
  std::vector<double>& model_spectrum_bands)
{
  model_spectrum_bands.assign(nb_observation_points, 0.0);
  
  std::vector<double>::iterator it = model_spectrum_bands.begin();

  for (size_t i=0; i<observations.size(); ++i)
  {
    const bool is_flux = false;

    std::vector<double> observation_bands = 
      observations[i].processModelSpectrum(model_spectrum, is_flux);

    //copy the band-integrated values for this observation into the global
    //vector of all band-integrated points, model_spectrum_bands
    std::copy(observation_bands.begin(), observation_bands.end(), it);
    it += observation_bands.size();
  }
}


//integrate the high-res spectrum to observational bands
//and convolve if necessary 
void FlatLine::postProcessSpectrumGPU(
  double* model_spectrum_gpu, 
  double* model_spectrum_bands)
{
  unsigned int start_index = 0;
  for (size_t i=0; i<observations.size(); ++i)
  {
    const bool is_flux = false;

    observations[i].processModelSpectrumGPU(
      model_spectrum_gpu, 
      model_spectrum_bands, 
      start_index, 
      is_flux);

    start_index += observations[i].spectral_bands.nbBands();
  }
}


//Postprocess
void FlatLine::postProcess(
  const std::vector< std::vector<double> >& model_parameter, 
  const std::vector< std::vector<double> >& model_spectrum_bands, 
  const size_t best_fit_model)
{
  //nothing to do here...

}


std::vector<double> FlatLine::convertSpectrumToModel(const std::vector<double>& spectrum)
{
  //the high-res spectrum is already a flat line
  std::vector<double> model_spectrum = spectrum;

  return model_spectrum;
}


FlatLine::~FlatLine()
{

}



}

