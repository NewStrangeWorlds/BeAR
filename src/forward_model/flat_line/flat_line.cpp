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
#include "../../observations/observations.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"
#include "../../CUDA_kernels/data_management_kernels.h"


namespace bear{


FlatLine::FlatLine (
  GlobalConfig* config_,
  SpectralGrid* spectral_grid_,
  std::vector<Observation>& observations_)
    : ForwardModel(config_, spectral_grid_, observations_)
{
  std::cout << "Forward model selected: Flat line\n\n"; 

  //this forward model has three free general parameters
  nb_general_param = 1;
}


//Runs the forward model on the CPU and calculates a high-resolution spectrum
bool FlatLine::calcModelCPU(
  const std::vector<double>& parameter, 
  std::vector<double>& spectrum, 
  std::vector<std::vector<double>>& spectrum_obs)
{
  bool neglect = false;

  const double spectrum_value = parameter[0];

  spectrum.assign(spectral_grid->nbSpectralPoints(), spectrum_value);

  convertSpectrumToObservation(spectrum, false, spectrum_obs);

  return neglect;
}



//run the forward model with the help of the GPU
//the atmospheric structure itself is still done on the CPU
bool FlatLine::calcModelGPU(
  const std::vector<double>& parameters, 
  double* spectrum, 
  std::vector<double*>& spectrum_obs)
{ 
  bool neglect = false;

  const double spectrum_value = parameters[0];

  std::vector<double> spectrum_cpu(spectral_grid->nbSpectralPoints(), spectrum_value);

  moveToDevice(spectrum, spectrum_cpu, false);

  convertSpectrumToObservationGPU(spectrum, false, spectrum_obs);

  return neglect;
}



FlatLine::~FlatLine()
{

}



}

