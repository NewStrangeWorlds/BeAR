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

  //this forward model has a single free general parameter
  nb_general_param = 1;

  //parameter_names is the source of truth; the flat line has one parameter
  parameter_names = {"spectrum_value"};
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
  float* spectrum,
  std::vector<float*>& spectrum_obs)
{
  bool neglect = false;

  const double spectrum_value = parameters[0];

  std::vector<double> spectrum_cpu(spectral_grid->nbSpectralPoints(), spectrum_value);

  moveToDevice(spectrum, spectrum_cpu);

  convertSpectrumToObservationGPU(spectrum, false, spectrum_obs);

  return neglect;
}



FlatLine::~FlatLine()
{

}



}

