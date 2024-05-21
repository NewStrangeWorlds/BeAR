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


#include "post_process.h"


#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>


#include "../observations/observations.h"
#include "../forward_model/forward_model.h"
#include "../additional/physical_const.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "retrieval.h"



namespace helios{



void PostProcess::postProcessSpectra(std::vector< std::vector<double> >& model_spectrum_bands)
{ 
  const size_t nb_models = model_parameter.size();

  
  model_spectrum_bands.resize(nb_models);
  

  for (size_t i=0; i<nb_models; ++i)
  {
    std::cout << "Postprocess spectrum, model " << i << "\n";

    if (config->use_gpu)
      calcSpectrumGPU(i, model_spectrum_bands[i]);
    else
      calcSpectrum(i, model_spectrum_bands[i]);
  }

}



//Save the high-res spectrum of the best-fit model
//the saved spectrum will have units of W m-2 mu-1 instead of W m-2 cm used internally 
void PostProcess::saveBestFitSpectrum(const std::vector<double>& spectrum)
{ 
  std::vector<double> model_spectrum = forward_model->convertSpectrumToModel(spectrum);

  std::string file_name = config->retrieval_folder_path + "/spectrum_best_fit_hr.dat";

  std::fstream file(file_name.c_str(), std::ios::out);

  
  for (size_t i=0; i<spectrum.size(); ++i)
    file << std::setprecision(10) << std::scientific
         << spectral_grid.wavelength_list[i] << "\t"
         << model_spectrum[i] << "\n";


  file.close();
}




//Calculate the high-res and low-res spectra for a specific set of model parameter from the posterior
void PostProcess::calcSpectrum(const unsigned int model_id, std::vector<double>& model_spectrum_bands)
{ 
  std::vector<double> spectrum_high_res(spectral_grid.nbSpectralPoints(), 0.0);
  model_spectrum_bands.assign(nb_observation_points, 0.0);

  forward_model->calcModel(model_parameter[model_id], spectrum_high_res, model_spectrum_bands);

  //if the current model is the best-fit one, save the high-res spectrum
  if (model_id == best_fit_model)
    saveBestFitSpectrum(spectrum_high_res);
}




void PostProcess::calcSpectrumGPU(const unsigned int model_id, std::vector<double>& model_spectrum_bands)
{
  size_t nb_points = spectral_grid.nbSpectralPoints();

  //pointer to the spectrum on the GPU
  double* spectrum_dev = nullptr;
  allocateOnDevice(spectrum_dev, nb_points);

  //integrate to observational bands and save the low-res spectrum
  double* spectrum_bands_dev = nullptr;
  allocateOnDevice(spectrum_bands_dev, nb_observation_points);

  intializeOnDevice(model_spectrum_gpu, nb_points);


  forward_model->calcModelGPU(model_parameter[model_id], model_spectrum_gpu, spectrum_bands_dev);
  //forward_model->calcModelGPU(model_parameter[model_id], spectrum_dev, spectrum_bands_dev);


  if (model_id == best_fit_model)
  {
    std::vector<double> spectrum(nb_points, 0.0);
    
    moveToHost(model_spectrum_gpu, spectrum);

    saveBestFitSpectrum(spectrum);
  }


  //copy data from the GPU
  model_spectrum_bands.assign(nb_observation_points, 0.0);
  moveToHost(spectrum_bands_dev, model_spectrum_bands);


  deleteFromDevice(spectrum_dev);
  deleteFromDevice(spectrum_bands_dev);
}



}
