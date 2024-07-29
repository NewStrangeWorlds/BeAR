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

#include "forward_model.h"

#include "../config/global_config.h"
#include "../spectral_grid/spectral_grid.h"
#include "../retrieval/priors.h"
#include "../observations/observations.h"
#include "../additional/exceptions.h"



namespace bear{

ForwardModel::ForwardModel (
  GlobalConfig* config_, 
  SpectralGrid* spectral_grid_,
  std::vector<Observation>& observations_) 
    : config(config_)
    , spectral_grid(spectral_grid_)
    , observations(observations_) 
{
  for (auto & i : observations)
  {
    nb_observation_points += i.nbPoints();
    nb_spectrum_modifier_param += i.nb_modifier_param;
  }

  nb_spectral_points = spectral_grid->nbSpectralPoints();
}



void ForwardModel::readPriorConfigFile(
  const std::string& file_path, 
  std::vector<std::string>& prior_type, 
  std::vector<std::string>& prior_description, 
  std::vector<std::vector<std::string>>& prior_parameter)
{
  std::fstream file;
  file.open(file_path.c_str(), std::ios::in);

  if (file.fail())  
    throw FileNotFound(std::string ("ForwardModel::readPriorConfigFile"), file_path);


  std::string line;

  while (std::getline(file, line))
  {
    std::istringstream input(line);

    std::string type, description;
    std::vector<std::string> parameter;

    input >> type >> description;

    std::string single_parameter;

    while (input >> single_parameter)
      parameter.push_back(single_parameter);

    prior_type.push_back(type);
    prior_description.push_back(description);
    prior_parameter.push_back(parameter);
  }

  file.close();
}



void ForwardModel::calcPostProcessSpectra(
  const std::vector< std::vector<double> >& model_parameter,
  const size_t best_fit_model,
  std::vector< std::vector<double> >& model_spectrum_bands)
{
  const size_t nb_models = model_parameter.size();

  model_spectrum_bands.resize(nb_models);
  
  std::cout << "\n";

  for (size_t i=0; i<nb_models; ++i)
  {
    std::cout << "\rPostprocess spectra, model " << i << " of " << nb_models << std::flush;
    
    std::vector<double> model_spectrum_high_res;

    calcPostProcessSpectrum(
      model_parameter[i],
      model_spectrum_high_res,
      model_spectrum_bands[i]);

    if (i == best_fit_model)
      saveBestFitSpectrum(model_spectrum_high_res);
  }

  std::cout << "\n";

  savePostProcessSpectra(model_spectrum_bands);
}



void ForwardModel::calcPostProcessSpectrum(
  const std::vector<double>& model_parameter,
  std::vector<double>& model_spectrum,
  std::vector<double>& model_spectrum_bands)
{
  model_spectrum.assign(nb_spectral_points, 0.0);
  model_spectrum_bands.assign(nb_observation_points, 0.0);

  if (config->use_gpu)
  {
    double* spectrum_dev = nullptr;
    allocateOnDevice(spectrum_dev, nb_spectral_points);

    double* spectrum_bands_dev = nullptr;
    allocateOnDevice(spectrum_bands_dev, nb_observation_points);

    calcModelGPU(model_parameter, spectrum_dev, spectrum_bands_dev);

    moveToHost(spectrum_dev, model_spectrum);
    moveToHost(spectrum_bands_dev, model_spectrum_bands);

    deleteFromDevice(spectrum_dev);
    deleteFromDevice(spectrum_bands_dev);
  }
  else
  {
    calcModel(model_parameter, model_spectrum, model_spectrum_bands);
  }
}



void ForwardModel::saveBestFitSpectrum(const std::vector<double>& spectrum)
{ 
  std::string file_name = config->retrieval_folder_path + "/spectrum_best_fit_hr.dat";

  std::fstream file(file_name.c_str(), std::ios::out);

  for (size_t i=0; i<spectrum.size(); ++i)
    file << std::setprecision(10) << std::scientific
         << spectral_grid->wavelength_list[i] << "\t"
         << spectrum[i] << "\n";

  file.close();
}


void ForwardModel::savePostProcessSpectra(
  const std::vector< std::vector<double> >& model_spectrum_bands)
{
  const size_t nb_models = model_spectrum_bands.size();

  unsigned int band_index = 0;

  for (size_t j=0; j<observations.size(); ++j)
  {
    std::string observation_name = observations[j].observationName();
    std::replace(observation_name.begin(), observation_name.end(), ' ', '_'); 
    
    std::string file_name = config->retrieval_folder_path + "/spectrum_post_" + observation_name + ".dat"; 
    std::fstream file(file_name.c_str(), std::ios::out);

    for (size_t i=0; i<observations[j].spectral_bands.nbBands(); ++i)
    {
      file << std::setprecision(10) << std::scientific << observations[j].spectral_bands.center_wavelengths[i];
      
      for (size_t k=0; k<nb_models; ++k)
        file << "\t" << model_spectrum_bands[k][band_index + i];
     
      file << "\n";
    }

    file.close();

    band_index += observations[j].spectral_bands.nbBands();
  }

}


}