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



ForwardModelOutput ForwardModel::calcModel(
  const std::vector<double>& physical_parameters,
  const bool return_high_res_spectrum)
{
  ForwardModelOutput output;

  output.spectrum.assign(
    spectral_grid->nbSpectralPoints(), 
    0.0);

  output.spectrum_obs.resize(observations.size());

  for (size_t i=0; i<output.spectrum_obs.size(); ++i)
    output.spectrum_obs[i].assign(
      observations[i].nbPoints(), 
      0.0);


  if (config->use_gpu)
  {
    double* spectrum = nullptr;
    
    allocateOnDevice(
      spectrum, 
      spectral_grid->nbSpectralPoints());

    std::vector<double*> spectrum_obs{
      observations.size(), 
      nullptr};

    for (size_t i=0; i<observations.size(); ++i)
      allocateOnDevice(
        spectrum_obs[i], 
        observations[i].nbPoints());

    output.neglect_model = calcModelGPU(
      physical_parameters, 
      spectrum, 
      spectrum_obs);
    
    if (return_high_res_spectrum)
      moveToHostAndDelete(spectrum, output.spectrum);
    else
      deleteFromDevice(spectrum);

    for (size_t i=0; i<observations.size(); ++i)
      moveToHostAndDelete(spectrum_obs[i], output.spectrum_obs[i]);
  }
  else
  {
    output.neglect_model = calcModelCPU(
      physical_parameters, 
      output.spectrum, 
      output.spectrum_obs);

    if (return_high_res_spectrum == false)
     output.spectrum.clear();
  }

  for (size_t i=0; i<observations.size(); ++i)
  {
    if (observations[i].ascending_wavelengths)
      std::reverse(output.spectrum_obs[i].begin(), output.spectrum_obs[i].end());
  }
  
  
  return output;
}



void ForwardModel::convertSpectrumToObservation(
  const std::vector<double>& spectrum, 
  const bool is_flux,
  std::vector<std::vector<double>>& spectrum_obs)
{
  for (size_t i=0; i<observations.size(); ++i)
    spectrum_obs[i] = observations[i].processModelSpectrum(spectrum, is_flux);
}


void ForwardModel::convertSpectrumToObservationGPU(
  double* spectrum, 
  const bool is_flux,
  std::vector<double*>& spectrum_obs)
{
  unsigned int start_index = 0;

  for (size_t i=0; i<observations.size(); ++i)
  {
    observations[i].processModelSpectrumGPU(
      spectrum, 
      spectrum_obs[i], 
      is_flux);

    start_index += observations[i].spectral_bands.nbBands();
  }
}



void ForwardModel::applyObservationModifier(
  const std::vector<double>& spectrum_modifier_param,
  std::vector<std::vector<double>>& spectrum_obs)
{
  auto param_it = spectrum_modifier_param.begin();

  for (size_t i=0; i<observations.size(); ++i)
  {
    if (observations[i].nb_modifier_param != 0)
    {
      const double spectrum_modifier = *param_it;

      if (spectrum_modifier != 0)
        for (auto & s : spectrum_obs[i])
          s += spectrum_modifier;

      param_it += observations[i].nb_modifier_param;
    }
  }
}



void ForwardModel::applyObservationModifierGPU(
  const std::vector<double>& spectrum_modifier_param,
  std::vector<double*>& spectrum_obs)
{
  auto param_it = spectrum_modifier_param.begin();
  
  for (size_t i=0; i<observations.size(); ++i)
  {
    if (observations[i].nb_modifier_param != 0)
    {
      const double spectrum_modifier = *param_it;

      if (spectrum_modifier != 0)
        observations[i].addShiftToSpectrumGPU(
          spectrum_obs[i], 
          spectrum_modifier);

      param_it += observations[i].nb_modifier_param;
    }
  }
}




void ForwardModel::calcPostProcessSpectra(
  const std::vector< std::vector<double> >& model_parameter,
  const size_t best_fit_model,
  std::vector<std::vector< std::vector<double>>>& model_spectra_obs,
  std::vector<double>& spectrum_best_fit)
{ 
  const size_t nb_models = model_parameter.size();

  model_spectra_obs.resize(nb_models);
  
  std::cout << "\n";

  for (size_t i=0; i<nb_models; ++i)
  {
    std::cout << "\rPostprocess spectra, model " << i << " of " << nb_models << std::flush;
    
    std::vector<double> model_spectrum_high_res;

    calcPostProcessSpectrum(
      model_parameter[i],
      model_spectrum_high_res,
      model_spectra_obs[i]);

    if (i == best_fit_model)
      spectrum_best_fit = model_spectrum_high_res;
  }

  std::cout << "\n";
}



void ForwardModel::calcPostProcessSpectrum(
  const std::vector<double>& model_parameter,
  std::vector<double>& spectrum,
  std::vector<std::vector<double>>& spectrum_obs)
{
  spectrum.assign(nb_spectral_points, 0.0);
  spectrum_obs.assign(observations.size(), std::vector<double>(0.0));

  if (config->use_gpu)
  {
    double* spectrum_gpu = nullptr;
    allocateOnDevice(spectrum_gpu, nb_spectral_points);

    std::vector<double*> spectrum_obs_gpu(observations.size(), nullptr);
    
    for (size_t i=0; i<observations.size(); ++i)
      allocateOnDevice(spectrum_obs_gpu[i], observations[i].nbPoints());

    calcModelGPU(model_parameter, spectrum_gpu, spectrum_obs_gpu);

    moveToHost(spectrum_gpu, spectrum);
    deleteFromDevice(spectrum_gpu);

    for (size_t i=0; i<observations.size(); ++i)
    {
      spectrum_obs[i].assign(observations[i].nbPoints(), 0.0);

      moveToHost(spectrum_obs_gpu[i], spectrum_obs[i]);
      deleteFromDevice(spectrum_obs_gpu[i]);
    }
  }
  else
  {
    calcModelCPU(model_parameter, spectrum, spectrum_obs);
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
  const std::vector<std::vector< std::vector<double>>>& model_spectra_obs)
{
  const size_t nb_models = model_spectra_obs.size();

  for (size_t j=0; j<observations.size(); ++j)
  {
    std::string observation_name = observations[j].observationName();
    std::replace(observation_name.begin(), observation_name.end(), ' ', '_'); 
    
    std::string file_name = config->retrieval_folder_path + "/spectrum_post_" + observation_name + ".dat"; 
    std::fstream file(file_name.c_str(), std::ios::out);
    
    
    std::vector<double> mu = observations[j].spectral_bands.center_wavelengths;
    std::vector<std::vector<double>> spectra = std::vector<std::vector<double>>(nb_models, std::vector<double>(0.0));
    
    for (size_t k=0; k<nb_models; ++k)
      spectra[k] = model_spectra_obs[k][j];
     
    if (observations[j].ascending_wavelengths == true)
    { 
      std::reverse(mu.begin(), mu.end());
      
      for (size_t k=0; k<nb_models; ++k)
        std::reverse(spectra[k].begin(), spectra[k].end());
    }


    for (size_t i=0; i<observations[j].nbPoints(); ++i)
    {
      file << std::setprecision(10) << std::scientific << mu[i];

      for (size_t k=0; k<nb_models; ++k)
         file << "\t" << spectra[k][i];
     
       file << "\n";
    }

    file.close();
  }
}



bool ForwardModel::testCPUvsGPU(const std::vector<double>& parameters)
{ 
  //first we calculate the model on the GPU
  std::cout << "Start test on GPU\n";

  double* spectrum_gpu_dev = nullptr;
  allocateOnDevice(spectrum_gpu_dev, spectral_grid->nbSpectralPoints());

  std::vector<double*> spectrum_obs_gpu_dev{observations.size(), nullptr};

  for (size_t i=0; i<observations.size(); ++i)
    allocateOnDevice(spectrum_obs_gpu_dev[i], observations[i].nbPoints());

  calcModelGPU(parameters, spectrum_gpu_dev, spectrum_obs_gpu_dev);
  
  std::vector<double> spectrum_gpu(spectral_grid->nbSpectralPoints(), 0);
  moveToHost(spectrum_gpu_dev, spectrum_gpu);
  deleteFromDevice(spectrum_gpu_dev);

  std::vector<std::vector<double>> spectrum_obs_gpu(
    observations.size(), 
    std::vector<double>{});

  for (size_t i=0; i<observations.size(); ++i)
  {
    spectrum_obs_gpu[i].assign(observations[i].nbPoints(), 0.0);
    moveToHost(spectrum_obs_gpu_dev[i], spectrum_obs_gpu[i]);
    deleteFromDevice(spectrum_obs_gpu_dev[i]);
  }

  //now we run it on the CPU
  std::cout << "Start test on CPU\n";
  std::vector<double> spectrum_cpu(spectral_grid->nbSpectralPoints(), 0);
  std::vector<std::vector<double>> spectrum_obs_cpu(
    observations.size(), 
    std::vector<double>{});

  calcModelCPU(parameters, spectrum_cpu, spectrum_obs_cpu);
  
  std::cout << "done.\n\n";
  
  std::cout << "Testing high-resolution spectrum:\n";
  std::vector<double> difference_hr(spectral_grid->nbSpectralPoints(), 0);

  for (size_t i=0; i<difference_hr.size(); ++i)
    difference_hr[i] = std::abs(spectrum_cpu[i] - spectrum_gpu[i])/spectrum_cpu[i];

  size_t max_diff_hr_index = std::max_element(difference_hr.begin(),difference_hr.end()) - difference_hr.begin();
  double max_diff_hr = *std::max_element(difference_hr.begin(), difference_hr.end());
  
  std::cout << "Maximum difference of CPU vs GPU: " << max_diff_hr << " at index " << max_diff_hr_index << "\n";
  
  bool test_hr_ok = true;
   if (max_diff_hr*100 > 0.1) test_hr_ok = false;

  std::cout << "Test ok: " << test_hr_ok << "\n\n";
  
  bool test_obs_ok = true;
  std::cout << "Testing binned spectra:\n";
  
  for (size_t i=0; i<observations.size(); ++i)
  {
    std::vector<double> difference_bands(observations[i].nbPoints(), 0);

    for (size_t j=0; j<difference_bands.size(); ++j)
      difference_bands[j] = std::abs(spectrum_obs_cpu[i][j] - spectrum_obs_gpu[i][j])/spectrum_obs_cpu[i][j];

    size_t max_diff_bands_index = std::max_element(difference_bands.begin(),difference_bands.end()) - difference_bands.begin();
    double max_diff_bands = *std::max_element(difference_bands.begin(), difference_bands.end());
    
    std::cout << "Maximum difference of CPU vs GPU for observation " << i << ": " << max_diff_bands << " at index " << max_diff_bands_index << "\n";
    
    bool test_bands_single_ok = true;
    if (max_diff_bands*100 > 0.1) test_bands_single_ok = false;
    
    if (test_bands_single_ok == false)
      test_obs_ok = false;
    
    std::cout << "Test ok: " << test_bands_single_ok << "\n\n";
  }


  bool test_ok = true;

  if (test_hr_ok == false || test_obs_ok == false)
    test_ok = false;
  
  return test_ok;
}


}