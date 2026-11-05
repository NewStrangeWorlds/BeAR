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
#include "../../retrieval/priors.h"
#include "../../observations/observations.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"
#include "../../CUDA_kernels/data_management_kernels.h"


namespace bear{

FlatLinePostProcessConfig::FlatLinePostProcessConfig (
  const std::string& folder_path,
  const std::string& file_name)
{
  readConfigFile(folder_path + file_name);
}


void FlatLinePostProcessConfig::readConfigFile(const std::string& file_name)
{
  std::ifstream test_file(file_name.c_str());

  if (!test_file.good())
  {
    std::cout << "\n Post-process config file not found. Using default options!\n\n";

    return;
  }
  test_file.close();

  std::cout << "\nParameters read from " << file_name << " :\n";

  toml::table cfg = parseConfigFile(file_name);

  delete_sampler_files = readBooleanParameter(cfg, "delete_sampler_files", delete_sampler_files);
  save_spectra = readBooleanParameter(cfg, "save_spectra", save_spectra);
}



void FlatLine::postProcess(
  GenericConfig* post_process_config_,
  const std::vector< std::vector<double> >& model_parameter,
  const size_t best_fit_model,
  bool& delete_unused_files)
{
  FlatLinePostProcessConfig post_process_config = 
    dynamic_cast<FlatLinePostProcessConfig&>(*post_process_config_);

  if (post_process_config.delete_sampler_files)
    delete_unused_files = true;

  postProcess(
    post_process_config,
    model_parameter,
    best_fit_model);
}



void FlatLine::postProcess(
  const std::vector< std::vector<double> >& model_parameter,
  const size_t best_fit_model,
  bool& delete_unused_files)
{
  FlatLinePostProcessConfig post_process_config(
    config->retrieval_folder_path, config->post_process_config_file);

  if (post_process_config.delete_sampler_files)
    delete_unused_files = true;

  if (post_process_config.save_spectra)
  {
    std::vector<std::vector<std::vector<double>>> model_spectra_obs;
    std::vector<double> model_spectrum_best_fit;
  
    calcPostProcessSpectra(
      model_parameter, 
      best_fit_model, 
      model_spectra_obs,
      model_spectrum_best_fit);
    
    saveBestFitSpectrum(model_spectrum_best_fit);
    savePostProcessSpectra(model_spectra_obs);
  }
}



void FlatLine::postProcess(
  const FlatLinePostProcessConfig& post_process_config,
  const std::vector< std::vector<double> >& model_parameter,
  const size_t best_fit_model)
{
  if (post_process_config.save_spectra)
  {
    std::vector<std::vector<std::vector<double>>> model_spectra_obs;
    std::vector<double> model_spectrum_best_fit;
  
    calcPostProcessSpectra(
      model_parameter, 
      best_fit_model, 
      model_spectra_obs,
      model_spectrum_best_fit);
    
    saveBestFitSpectrum(model_spectrum_best_fit);
    savePostProcessSpectra(model_spectra_obs);
  }
}



}

