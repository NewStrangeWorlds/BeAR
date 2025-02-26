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
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "secondary_eclipse_bb.h"

#include "../../chemistry/chem_species.h"
#include "../../CUDA_kernels/data_management_kernels.h"
#include "../../CUDA_kernels/contribution_function_kernels.h"
#include "../atmosphere/atmosphere.h"


namespace bear{

OccultationBlackBodyPostConfig::OccultationBlackBodyPostConfig (const std::string& folder_path)
{
  const std::string config_file_name = folder_path + "post_process.config";

  readConfigFile(config_file_name);
}


void OccultationBlackBodyPostConfig::readConfigFile(const std::string& file_name)
{
  std::fstream file;
  file.open(file_name.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "\n Post-process config file not found. Using default options!\n\n";

    return;
  }

  std::cout << "\nParameters read from " << file_name << " :\n";

  delete_sampler_files = readBooleanParameter(file, "Delete sampler files");

  save_spectra = readBooleanParameter(file, "Save posterior spectra");

  file.close();
}


//calls the model specific posterior calculations
void OccultationBlackBodyModel::postProcess(
  const std::vector< std::vector<double> >& model_parameter,
  const size_t best_fit_model,
  bool& delete_unused_files)
{
  OccultationBlackBodyPostConfig post_process_config(config->retrieval_folder_path);

  if (post_process_config.delete_sampler_files)
    delete_unused_files = true;

   postProcess(
    post_process_config,
    model_parameter,
    best_fit_model);
}



void OccultationBlackBodyModel::postProcess(
  GenericConfig* post_process_config_,
  const std::vector< std::vector<double> >& model_parameter,
  const size_t best_fit_model,
  bool& delete_unused_files)
{
  OccultationBlackBodyPostConfig* post_process_config = dynamic_cast<OccultationBlackBodyPostConfig*>(post_process_config_);

  if (post_process_config->delete_sampler_files)
    delete_unused_files = true;

  postProcess(
    post_process_config,
    model_parameter,
    best_fit_model);
}



void OccultationBlackBodyModel::postProcess(
  const OccultationBlackBodyPostConfig& post_process_config,
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

