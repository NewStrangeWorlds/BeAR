// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef GLOBAL_CONFIG_H
#define GLOBAL_CONFIG_H

#include <string>

namespace bear {

struct GlobalConfig {
  GlobalConfig() {};
  GlobalConfig(const std::string retrieval_folder) {
    loadConfigFile(retrieval_folder);};

  GlobalConfig(
    const bool use_gpu_,
    const std::string forward_model_type_,
    const std::string cross_section_file_path_,
    const std::string spectral_disecretisation_,
    const double resolution_,
    const std::string output_path_,
    const std::string post_output_path_);

  bool loadConfigFile(std::string retrieval_folder);
  
  std::string forward_model_type = "";

  std::string cross_section_file_path = "";
  std::string wavenumber_file_path = "";
  std::string retrieval_folder_path = "";
  std::string output_path = "";
  std::string post_output_path = "";

  //optional overrides for the per-run config file names (see retrieval.toml);
  //let alternate setups live in one directory without renaming
  std::string forward_model_config_file = "forward_model.toml";
  std::string priors_config_file = "priors.config";
  std::string post_process_config_file = "post_process.toml";

  unsigned int spectral_disecretisation = 0;
  double spectral_resolution = 0;
  double spectral_resolution_highres = 0;

  unsigned int nb_mpi_processes = 1;
  unsigned int nb_omp_processes = 0;
  
  unsigned int nb_disort_streams = 4;
  bool use_gpu = false;

  bool use_error_inflation = false;
};



}


#endif
