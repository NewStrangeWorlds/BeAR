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


#include "global_config.h"

#include <toml++/toml.hpp>

#include "../additional/exceptions.h"

#include <exception>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <omp.h>
#include <stdlib.h>

namespace bear {


GlobalConfig::GlobalConfig(
  const bool use_gpu_,
  const std::string forward_model_type_,
  const std::string cross_section_file_path_,
  const std::string spectral_disecretisation_,
  const double resolution_,
  const std::string output_path_,
  const std::string post_output_path_)
  : use_gpu(use_gpu_),
    forward_model_type(forward_model_type_),
    cross_section_file_path(cross_section_file_path_),
    output_path(output_path_),
    post_output_path(post_output_path_)
{
  if (spectral_disecretisation_ != "const_wavenumber" 
   && spectral_disecretisation_ != "const_wavelength" 
   && spectral_disecretisation_ != "const_resolution")
  {
    std::string error_message = "Spectral discretisation parameter: "
      + spectral_disecretisation_ + " in retrieval.toml unknown!\n";
    throw InvalidInput(std::string ("GlobalConfig::GlobalConfig"), error_message);
  }

  if (spectral_disecretisation_ == "const_wavenumber")
  {
    spectral_disecretisation = 0;
  }
  else if (spectral_disecretisation_ == "const_wavelength")
  {
    spectral_disecretisation = 1;
  }
  else if (spectral_disecretisation_ == "const_resolution")
  {
    spectral_disecretisation = 2;
  }

  spectral_resolution = resolution_;
}



bool GlobalConfig::loadConfigFile(std::string retrieval_folder)
{
  if (retrieval_folder.back() != '/')
    retrieval_folder.append("/");

  retrieval_folder_path = retrieval_folder;

  const std::string file_path = retrieval_folder + "retrieval.toml";

  toml::table cfg;

  try
  {
    cfg = toml::parse_file(file_path);
  }
  catch (const toml::parse_error& e)
  {
    std::ostringstream error_message;
    error_message << "Error parsing " << file_path << ": "
                  << e.description() << " (at " << e.source().begin << ")\n";
    throw InvalidInput(std::string ("GlobalConfig::loadConfigFile"), error_message.str());
  }

  std::cout << "\nParameters found in retrieval.toml:\n";
  std::cout << "General Program Parameters\n";

  use_gpu = cfg["general"]["use_gpu"].value_or(false);
  std::cout << "- Use GPU: " << use_gpu << "\n";

  nb_omp_processes = cfg["general"]["nb_omp_threads"].value_or(0);

  if (nb_omp_processes == 0)
    nb_omp_processes = omp_get_max_threads();

  std::cout << "- #OpenMP threads: " << nb_omp_processes << "\n";


  std::cout << "\n" << "General Retrieval Parameters\n";

  auto retrieval = cfg["retrieval"];

  forward_model_type = retrieval["forward_model_type"].value_or(std::string(""));
  std::cout << "- Forward model type: " << forward_model_type << "\n";

  const std::string discretisation =
    retrieval["spectral_discretisation"].value_or(std::string(""));
  spectral_resolution = retrieval["spectral_resolution"].value_or(0.0);

  if (discretisation == "const_wavenumber")
    spectral_disecretisation = 0;
  else if (discretisation == "const_wavelength")
    spectral_disecretisation = 1;
  else if (discretisation == "const_resolution")
    spectral_disecretisation = 2;
  else
  {
    std::string error_message = "Spectral discretisation parameter: " + discretisation
      + " in retrieval.toml unknown!\n";
    throw InvalidInput(std::string ("GlobalConfig::loadConfigFile"), error_message);
  }

  std::cout << "- Spectral grid disretisation: " << discretisation
            << "  " << spectral_resolution << "\n";

  cross_section_file_path = retrieval["opacity_data_folder"].value_or(std::string(""));
  std::cout << "- Opacity data folder: " << cross_section_file_path << "\n";

  if (cross_section_file_path.empty() || cross_section_file_path.back() != '/')
    cross_section_file_path.append("/");

  wavenumber_file_path = cross_section_file_path + "wavenumber_full.dat";

  use_error_inflation = retrieval["use_error_inflation"].value_or(false);
  std::cout << "- Use error inflation prior: " << use_error_inflation << "\n";

  //optional high-resolution spectral grid
  spectral_resolution_highres =
    retrieval["spectral_resolution_highres"].value_or(0.0);

  if (spectral_resolution_highres > 0)
    std::cout << "- High-resolution spectral grid resolution: "
              << spectral_resolution_highres << "\n";

  //optional overrides for the per-run config file names
  forward_model_config_file =
    retrieval["forward_model_config"].value_or(forward_model_config_file);
  priors_config_file =
    retrieval["priors_config"].value_or(priors_config_file);
  post_process_config_file =
    retrieval["post_process_config"].value_or(post_process_config_file);

  std::cout << "\n";

  output_path = retrieval_folder;
  post_output_path = retrieval_folder;

  return true;
}


}
