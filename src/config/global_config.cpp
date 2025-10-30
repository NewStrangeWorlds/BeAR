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
  const std::string multinest_output_path_,
  const std::string post_output_path_)
  : use_gpu(use_gpu_),
    forward_model_type(forward_model_type_),
    cross_section_file_path(cross_section_file_path_),
    multinest_output_path(multinest_output_path_),
    post_output_path(post_output_path_)
{
  if (spectral_disecretisation_ != "const_wavenumber" 
   && spectral_disecretisation_ != "const_wavelength" 
   && spectral_disecretisation_ != "const_resolution")
  {
    std::string error_message = "Spectral discretisation parameter: " 
      + spectral_disecretisation_ + " in retrieval.config unknown!\n";
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

  
  std::string file_path = retrieval_folder;
  file_path.append("retrieval.config");

  
  std::fstream file;
  file.open(file_path.c_str(), std::ios::in);

  if (file.fail()) 
  {
    std::cout << "Couldn't open retrieval options file " << file_path << "\n";
    
    return false;
  }

  
  std::cout << "\nParameters found in retrieval.config:\n";

  std::string line;
  std::string input;

  //Header General Config
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  std::cout << "General Program Parameters\n";

  std::getline(file, line);

  file >> input >> line;
  if (input == "Y" || input == "1") use_gpu = true;
  std::cout << "- Use GPU: " << use_gpu << "\n";


  std::getline(file, line);

  file >> nb_omp_processes >> line;

  if (nb_omp_processes == 0)
    nb_omp_processes = omp_get_max_threads();

  std::cout << "- #OpenMP threads: " << nb_omp_processes << "\n";


  //Header Forward Model
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  std::cout << "\n" <<  "General Retrieval Parameters\n";


  std::getline(file, line);
  std::getline(file, line);
  file >> input >> line;
  std::cout << "- Forward model type: " << input << "\n";
  forward_model_type = input;
  
  double spectral_param = 0;
  
  std::getline(file, line);
  file >> input >> spectral_param >> line;
  
  if (input != "const_wavenumber" && input != "const_wavelength" && input != "const_resolution")
  {
    std::string error_message = "Spectral discretisation parameter: " + input + " in retrieval.config unknown!\n";
    throw InvalidInput(std::string ("GlobalConfig::loadConfigFile"), error_message);

    return false;
  }

  if (input == "const_wavenumber")
  {
    spectral_disecretisation = 0;
  }
  else if (input == "const_wavelength")
  {
    spectral_disecretisation = 1;
  }
  else if (input == "const_resolution")
  {
    spectral_disecretisation = 2;
  }

  spectral_resolution = spectral_param;

  std::cout << "- Spectral grid disretisation: " << input << "  " << spectral_param << "\n";


  std::getline(file, line);

  file >> input >> line;
  std::cout << "- Opacity data folder: " << input << "\n";
  cross_section_file_path = input;
  
  if (cross_section_file_path.back() != '/')
    cross_section_file_path.append("/");
  
  wavenumber_file_path = cross_section_file_path + "wavenumber_full.dat";


  std::getline(file, line);

  file >> input >> line;
  if (input == "Y" || input == "1") use_error_inflation = true;
  std::cout << "- Use error inflation prior: " << use_error_inflation << "\n";


  //Header Multinest
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  std::cout << "\n" <<  "Multinest Parameters\n";
  
  std::getline(file, line);

  file >> input >> line;
  if (input == "Y" || input == "1") multinest_ins = true;
  std::cout << "- Importance Nested Sampling: " << multinest_ins << "\n";


  std::getline(file, line);

  file >> input >> line;
  if (input == "Y" || input == "1") multinest_mode_sep = true;
  std::cout << "- Mode separation: " << multinest_mode_sep << "\n";


  std::getline(file, line);

  file >> multinest_nb_living_points >> line;
  std::cout << "- #Living points: " << multinest_nb_living_points << "\n";


  std::getline(file, line);

  file >> multinest_efficiency >> line;
  std::cout << "- Efficiency: " << multinest_efficiency << "\n";

  
  std::getline(file, line);

  file >> multinest_nb_iterations >> line;
  std::cout << "- #Iterations: " << multinest_nb_iterations << "\n";


  std::getline(file, line);

  file >> input >> line;
  if (input == "Y" || input == "1") multinest_resume = true;
  std::cout << "- Resume: " << multinest_resume << "\n";


  std::getline(file, line);

  file >> input >> line;
  if (input == "Y" || input == "1") multinest_feedback = true;
  std::cout << "- Console feedback: " << multinest_feedback << "\n";


  std::getline(file, line);

  file >> input >> line;
  if (input == "Y" || input == "1") multinest_print_iter_values = true;
  std::cout << "- Print iteration values: " << multinest_print_iter_values << "\n";
  
  
  std::cout << "\n";
  
  multinest_output_path = retrieval_folder;
  post_output_path = retrieval_folder;

  return true;
}


}
