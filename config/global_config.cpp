/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2019 Daniel Kitzmann
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


#include "global_config.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include <stdlib.h>

namespace helios {


bool GlobalConfig::loadGlobalConfig(std::string retrieval_folder)
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

  



  std::string line;
  std::string input;

  //Header General Config
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  std::cout << "\n" << "General Parameter\n";

  std::getline(file, line);

  file >> input >> line;
  if (input == "Y" || input == "1") use_gpu = true;
  std::cout << "- Use GPU: " << use_gpu << "\n";


  std::getline(file, line);

  file >> nb_omp_processes >> line;
  std::cout << "- #OpenMP threads: " << nb_omp_processes << "\n";

  //Header Forward Model
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  std::cout << "\n" << "Forward model Parameter\n";
  
  std::getline(file, line);
  
  file >> nb_atmosphere_levels >> line;
  std::cout << "- #Atmosphere levels: " << nb_atmosphere_levels << "\n";


  std::getline(file, line);

  file >> atmosphere_top_pressure >> line;
  std::cout << "- Top of atmosphere pressure: " << atmosphere_top_pressure << "\n";


  std::getline(file, line);

  file >> atmosphere_bottom_pressure >> line;
  std::cout << "- Bottom of atmosphere pressure: " << atmosphere_bottom_pressure << "\n";


  std::getline(file, line);

  file >> spectral_resolution >> line;
  std::cout << "- Spectral resolution: " << spectral_resolution << "\n";

  std::getline(file, line);

  file >> input >> line;
  std::cout << "- Opacity data folder: " << input << "\n";
  cross_section_file_path = input + "molecules/";
  wavenumber_file_path = input + "wavenumber_full.dat";


  std::getline(file, line);

  file >> fastchem_parameter_file >> line;
  std::cout << "- FastChem parameter file: " << fastchem_parameter_file << "\n";
  fastchem_parameter_file = retrieval_folder_path + "/" + fastchem_parameter_file;


  std::getline(file, line);

  file >> nb_temperature_elements >> line;
  std::cout << "- Number of Temperature elements: " << nb_temperature_elements << "\n";


  std::getline(file, line);

  file >> temperature_poly_degree >> line;
  std::cout << "- Temperature polynomial degree: " << temperature_poly_degree << "\n";


  //Header Multinest
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  std::cout << "\n" <<  "Multinest Parameter\n";
  
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


  return true;
}


}
