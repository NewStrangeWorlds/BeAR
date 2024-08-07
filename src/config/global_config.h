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


#ifndef GLOBAL_CONFIG_H
#define GLOBAL_CONFIG_H

#include <string>

namespace bear {

struct GlobalConfig {
  bool loadConfigFile(std::string retrieval_folder);
  
  std::string forward_model_type = "";

  std::string cross_section_file_path = "";
  std::string wavenumber_file_path = "";
  std::string retrieval_folder_path = "";

  //double spectral_resolution = 0;
  
  unsigned int spectral_disecretisation = 0;
  double const_wavenumber_step = 0;
  double const_wavelength_step = 0;
  double const_spectral_resolution = 0;
  
  bool multinest_ins = false;
  bool multinest_mode_sep = false;
  unsigned int multinest_nb_living_points = 100;
  double multinest_efficiency = 0.8;
  unsigned int multinest_nb_iterations = 0;
  bool multinest_resume = false;
  bool multinest_feedback = false;
  bool multinest_print_iter_values = false;

  unsigned int nb_mpi_processes = 1;
  unsigned int nb_omp_processes = 1;
  
  unsigned int nb_disort_streams = 4;
  bool use_gpu = false;

  bool use_error_inflation = false;
};



}


#endif
