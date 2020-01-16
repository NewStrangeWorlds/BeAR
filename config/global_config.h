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


#ifndef GLOBAL_CONFIG_H
#define GLOBAL_CONFIG_H

#include <string>

namespace helios {

class GlobalConfig {
  public:
    GlobalConfig () {}

    std::string crossSectionFilePath() {return cross_section_file_path;}
	  std::string wavenumberFilePath() {return wavenumber_file_path;}

	  std::string retrievalFolder() {return retrieval_folder_path;}

    unsigned int nbOMPProcesses() {return nb_omp_processes;}
	  bool loadGlobalConfig(std::string retrieval_folder);
	  bool useGPU() {return use_gpu;}

    std::string cross_section_file_path = "/gpfs/homefs/csh/kitzmann/data/cross_sections_helios_new/molecules/";
    std::string wavenumber_file_path = "/gpfs/homefs/csh/kitzmann/data/cross_sections_helios_new/wavenumber_full.dat";
    std::string retrieval_folder_path = "";

    unsigned int nb_atmosphere_levels = 0;
    double spectral_resolution = 0;
    std::string fastchem_parameter_file;
    unsigned int nb_temperature_elements = 4;
    unsigned int temperature_poly_degree = 1;
    double atmosphere_top_pressure = 100;
    double atmosphere_bottom_pressure = 1e-3;

    bool multinest_ins = false;
    bool multinest_mode_sep = false;
    unsigned int multinest_nb_living_points = 100;
    double multinest_efficiency = 0.8;
    unsigned int multinest_nb_iterations = 0;
    bool multinest_resume = false;
    bool multinest_feedback = false;
    bool multinest_print_iter_values = false;

  private:
    unsigned int nb_mpi_processes = 1;
    unsigned int nb_omp_processes = 1;

  
    unsigned int nb_disort_streams = 4;
    bool use_gpu = false;
};



}


#endif
