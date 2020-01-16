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


#include <iostream>


#include "../spectral_grid/spectral_grid.h"
#include "../config/global_config.h"
#include "../retrieval/retrieval.h"
#include "../retrieval/post_process.h"
#include <omp.h>


void doRetrieval(helios::GlobalConfig& config)
{
  
  std::cout << "Starting retrieval\n\n";
  
  helios::Retrieval retrieval(&config);
  retrieval.doRetrieval();

}


void doPostProcess(helios::GlobalConfig& config)
{
  
  std::cout << "Starting post processing\n\n";
  
  helios::PostProcess post_process(&config);
  post_process.doRetrieval();

}


int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    std::cout << "config file command line parameter missing!\n";

    return 1;
  }

  
  std::string retrieval_folder = argv[1];

  
  helios::GlobalConfig config;
  if (config.loadGlobalConfig(retrieval_folder) == false)
    return 1;


  omp_set_num_threads(config.nbOMPProcesses());

 
  doRetrieval(config);

  doPostProcess(config);


  std::cout << "Helios-r finished" << std::endl;

  return 0;
}
