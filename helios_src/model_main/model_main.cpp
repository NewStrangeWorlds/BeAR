/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
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


#include "command_line_param.h"
#include "../spectral_grid/spectral_grid.h"
#include "../config/global_config.h"
#include "../retrieval/retrieval.h"
#include "../retrieval/post_process.h"
#include <omp.h>
#include <csignal>
#include <cstdlib>



bool doRetrieval(helios::GlobalConfig& config)
{
  std::cout << "Starting retrieval\n\n";
  
  helios::Retrieval retrieval(&config);
  return retrieval.doRetrieval();
}



bool doPostProcess(helios::GlobalConfig& config)
{
  std::cout << "Starting post processing\n\n";
  
  helios::PostProcess post_process(&config);
  return post_process.doRetrieval();
}




int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    std::cout << "config file command line parameter missing!\n";

    return 1;
  }
    
  
  std::string retrieval_folder = argv[1];


  CommandLineParam command_line_param;
  
  parseCommandLine(argc, argv, command_line_param);

  
  helios::GlobalConfig config;
  if (config.loadConfigFile(retrieval_folder) == false)
    return 1;
  

  if (command_line_param.multinest_resume == true)
    config.multinest_resume = true; 
  

  omp_set_num_threads(config.nb_omp_processes);

  
  bool retrieval_success = true;
  
  if (!command_line_param.only_post_process) 
    retrieval_success = doRetrieval(config);

  if (retrieval_success)
    retrieval_success = doPostProcess(config);



  
  if (retrieval_success)
  {
    std::cout << "\nHelios-r finished\n" << std::endl;
    return 0;
  }
  else
  {
    std::cout << "\nHelios-r finished with errors :(\n" << std::endl;
    return 1;
  }
  
  
  return 0;
}
