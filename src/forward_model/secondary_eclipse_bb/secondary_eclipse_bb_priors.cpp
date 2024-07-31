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
#include <vector>
#include <sstream>

#include "secondary_eclipse_bb.h"

#include "../../retrieval/priors.h"
#include "../../additional/exceptions.h"


namespace bear{


void SecondaryEclipseBlackBodyModel::setPriors(Priors* priors)
{
  const std::string file_name = config->retrieval_folder_path + "priors.config";

  std::vector<std::string> prior_type; 
  std::vector<std::string> prior_description; 
  std::vector<std::vector<std::string>> prior_parameter;


  readPriorConfigFile(file_name, prior_type, prior_description, prior_parameter);

  //check if we have the correct number of piors
  if (prior_type.size() != nb_total_param())
  {
    std::string error_message = 
      "Found " + std::to_string(prior_type.size()) 
      + " priors in priors.config but expected " + std::to_string(nb_total_param()) + "\n";

    throw InvalidInput(std::string ("SecondaryEclipseBlackBodyModel::setPriors"), error_message);
  }


  priors->add(prior_type, prior_description, prior_parameter);
}


}

