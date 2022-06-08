/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2022 Daniel Kitzmann
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


#ifndef _select_cloud_model_h
#define _select_cloud_model_h


#include "cloud_model.h"
#include "grey_cloud_model.h"

#include "../config/global_config.h"
#include "../additional/exceptions.h"

#include <vector>
#include <string>
#include <algorithm>


namespace helios {

//definition of the different chemistry modules with an
//identifier, a keyword to be located in the config file
namespace cloud_modules{
  enum id {none, grey, KHnongrey}; 
  const std::vector<std::string> description {"none", "grey", "KHnongrey"};
}



inline CloudModel* selectCloudModel(const std::string type, const std::vector<std::string>& parameters)
{
  //find the corresponding cloud module to the supplied type string
  auto it = std::find(cloud_modules::description.begin(), cloud_modules::description.end(), type);


  //no module is found
  if (it == cloud_modules::description.end())
  {
    std::string error_message = "Cloud module type " + type + " unknown!\n";
    throw ExceptionInvalidInput(std::string ("forward_model.config"), error_message);
  }


  //get the id of the chosen module
  cloud_modules::id module_id = static_cast<cloud_modules::id>(std::distance(cloud_modules::description.begin(), it));


  //create the radiative transfer object based on the chosen module
  CloudModel* cloud_model = nullptr;

  switch (module_id)
  {
    case cloud_modules::none :
      break;

    case cloud_modules::grey :
      {
        GreyCloudModel* model = new GreyCloudModel(parameters);
        cloud_model = model;  
      }
      break;
  }


  return cloud_model;
}


}
#endif

