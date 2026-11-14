// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _select_cloud_model_h
#define _select_cloud_model_h


#include "cloud_model.h"
#include "grey_cloud_model.h"
#include "kh_cloud_model.h"
#include "power_law_cloud_model.h"
#include "mie_cloud_model.h"

#include "../config/global_config.h"
#include "../additional/exceptions.h"

#include <vector>
#include <string>
#include <algorithm>
#include <memory>


namespace bear {

//definition of the different chemistry modules with an
//identifier, a keyword to be located in the config file
namespace cloud_modules{
  enum id {none, grey, KHnongrey, power_law, mie};
  const std::vector<std::string> description {"none", "grey", "KHnongrey", "power_law", "mie"};
}



inline std::unique_ptr<CloudModel> selectCloudModel(
  const std::string type,
  const std::vector<std::string>& parameters,
  GlobalConfig* config)
{
  //find the corresponding cloud module to the supplied type string
  auto it = std::find(cloud_modules::description.begin(), cloud_modules::description.end(), type);


  //no module is found
  if (it == cloud_modules::description.end())
  {
    std::string error_message = "Cloud module type " + type + " unknown!\n";
    throw InvalidInput(std::string ("forward_model.toml"), error_message);
  }


  //get the id of the chosen module
  cloud_modules::id module_id = static_cast<cloud_modules::id>(std::distance(cloud_modules::description.begin(), it));


  //create the cloud object based on the chosen module
  switch (module_id)
  {
    case cloud_modules::none :
      break;

    case cloud_modules::grey :
      return std::make_unique<GreyCloudModel>(parameters);

    case cloud_modules::KHnongrey :
      return std::make_unique<KHCloudModel>(parameters);

    case cloud_modules::power_law :
      return std::make_unique<PowerLawCloudModel>(parameters);

    case cloud_modules::mie :
      return std::make_unique<MieCloudModel>(parameters, config->retrieval_folder_path);
  }


  return nullptr;
}


}
#endif
