// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "secondary_eclipse_bb.h"

#include "../../additional/exceptions.h"


namespace bear{


OccultationBlackBodyConfig::OccultationBlackBodyConfig (
  const std::string& folder_path,
  const std::string& file_name)
{
  readConfigFile(folder_path + file_name);
}


OccultationBlackBodyConfig::OccultationBlackBodyConfig (
  const std::string stellar_spectrum_model_,
  const std::vector<std::string>& stellar_model_parameters_)
{
  stellar_spectrum_model = stellar_spectrum_model_;
  stellar_model_parameters = stellar_model_parameters_;
}



void OccultationBlackBodyConfig::readConfigFile(const std::string& file_name)
{
  std::cout << "Parameters read from " << file_name << " :\n";

  toml::table cfg = parseConfigFile(file_name);

  readModelBlock(cfg, "stellar_spectrum",
    stellar_spectrum_model, stellar_model_parameters, "Stellar spectrum model");
}



}

