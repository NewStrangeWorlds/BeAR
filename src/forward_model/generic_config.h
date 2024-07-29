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


#ifndef _generic_config_h
#define _generic_config_h


#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>


namespace bear {


class GenericConfig{
  public:
    virtual void readConfigFile(const std::string& file_name) = 0;
  protected:
    void readAtmosphereConfig(
      std::fstream& file,
      size_t& nb_grid_points,
      std::vector<double>& pressure_boundaries);
    void readTemperatureConfig(
      std::fstream& file,
      std::string& model,
      std::vector<std::string>& parameters);
    void readCloudConfig(
      std::fstream& file,
      std::vector<std::string>& models,
      std::vector<std::vector<std::string>>& parameters);
    void readModuleConfig(
      std::fstream& file,
      std::vector<std::string>& modules,
      std::vector<std::vector<std::string>>& parameters);
    void readChemistryConfig(
      std::fstream& file,
      std::vector<std::string>& models,
      std::vector<std::vector<std::string>>& parameters);
    void readOpacityConfig(
      std::fstream& file,
      std::vector<std::string>& opacity_species_symbol,
      std::vector<std::string>& opacity_species_folder);
    bool readBooleanParameter(
      std::fstream& file,
      const std::string& parameter_name);
    std::string readParameter(
      std::fstream& file,
      const std::string& parameter_name,
      const std::vector<std::string>& allowed_values);
};



}


#endif

