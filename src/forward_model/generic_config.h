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
#include <string>

//toml++ is host-only and is not parseable by nvcc; the config parsing helpers
//below are only used from host (.cpp) translation units, so hide them (and the
//toml++ include) from the CUDA compiler. They are non-virtual, so omitting them
//under nvcc does not change the class layout or vtable.
#ifndef __CUDACC__
  #include <toml++/toml.hpp>
#endif

#include "../chemistry/chem_species.h"

namespace bear {


class GenericConfig{
  public:
    virtual void readConfigFile(const std::string& file_name) = 0;

#ifndef __CUDACC__
  protected:
    //parse a TOML file, rethrowing parse errors as InvalidInput with the file path
    static toml::table parseConfigFile(const std::string& file_name);
    //stringify a scalar node (string returned as-is, numbers/bools converted) so
    //that parameter tokens stay std::string exactly like the old parser produced
    static std::string nodeToString(const toml::node& node);
    //error text for a required section that is missing/empty, explaining the
    //common cause (a preceding [section] header capturing the intended top-level key)
    static std::string missingSectionMessage(const std::string& key);
    //fetch a required numeric key that must be present and positive, throwing a
    //clear error (naming the key + section) on a missing/mistyped or invalid value
    static double requirePositive(
      const toml::table& tbl,
      const std::string& section,
      const std::string& key);

    void readAtmosphereConfig(
      const toml::table& cfg,
      size_t& nb_grid_points,
      std::vector<double>& pressure_boundaries);
    //single "model + params" block (temperature, radiative_transfer, stellar_spectrum)
    void readModelBlock(
      const toml::table& cfg,
      const std::string& key,
      std::string& model,
      std::vector<std::string>& parameters,
      const std::string& log_label);
    //variable-length array-of-tables list (chemistry, clouds, modules);
    //an absent or empty array means "none"
    void readModelList(
      const toml::table& cfg,
      const std::string& key,
      std::vector<std::string>& models,
      std::vector<std::vector<std::string>>& parameters,
      bool skip_none,
      bool required,
      const std::string& log_label);
    void readOpacityConfig(
      const toml::table& cfg,
      const std::string& key,
      std::vector<std::string>& opacity_species_symbol,
      std::vector<std::string>& opacity_species_folder,
      bool required);
    bool readBooleanParameter(
      const toml::table& cfg,
      const std::string& key,
      bool default_value);
    std::string readParameter(
      const toml::table& cfg,
      const std::string& key,
      const std::vector<std::string>& allowed_values);
    std::vector<chemical_species_id> readChemicalSpecies(
      const toml::table& cfg,
      const std::string& key);
#endif
};



}


#endif
