// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "sampled_data.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

#include "../additional/exceptions.h"


namespace bear{

//reads the binary single-precision cross-section file in one go
//the values are kept as stored in the file (linear or log10, see is_data_log);
//the conversion to linear values is done on the sampled subset only
void CrossSectionFile::loadFile()
{
  std::ifstream file(filename.c_str(), std::ios::binary | std::ios::ate);

  if (file.fail())
  {
    throw FileNotFound(std::string ("CrossSectionFile::loadFile"), filename);
  }

  const std::streamoff nb_bytes = file.tellg();
  const size_t nb_data_points = static_cast<size_t>(nb_bytes) / sizeof(float);

  file.seekg(0, std::ios::beg);

  cross_sections.resize(nb_data_points);

  if (nb_data_points > 0)
    file.read(reinterpret_cast<char*>(cross_sections.data()), nb_data_points * sizeof(float));

  if (file.fail())
  {
    throw std::runtime_error("CrossSectionFile::loadFile: error while reading " + filename + "\n");
  }

  is_loaded = true;
}


void CrossSectionFile::unloadData()
{
  std::vector<float>().swap (cross_sections);

  is_loaded = false;
}


}
