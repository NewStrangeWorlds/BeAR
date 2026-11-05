// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <assert.h>

#include "stellar_spectrum_grid.h"


namespace bear{


void SpectrumFile::loadFile()
{
  std::fstream file;


  file.open(file_path.c_str(), std::ios::binary | std::ios::in | std::ios::ate);

  if (file.fail()) std::cout << "cross section file " << file_path << " not found! :(((( \n";
  
  assert (!file.fail( ));


  int nb_data_points = file.tellg()/sizeof(float);
  file.seekg(0, std::ios::beg);

  
  spectrum.resize(nb_data_points);


  for (int i=0; i<nb_data_points; ++i)
  {
    float x;

    file.read((char *) &x, sizeof x);

    spectrum[i] = x;
  }


  file.close();

  is_loaded = true;
}


void SpectrumFile::unloadData()
{
  std::vector<double>().swap (spectrum);

  is_loaded = false;
}


} 
