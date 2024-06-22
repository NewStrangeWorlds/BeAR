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
