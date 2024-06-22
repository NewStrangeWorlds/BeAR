
#include "sampled_data.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <assert.h>
#include <cmath>


namespace bear{

void CrossSectionFile::loadFile()
{
  std::fstream file;


  file.open(filename.c_str(), std::ios::binary | std::ios::in | std::ios::ate);

  if (file.fail()) std::cout << "cross section file " << filename << " not found! :(((( \n";
  
  assert (!file.fail( ));


  int nb_data_points = file.tellg()/sizeof(float);
  file.seekg(0, std::ios::beg);

  
  cross_sections.resize(nb_data_points);


  for (int i=0; i<nb_data_points; ++i)
  {
    float x;

    file.read((char *) &x, sizeof x);

    cross_sections[i] = x;

    if (is_data_log) cross_sections[i] = std::pow(10.0, cross_sections[i]);
  }


  file.close();

  is_loaded = true;
}


void CrossSectionFile::unloadData()
{
  std::vector<double>().swap (cross_sections);

  is_loaded = false;
}


}