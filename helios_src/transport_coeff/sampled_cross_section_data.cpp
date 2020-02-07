/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
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


#include "transport_coeff_single_species.h"


#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <cmath>


#include "../CUDA_kernels/cross_section_kernels.h"
#include "../CUDA_kernels/data_management_kernels.h"



namespace helios{


SampledCrossSectionData::~SampledCrossSectionData()
{
  
  if (use_gpu)
    deleteFromDevice(cross_sections_device);

}



double SampledCrossSectionData::getCrossSection(unsigned int wavelength_index)
{
  if (is_sampled)
    return cross_sections[wavelength_index];
  else
    std::cout << "sampled data not available: T: " << temperature << " P: " << pressure << std::endl;

  return 0;
}


void SampledCrossSectionData::getCrossSections(std::vector<double>& data_vector)
{
  if (is_sampled)
    data_vector = cross_sections;
  else
    std::cout << "sampled data not available: T: " << temperature << " P: " << pressure << std::endl;
}


void SampledCrossSectionData::init(const double temperature_data, const double pressure_data, const bool gpu_usage)
{
  temperature = temperature_data;
  pressure = pressure_data;
  is_sampled = false;
  
  use_gpu = gpu_usage;
}


void SampledCrossSectionData::deleteSampledData()
{
  std::vector<double>().swap (cross_sections);

  is_sampled = false;
}




void SampledCrossSectionData::sampleCrossSections(CrossSectionFile& input_file, const std::vector<size_t>& sampling_list_indices)
{
  if (!input_file.is_loaded) input_file.loadFile();

  pressure = input_file.pressure;
  temperature = input_file.temperature;


  cross_sections.assign(sampling_list_indices.size(), 0.0);


  for(size_t i=0; i<sampling_list_indices.size(); ++i)
  {
    if (sampling_list_indices[i] > input_file.cross_sections.size()-1) break;


    cross_sections[i] = input_file.cross_sections[sampling_list_indices[i]];
  }


  //apply a small minimum value to allow for interpolation in log later
  for (size_t i=0; i<cross_sections.size(); ++i)
    if (cross_sections[i] < 1e-200) cross_sections[i] = 1e-200;

  
  if (use_gpu)
  {
    //for the GPU, we directly convert the cross-section to log10
    //in order to avoid having to do that constantly on the GPU
    std::vector<double> log_cross_sections = cross_sections;
    for (size_t i=0; i<log_cross_sections.size(); ++i)
      log_cross_sections[i] = std::log10(log_cross_sections[i]);


    moveToDevice(cross_sections_device, log_cross_sections);
  }


  is_sampled = true;


  input_file.unloadData();


  //for (unsigned int i=0; i<cross_sections.size(); i++)
    //std::cout << sampling_list_indices[i] << "\t" << cross_sections[i] << "\n"; 
}




}


