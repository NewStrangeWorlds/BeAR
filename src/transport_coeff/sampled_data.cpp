// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "sampled_data.h"

#include "../CUDA_kernels/cross_section_kernels.h"
#include "../CUDA_kernels/data_management_kernels.h"


namespace bear{


void SampledData::deleteSampledData()
{
  std::vector<double>().swap (cross_sections);

  is_sampled = false;
}



void SampledData::sampleCrossSections(
  const std::vector<size_t>& sampling_list_indices, 
  const double species_mass)
{
  if (!data_file.is_loaded) data_file.loadFile();


  cross_sections.assign(sampling_list_indices.size(), 0.0);


  const std::vector<float>& file_data = data_file.cross_sections;

  for(size_t i=0; i<sampling_list_indices.size(); ++i)
  {
    if (sampling_list_indices[i] >= file_data.size()) break;

    cross_sections[i] = file_data[sampling_list_indices[i]];

    //CIA data are stored as log10 values, see CrossSectionFile
    if (data_file.is_data_log)
      cross_sections[i] = std::pow(10.0, cross_sections[i]);
  }


  //convert from opacity in cm2/g to cm2 if neccessary
  //also, apply a small minimum value to allow for interpolation in log later
  for (auto & i : cross_sections)
  {
    if (species_mass > 0) i *= species_mass/6.022140857e23;
    
    if (i < 1e-200) i = 1e-200;
  }
  
  //convert the cross-section to log for performance reasons
  for (auto & i : cross_sections)
    i = std::log10(i);

  
  if (use_gpu)
    moveToDevice(cross_sections_device, cross_sections);


  is_sampled = true;

  data_file.unloadData();
}



SampledData::~SampledData()
{
  if (use_gpu)
    deleteFromDevice(cross_sections_device);
}


}


