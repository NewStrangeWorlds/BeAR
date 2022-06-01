

#include "sampled_data.h"


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


void SampledData::deleteSampledData()
{
  std::vector<double>().swap (cross_sections);

  is_sampled = false;
}



void SampledData::sampleCrossSections(const std::vector<size_t>& sampling_list_indices, const double species_mass)
{
  if (!data_file.is_loaded) data_file.loadFile();


  cross_sections.assign(sampling_list_indices.size(), 0.0);


  for(size_t i=0; i<sampling_list_indices.size(); ++i)
  {
    if (sampling_list_indices[i] > data_file.cross_sections.size()-1) break;

    cross_sections[i] = data_file.cross_sections[sampling_list_indices[i]];
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

}


