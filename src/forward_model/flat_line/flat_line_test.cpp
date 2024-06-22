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


#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "flat_line.h"

#include "../../CUDA_kernels/data_management_kernels.h"


namespace bear{


bool FlatLine::testModel(const std::vector<double>& parameter, double* model_spectrum_gpu)
{
  bool test_ok = false;
  
  test_ok = testCPUvsGPU(parameter, model_spectrum_gpu);

  return test_ok;
}



bool FlatLine::testCPUvsGPU(const std::vector<double>& parameter, double* model_spectrum_gpu)
{ 
  //first we calculate the model on the GPU
  std::cout << "Start test on GPU\n";
  //pointer to the spectrum on the GPU
  double* spectrum_bands_dev = nullptr;
  allocateOnDevice(spectrum_bands_dev, nb_observation_points);

  //intialise the high-res spectrum on the GPU (set it to 0)
  intializeOnDevice(model_spectrum_gpu, spectral_grid->nbSpectralPoints());

  calcModelGPU(parameter, model_spectrum_gpu, spectrum_bands_dev);
  
  std::vector<double> spectrum_bands_gpu(nb_observation_points, 0);
  moveToHost(spectrum_bands_dev, spectrum_bands_gpu);


  //now we run it on the CPU
  std::cout << "Start test on CPU\n";
  std::vector<double> spectrum_cpu, spectrum_bands_cpu;

  calcModel(parameter, spectrum_cpu, spectrum_bands_cpu);
  
  std::cout << "done.\n";
  std::vector<double> difference(nb_observation_points, 0);

  for (size_t i=0; i<difference.size(); ++i)
    difference[i] = std::abs(spectrum_bands_cpu[i] - spectrum_bands_gpu[i])/spectrum_bands_cpu[i];

  size_t max_diff_index = std::max_element(difference.begin(),difference.end()) - difference.begin();
  double max_diff = *std::max_element(difference.begin(), difference.end());
  

  std::cout << "Maximum difference of CPU vs GPU: " << max_diff << " at index " << max_diff_index << "\n";
  
  bool test_ok = true;
   if (max_diff*100 > 0.1) test_ok = false;

  std::cout << "Test ok: " << test_ok << "\n\n";

  //for (size_t i=0; i<spectrum_bands_cpu.size(); ++i)
    //std::cout << i << "\t" << spectrum_bands_cpu[i] << "\t" << spectrum_bands_gpu[i] << "\n";
  
  return test_ok;
}





}