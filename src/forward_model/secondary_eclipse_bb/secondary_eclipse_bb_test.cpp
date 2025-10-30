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

#include "secondary_eclipse_bb.h"

#include "../../CUDA_kernels/data_management_kernels.h"


namespace bear{


bool OccultationBlackBodyModel::testModel(
  const std::vector<double>& parameters)
{
  bool test_ok = false;
  
  test_ok = testCPUvsGPU(parameters);

  return test_ok;
}


}