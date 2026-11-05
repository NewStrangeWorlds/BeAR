// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "emission.h"

#include "../../retrieval/retrieval.h"
#include "../../CUDA_kernels/data_management_kernels.h"


namespace bear{


bool EmissionModel::testModel(const std::vector<double>& parameters)
{
  bool test_ok = false;
  
  test_ok = testCPUvsGPU(parameters);

  return test_ok;
}



}