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


#ifndef _data_management_kernels_h
#define _data_management_kernels_h


#include <vector>

namespace bear{


void allocateOnDevice(double*& device_data, size_t nb_double_values);
void deleteFromDevice(double*& device_data);
void moveToDevice(double*& device_data, std::vector<double>& host_data, const bool alloc_memory);
void moveToDevice(double*& device_data, std::vector<double>& host_data);
void moveToHost(double*& device_data, std::vector<double>& host_data);



void moveToDevice(double**& device_data, std::vector<double*>& host_data);
void deleteFromDevice(double**& device_data);


void deleteFromDevice(int*& device_data);
void moveToDevice(int*& device_data, std::vector<int>& host_data);


void intializeOnDevice(double* device_data, const size_t nb_points);
void intializeOnDevice2D(std::vector< double* >& device_data, const size_t nb_rows);


}


#endif
