// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "data_management_kernels.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <new>
#include <algorithm>

#include "error_check.h"



namespace bear{

template <typename T>
__host__ void deleteFromDevice(T*& device_data)
{
  if (device_data != nullptr)
    gpuErrchk(cudaFree(device_data));

  device_data = nullptr;
}


template <typename T>
__host__ void moveToHost(T*& device_data, std::vector<T>& host_data)
{
  const int bytes = host_data.size()*sizeof(T);

  gpuErrchk(cudaMemcpy(host_data.data(), device_data, bytes, cudaMemcpyDeviceToHost));
}



template <typename T>
__host__ void moveToHostAndDelete(T*& device_data, std::vector<T>& host_data)
{
  const int bytes = host_data.size()*sizeof(T);

  gpuErrchk(cudaMemcpy(host_data.data(), device_data, bytes, cudaMemcpyDeviceToHost));

  deleteFromDevice(device_data);
}



template <typename T>
__host__ void allocateOnDevice(T*& device_data, const size_t nb_values)
{
  const int bytes = nb_values*sizeof(T);

  gpuErrchk(cudaMalloc((void**)&device_data, bytes));
}



template <typename T>
__host__ void moveToDevice(T*& device_data, std::vector<T>& host_data, const bool alloc_memory)
{
  if (alloc_memory)
    allocateOnDevice(device_data, host_data.size());

  moveToDevice(device_data, host_data);
}



template <typename T>
__host__ void moveToDevice(T*& device_data, std::vector<T>& host_data)
{
  const int bytes = host_data.size()*sizeof(T);

  if (device_data == nullptr)
    allocateOnDevice(device_data, host_data.size());

  gpuErrchk(cudaMemcpy(device_data, &host_data[0], bytes, cudaMemcpyHostToDevice));
}


__host__
void moveToDevice(float*& device_data, std::vector<double>& host_data)
{
  if (device_data == nullptr)
    allocateOnDevice(device_data, host_data.size());

  const int bytes = host_data.size()*sizeof(float);

  std::vector<float> host_data_f(host_data.size());

  std::transform(host_data.begin(), host_data.end(), host_data_f.begin(),
               [](double x) { return static_cast<float>(x); });

  gpuErrchk(cudaMemcpy(device_data, &host_data_f[0], bytes, cudaMemcpyHostToDevice));
}


//sets all entries of a 1D double array on the GPU to 0
//device_data holds the pointer to the GPU array
//nb_points is length of the array
template <typename T>
__host__ void initializeOnDevice(T*& device_data, const size_t nb_points)
{
  if (device_data != nullptr)
    gpuErrchk(cudaMemset(device_data, 0, nb_points*sizeof(T)));
}


template void moveToHost<double>(double*&, std::vector<double>&);
template void moveToHost<float>(float*&, std::vector<float>&);
template void moveToHost<int>(int*&, std::vector<int>&);

template void moveToHostAndDelete<double>(double*&, std::vector<double>&);
template void moveToHostAndDelete<float>(float*&, std::vector<float>&);
template void moveToHostAndDelete<int>(int*&, std::vector<int>&);

template void deleteFromDevice<double>(double*&);
template void deleteFromDevice<int>(int*&);
template void deleteFromDevice<float>(float*&);
template void deleteFromDevice<float*>(float**&);

template void initializeOnDevice<double>(double*&, const size_t);
template void initializeOnDevice<int>(int*&, const size_t);
template void initializeOnDevice<float>(float*&, const size_t);

template void moveToDevice<double>(double*&, std::vector<double>&);
template void moveToDevice<int>(int*&, std::vector<int>&);
template void moveToDevice<float>(float*&, std::vector<float>&);

template void moveToDevice<double>(double*&, std::vector<double>&, const bool);
template void moveToDevice<int>(int*&, std::vector<int>&, const bool);
template void moveToDevice<float>(float*&, std::vector<float>&, const bool);

template void allocateOnDevice<double>(double*&, size_t);
template void allocateOnDevice<int>(int*&, size_t);
template void allocateOnDevice<float>(float*&, size_t);
template void allocateOnDevice<float*>(float**&, size_t);

template <typename T>
__host__ void copyOnDevice(T* dst, const T* src, size_t nb_elements)
{
  gpuErrchk(cudaMemcpy(dst, src, nb_elements * sizeof(T), cudaMemcpyDeviceToDevice));
}

template void copyOnDevice<float>(float*, const float*, size_t);
template void copyOnDevice<double>(double*, const double*, size_t);
}
