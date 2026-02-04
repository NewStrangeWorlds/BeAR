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
    cudaFree(device_data);

  device_data = nullptr;

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


template <typename T> 
__host__ void moveToHost(T*& device_data, std::vector<T>& host_data)
{
  const int bytes = host_data.size()*sizeof(T);

  cudaMemcpy(host_data.data(), device_data, bytes, cudaMemcpyDeviceToHost);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() ); 
}



template <typename T> 
__host__ void moveToHostAndDelete(T*& device_data, std::vector<T>& host_data)
{
  const int bytes = host_data.size()*sizeof(T);

  cudaMemcpy(host_data.data(), device_data, bytes, cudaMemcpyDeviceToHost);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );

  deleteFromDevice(device_data);
}



template <typename T> 
__host__ void allocateOnDevice(T*& device_data, const size_t nb_values)
{
  const int bytes = nb_values*sizeof(T);

  auto ret = cudaMalloc((void**)&device_data, bytes);

  if (ret == cudaErrorMemoryAllocation)
  {
    std::cerr << "Error: Could not allocate memory on device\n";
    throw std::bad_alloc();
  }

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
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

  cudaMemcpy(device_data, &host_data[0], bytes, cudaMemcpyHostToDevice);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
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

  cudaMemcpy(device_data, &host_data_f[0], bytes, cudaMemcpyHostToDevice);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


//sets all entries of a 1D double array on the GPU to 0
//device_data holds the pointer to the GPU array
//nb_points is length of the array
template <typename T> 
__host__ void initializeOnDevice(T*& device_data, const size_t nb_points)
{
  if (device_data != nullptr)
    cudaMemset(device_data, 0, nb_points*sizeof(T));

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
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
}
