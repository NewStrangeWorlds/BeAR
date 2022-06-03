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


#include "data_management_kernels.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <new>

#include "error_check.h"



namespace helios{

//a lot of the functions here can be merged using template programming
//might be added in a future version


//allocates a one-dimensional array of size nb_double_values on the GPU
__host__ void allocateOnDevice(double*& device_data, size_t nb_double_values)
{
  const int bytes = nb_double_values*sizeof(double);

  cudaMalloc((void**)&device_data, bytes);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


//deletes a 1D array of type double on the GPU and sets the pointer back to a null pointer
__host__ void deleteFromDevice(double*& device_data)
{
  if (device_data != nullptr)
    cudaFree(device_data);

  device_data = nullptr;

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}



//deletes a 1D array of type int on the GPU and sets the pointer back to a null pointer
__host__ void deleteFromDevice(int*& device_data)
{
  if (device_data != nullptr)
    cudaFree(device_data);

  device_data = nullptr;

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}




//deletes an array of pointers on the GPU and sets the pointer back to a null pointer
__host__ void deleteFromDevice(double**& device_data)
{
  if (device_data != nullptr)
    cudaFree(device_data);

  device_data = nullptr;

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}




//moves data array host_data of type double to the GPU
//returns the pointer *device_data to the data on the GPU
//if device_data already exists on the GPU it will be freed first
//if device_data has not been previously allocated, *device_data must be a null pointer!
__host__ void moveToDevice(double*& device_data, std::vector<double>& host_data)
{
  //delete the array if it has been previously allocated on the GPU
  //if (*device_data != nullptr)
    //cudaFree(*device_data);

  const int bytes = host_data.size()*sizeof(double);


  cudaMalloc((void**)&device_data, bytes);

  cudaMemcpy(device_data, &host_data[0], bytes, cudaMemcpyHostToDevice);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}



//moves data array host_data of type double to the GPU
//returns the pointer *device_data to the data on the GPU
//if device_data already exists on the GPU it will be freed first
//if device_data has not been previously allocated, *device_data must be a null pointer!
__host__ void moveToDevice(double*& device_data, std::vector<double>& host_data, const bool alloc_memory)
{
  //delete the array if it has been previously allocated on the GPU
  //if (*device_data != nullptr)
    //cudaFree(*device_data);

  const int bytes = host_data.size()*sizeof(double);

  //alloc memory if required
  if (alloc_memory) cudaMalloc((void**)&device_data, bytes);

  cudaMemcpy(device_data, &host_data[0], bytes, cudaMemcpyHostToDevice);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}



//moves data array host_data of pointers to double arrays to the GPU
//returns the pointer **device_data to the data on the GPU
//if device_data already exists on the GPU it will be freed first
//if device_data has not been previously allocated, *device_data must be a null pointer!
__host__ void moveToDevice(double**& device_data, std::vector<double*>& host_data)
{
  //delete the array if it has been previously allocated on the GPU
  //if (*device_data != nullptr)
    //cudaFree(*device_data);

  const int bytes = host_data.size()*sizeof(double*);


  cudaMalloc((void***)&device_data, bytes);

  cudaMemcpy(device_data, &host_data[0], bytes, cudaMemcpyHostToDevice);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}



//moves data array device_data to the host
__host__ void moveToHost(double*& device_data, std::vector<double>& host_data)
{
  const int bytes = host_data.size()*sizeof(double);

 
  cudaMemcpy(host_data.data(), device_data, bytes, cudaMemcpyDeviceToHost);

  
  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() ); 
}



//moves data array host_data of type int to the GPU
//returns the pointer *device_data to the data on the GPU
//if device_data already exists on the GPU it will be freed first
//if device_data has not been previously allocated, *device_data must be a null pointer!
__host__ void moveToDevice(int*& device_data, std::vector<int>& host_data)
{
  //delete the array if it has been previously allocated on the GPU
  //if (*device_data != nullptr)
    //cudaFree(*device_data);

  const int bytes = host_data.size()*sizeof(int);

  cudaMalloc((void**)&device_data, bytes);

  cudaMemcpy(device_data, &host_data[0], bytes, cudaMemcpyHostToDevice);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


//sets all entries of a 2D double array on the GPU to 0
//device_data holds the pointers to each 1D array
//nb_rows is length of each 1D array
__host__ void intializeOnDevice2D(std::vector< double* >& device_data, const size_t nb_rows)
{
  //initialize each row
  for (size_t i=0; i<device_data.size(); ++i)
  {

    if (device_data[i] != nullptr)
      cudaMemset(device_data[i], 0, nb_rows*sizeof(double));

  }


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}



//sets all entries of a 1D double array on the GPU to 0
//device_data holds the pointer to the GPU array
//nb_rows is length of the array
__host__ void intializeOnDevice(double* device_data, const size_t nb_points)
{
  if (device_data != nullptr)
    cudaMemset(device_data, 0, nb_points*sizeof(double));


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}




}

