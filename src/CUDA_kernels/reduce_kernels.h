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


#ifndef _reduce_kernels_h
#define _reduce_kernels_h


#pragma once

#include <vector>

namespace bear{



__inline__ __device__
double warpReduceSum(double val)
{
  for (int offset = warpSize/2; offset > 0; offset /= 2)
    //val += __shfl_down(val, offset);
    val += __shfl_down_sync(0xFFFFFFFF, val, offset); 


  return val;
}


__inline__ __device__
double blockReduceSum(double val)
{

  static __shared__ double shared[32]; // Shared mem for 32 partial sums
  int lane = threadIdx.x % warpSize;
  int wid = threadIdx.x / warpSize;

  // Each warp performs partial reduction
  val = warpReduceSum(val);


  // Write reduced value to shared memory
  if (lane==0) shared[wid]=val;


   // Wait for all partial reductions
  __syncthreads();


  //read from shared memory only if that warp existed
  val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

  if (wid==0) val = warpReduceSum(val); //Final reduce within first warp


  return val;
}


}


#endif

