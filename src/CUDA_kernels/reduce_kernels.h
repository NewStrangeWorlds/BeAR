// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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

  // Make sure no warp is still reading `shared` from a previous
  // blockReduceSum call before it is overwritten below
  __syncthreads();

  // Write reduced value to shared memory
  if (lane==0) shared[wid]=val;


   // Wait for all partial reductions
  __syncthreads();


  //read from shared memory only if that warp existed
  val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

  if (wid==0) val = warpReduceSum(val); //Final reduce within first warp


  return val;
}



__inline__ __device__
float warpReduceSum(float val)
{
  for (int offset = warpSize/2; offset > 0; offset /= 2)
    //val += __shfl_down(val, offset);
    val += __shfl_down_sync(0xFFFFFFFF, val, offset); 


  return val;
}


__inline__ __device__
float blockReduceSum(float val)
{
  static __shared__ float shared[32]; // Shared mem for 32 partial sums
  int lane = threadIdx.x % warpSize;
  int wid = threadIdx.x / warpSize;

  // Each warp performs partial reduction
  val = warpReduceSum(val);

  // Make sure no warp is still reading `shared` from a previous
  // blockReduceSum call before it is overwritten below
  __syncthreads();

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

