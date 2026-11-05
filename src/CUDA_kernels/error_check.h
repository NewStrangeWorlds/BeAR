// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _error_check_h
#define _error_check_h

#include <cstdio>
#include <cstdlib>


namespace bear{


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


#ifdef DEBUG_CUDA
  #define CUDA_CHECK_AFTER_KERNEL() \
    do { gpuErrchk(cudaDeviceSynchronize()); } while(0)
#else
  #define CUDA_CHECK_AFTER_KERNEL() \
    gpuErrchk(cudaPeekAtLastError())
#endif


}


#endif
