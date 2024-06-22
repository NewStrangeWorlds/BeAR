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
#include <cmath>
#include <stdio.h>

#include "../forward_model/modules/stellar_contamination.h"

#include "error_check.h"
#include "planck_function.h"
#include "data_management_kernels.h"
#include "../additional/physical_const.h"


namespace bear{


__global__ void addStellarContamination(
  double* spectrum,
  double* spectrum_phot,
  double* spectrum_fac,
  double* spectrum_spot,
  const double frac_fac,
  const double frac_spot,
  int nb_wavenumbers)
{
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_wavenumbers; tid += blockDim.x * gridDim.x)
  {
    double flux_spot = 0;
    if (frac_spot > 0) flux_spot = spectrum_spot[tid];

    double flux_fac = 0;
    if (frac_fac > 0) flux_fac = spectrum_fac[tid];

    const double stellar_activity = 1.0 
      - frac_spot * (1.0 - flux_spot/spectrum_phot[tid])
      - frac_fac * (1.0 - flux_fac/spectrum_phot[tid]);

    spectrum[tid] /= stellar_activity;

    // if (tid == 0) printf("%d %f %f %f %1.6e %1.6e %1.6e\n",
    //   tid, stellar_activity, frac_spot, frac_fac, spectrum_spot[tid], spectrum_fac[tid], spectrum_phot[tid]);
  }
}



__host__ void StellarContamination::modifySpectrumGPU(
  const std::vector<double>& parameter,
  Atmosphere* atmosphere,
  double* spectrum_gpu)
{
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  std::vector<double> stellar_param(
    parameter.begin(), 
    parameter.begin()+nb_stellar_model_param);

  const double temperature_phot = parameter[0];

  const double temperature_fac = temperature_phot + parameter[nb_stellar_model_param];
  const double temperature_spot = temperature_phot - parameter[nb_stellar_model_param+1];
  const double fraction_fac = parameter[nb_stellar_model_param+2];
  const double fraction_spot = parameter[nb_stellar_model_param+3];
  
  if (spectrum_phot_gpu == nullptr)
    allocateOnDevice(spectrum_phot_gpu, nb_spectral_points);

  stellar_model->calcFluxGPU(stellar_param, spectrum_phot_gpu);

  if (fraction_fac > 0)
  {
    if (spectrum_fac_gpu == nullptr)
      allocateOnDevice(spectrum_fac_gpu, nb_spectral_points);

    stellar_param[0] = temperature_fac;
    stellar_model->calcFluxGPU(stellar_param, spectrum_fac_gpu);
  }

  if (fraction_spot > 0)
  {
    if (spectrum_spot_gpu == nullptr)
      allocateOnDevice(spectrum_spot_gpu, nb_spectral_points);

    stellar_param[0] = temperature_spot;
    stellar_model->calcFluxGPU(stellar_param, spectrum_spot_gpu);
  }

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;


  addStellarContamination<<<blocks,threads>>>(
    spectrum_gpu,
    spectrum_phot_gpu,
    spectrum_fac_gpu,
    spectrum_spot_gpu,
    fraction_fac,
    fraction_spot,
    nb_spectral_points);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


}
 
