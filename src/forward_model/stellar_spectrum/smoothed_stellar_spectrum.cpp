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


#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

#include "smoothed_stellar_spectrum.h"
#include "../../CUDA_kernels/data_management_kernels.h"


namespace bear {


SmoothedStellarSpectrum::SmoothedStellarSpectrum(
  std::unique_ptr<StellarSpectrumModel> inner,
  double sigma_px,
  bool use_gpu)
  : inner_(std::move(inner))
  , sigma_px_(sigma_px)
  , use_gpu_(use_gpu)
{
  nb_parameters = inner_->nbParameters();
  std::cout << "  Stellar spectrum smoothing enabled: sigma = "
            << sigma_px_ << " px\n";
}


SmoothedStellarSpectrum::~SmoothedStellarSpectrum()
{
  if (gpu_cache_ != nullptr)
    deleteFromDevice(gpu_cache_);
}


// O(N) box filter pass using prefix sums.
// At array edges the window is clipped and the mean is over the available elements.
void SmoothedStellarSpectrum::boxFilterPass(
  const std::vector<double>& in,
  std::vector<double>& out,
  int half_w)
{
  const int n = static_cast<int>(in.size());
  out.resize(n);

  std::vector<double> prefix(n + 1, 0.0);
  for (int i = 0; i < n; ++i)
    prefix[i + 1] = prefix[i] + in[i];

  for (int i = 0; i < n; ++i)
  {
    const int lo = std::max(0, i - half_w);
    const int hi = std::min(n - 1, i + half_w);
    out[i] = (prefix[hi + 1] - prefix[lo]) / static_cast<double>(hi - lo + 1);
  }
}


// 3-pass box-filter approximation of a Gaussian — O(N) cost, errors < 0.1% for sigma > 5 px.
// Box half-width derived from: 3 * (k*(k+2)/3) = sigma^2 => k = sqrt(sigma^2+1) - 1.
std::vector<double> SmoothedStellarSpectrum::smooth(const std::vector<double>& in)
{
  const int half_w = static_cast<int>(
    std::round(std::sqrt(sigma_px_ * sigma_px_ + 1.0) - 1.0));

  std::vector<double> tmp1, tmp2, out;
  boxFilterPass(in,   tmp1, half_w);
  boxFilterPass(tmp1, tmp2, half_w);
  boxFilterPass(tmp2, out,  half_w);
  return out;
}


std::vector<double> SmoothedStellarSpectrum::calcFlux(
  const std::vector<double>& parameter)
{
  if (cache_valid_)
    return cpu_cache_;

  std::vector<double> raw = inner_->calcFlux(parameter);
  std::vector<double> smoothed = smooth(raw);

  if (nb_parameters == 0)
  {
    cpu_cache_   = smoothed;
    cache_valid_ = true;
  }

  return smoothed;
}


void SmoothedStellarSpectrum::calcFluxGPU(
  const std::vector<double>& parameter,
  float* spectrum_gpu)
{
  // Fast path: cached GPU spectrum is available.
  if (cache_valid_ && gpu_cache_ != nullptr)
  {
    copyOnDevice(spectrum_gpu, gpu_cache_, cpu_cache_.size());
    return;
  }

  // Compute smoothed spectrum on CPU (uses cpu_cache_ if already filled).
  std::vector<double> smoothed = calcFlux(parameter);

  // Upload to the caller-provided device buffer.
  // moveToDevice(float*&, vector<double>&) converts and copies without reallocating
  // because spectrum_gpu is already non-null (allocated by the caller).
  float* out_ptr = spectrum_gpu;
  moveToDevice(out_ptr, smoothed);

  // Cache the GPU result for parameter-free models so subsequent calls are instant.
  if (nb_parameters == 0 && gpu_cache_ == nullptr)
    moveToDevice(gpu_cache_, smoothed);  // gpu_cache_ is null → allocates and copies
}


}
