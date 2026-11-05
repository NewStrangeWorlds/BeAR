// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _smoothed_stellar_spectrum_h
#define _smoothed_stellar_spectrum_h

#include <vector>
#include <memory>

#include "stellar_spectrum.h"
#include "../../CUDA_kernels/data_management_kernels.h"


namespace bear {


// Decorator that applies Gaussian pixel-space smoothing to any StellarSpectrumModel.
// For parameter-free models (e.g. StarSpectrumFile) the smoothed result is
// computed once and cached; for parametric models (e.g. StellarSpectrumGrid)
// it is recomputed every call.
// Smoothing uses a fast 3-pass box-filter approximation of the Gaussian,
// which is O(N) regardless of sigma.
class SmoothedStellarSpectrum : public StellarSpectrumModel {
  public:
    SmoothedStellarSpectrum(
      std::unique_ptr<StellarSpectrumModel> inner,
      double sigma_px,
      bool use_gpu);
    virtual ~SmoothedStellarSpectrum();

    virtual std::vector<double> calcFlux(
      const std::vector<double>& parameter) override;
    virtual void calcFluxGPU(
      const std::vector<double>& parameter,
      float* spectrum_gpu) override;

  private:
    std::unique_ptr<StellarSpectrumModel> inner_;
    double sigma_px_;
    bool use_gpu_;

    bool cache_valid_ = false;
    std::vector<double> cpu_cache_;
    float* gpu_cache_ = nullptr;

    std::vector<double> smooth(const std::vector<double>& in);
    void boxFilterPass(
      const std::vector<double>& in,
      std::vector<double>& out,
      int half_w);
};


}

#endif
