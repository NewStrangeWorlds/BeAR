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


#ifndef HIGHRES_MODEL_TRANSFORM_H
#define HIGHRES_MODEL_TRANSFORM_H

#include <cmath>


// Included by both the CPU likelihood (host compiler) and the CUDA kernels
// (nvcc), so that both paths apply exactly the same per-exposure weights and
// model conversions.
#ifdef __CUDACC__
  #define BEAR_HOST_DEVICE __host__ __device__
#else
  #define BEAR_HOST_DEVICE
#endif


namespace bear {


// Dayside phase function of Pelletier et al. (2025, Sect. 3.3, step 4) and
// Herman et al. (2022):  w = [0.5 (1 + cos(2 pi phi - pi))]^2,  phi = 0 at
// transit.  w = 1 at secondary eclipse, 1/4 at quadrature, 0 at transit.
// The phase passed in must already include the retrieved offset dphi.
inline BEAR_HOST_DEVICE double highResPhaseWeight(const double phase)
{
  const double pi = 3.14159265358979323846;
  const double half_pf = 0.5 * (1.0 + cos(2.0 * pi * phase - pi));

  return half_pf * half_pf;
}


// Transit depth in ppm -> relative change of the in-transit flux, -1e-6 depth.
// The constant 1 of the in-transit flux 1 - depth is dropped on purpose: every
// likelihood path removes per-exposure or per-pixel constants from the model,
// so keeping it would only cost float precision.
constexpr double transit_depth_to_flux = -1e-6;


// Below this intra-exposure velocity smear (km/s) the model is point-interpolated
// instead of box-averaged.
constexpr double exposure_blur_min_kms = 1e-3;


// Coefficients of the least-squares quadratic a + b x + c x^2 through N points at
// the centred pixel index x = p - (N-1)/2, from the sums sy = Σy, sxy = Σxy,
// sx2y = Σx²y (the odd moments of the symmetric x vanish).
inline BEAR_HOST_DEVICE void quadraticFitFromSums(
  const double N,
  const double sy,
  const double sxy,
  const double sx2y,
  double& a,
  double& b,
  double& c)
{
  const double sigma2 = N * (N * N - 1.0) / 12.0;
  const double sigma4 = N * (N * N - 1.0) * (3.0 * N * N - 7.0) / 240.0;
  const double det    = N * sigma4 - sigma2 * sigma2;

  a = (det > 0.0)    ? (sigma4 * sy - sigma2 * sx2y) / det : sy / N;
  b = (sigma2 > 0.0) ? sxy / sigma2 : 0.0;
  c = (det > 0.0)    ? (N * sx2y - sigma2 * sy) / det : 0.0;
}


// Brogi & Line (2019) log-likelihood, normalised by the data variance so that a
// model uncorrelated with the data gives zero.
//   rxf = Σ d m, rff = Σ m² (detrended model), sf2 = Σ d²/N (detrended data).
//   alpha < 0:  alpha set to its maximum-likelihood value, -N/2 ln(1 - r²)
//               (Zucker 2003; the marginalized-alpha mode)
//   alpha >= 0: -N/2 ln[1 + (alpha² s_g² - 2 alpha R)/s_f²]  (free-alpha mode)
// Returns 0 for a featureless model or data, and -1e30 for an invalid argument.
inline BEAR_HOST_DEVICE double brogiLineLogLike(
  const double rxf,
  const double rff,
  const double sf2,
  const double N,
  const double alpha)
{
  if (!(rff > 0.0) || !(sf2 > 0.0))
    return 0.0;

  const double arg = (alpha < 0.0)
    ? 1.0 - (rxf * rxf) / (N * rff * sf2)
    : 1.0 + (alpha * alpha * rff / N - 2.0 * alpha * rxf / N) / sf2;

  return (arg > 0.0) ? -0.5 * N * log(arg) : -1e30;
}


// Gibson et al. (2022) Eq. 4 with beta marginalised, normalised by chi² at
// alpha = 0.  The S_f/S_1 and S_m/S_1 terms remove the inverse-variance weighted
// means of data and model.
inline BEAR_HOST_DEVICE double gibsonLogLike(
  const double S1, const double Sf, const double Sff,
  const double Sm, const double Sfm, const double Smm,
  const double N,
  const double alpha)
{
  const double chi2_data = Sff - Sf * Sf / S1;
  const double chi2 = chi2_data
                    + alpha * alpha * (Smm - Sm * Sm / S1)
                    - 2.0 * alpha * (Sfm - Sf * Sm / S1);

  return (chi2 > 0.0 && chi2_data > 0.0) ? -0.5 * N * log(chi2 / chi2_data) : -1e30;
}


}

#endif
