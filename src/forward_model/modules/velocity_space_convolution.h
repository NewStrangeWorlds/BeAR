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


#ifndef _velocity_space_convolution_h
#define _velocity_space_convolution_h

#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>

#include "../../additional/physical_const.h"


namespace bear {


// Velocity steps between neighbouring points of a spectral grid,
// dv[i] = c |ln(λ_i / λ_{i+1})| in km/s (n - 1 values).
//
// The high-resolution grid is a subset of the opacity wavenumber points, so these
// steps are not constant: at R = 2.5e5 on the 0.01 cm^-1 HELIOS-k grid they vary
// by a factor of ~3 across 1.45-2.45 µm.  Broadening kernels therefore have to
// use the actual steps rather than one value for the whole grid.
inline std::vector<double> gridVelocitySteps(const std::vector<double>& wavelengths)
{
  const double c_kms = constants::light_c * 1e-5;
  const size_t n = wavelengths.size();

  std::vector<double> dv(n > 1 ? n - 1 : 0, 0.0);

  for (size_t i = 0; i + 1 < n; ++i)
    dv[i] = c_kms * std::fabs(std::log(wavelengths[i] / wavelengths[i + 1]));

  return dv;
}


// Normalised convolution in velocity space on an irregular grid:
//
//   out_i = Σ_j K(u_ij) w_j in_j / Σ_j K(u_ij) w_j,     |u_ij| <= half_width,
//
// where u_ij is the velocity offset of point j from point i (the sum of the steps
// between them, positive towards shorter wavelengths) and w_j is the trapezoidal
// weight of point j.  Normalising per output point keeps the kernel at unit area
// however coarsely it is sampled, and also at the edges of the grid.  A kernel
// narrower than the local step leaves the spectrum unchanged.
template <typename Profile>
void convolveVelocitySpace(
  const std::vector<double>& in,
  std::vector<double>& out,
  const std::vector<double>& dv,
  const double half_width,
  const Profile& profile)
{
  const long n = static_cast<long>(in.size());
  out.resize(n);

  #pragma omp parallel for schedule(static)
  for (long i = 0; i < n; ++i)
  {
    const double w_i = 0.5 * ((i > 0 ? dv[i - 1] : 0.0) + (i < n - 1 ? dv[i] : 0.0));

    double k = profile(0.0) * w_i;
    double sum = k * in[i];
    double norm = k;

    double u = 0.0;

    for (long j = i + 1; j < n; ++j)
    {
      u += dv[j - 1];
      if (u > half_width) break;

      const double w_j = 0.5 * (dv[j - 1] + (j < n - 1 ? dv[j] : 0.0));
      k = profile(u) * w_j;
      sum += k * in[j];
      norm += k;
    }

    u = 0.0;

    for (long j = i - 1; j >= 0; --j)
    {
      u += dv[j];
      if (u > half_width) break;

      const double w_j = 0.5 * ((j > 0 ? dv[j - 1] : 0.0) + dv[j]);
      k = profile(-u) * w_j;
      sum += k * in[j];
      norm += k;
    }

    out[i] = norm > 0.0 ? sum / norm : in[i];
  }
}


// Rotational broadening profile of Gray (2005, Eq. 18.14), without its
// normalisation constant (the convolution normalises per point).
struct RotationalProfile {
  double vsini = 0, c1 = 0, c2 = 0;

  RotationalProfile(const double vsini_kms, const double epsilon)
    : vsini(vsini_kms)
    , c1(2.0 * (1.0 - epsilon))
    , c2(0.5 * constants::pi * epsilon) {}

  double operator()(const double u) const {
    const double x = u / vsini;
    if (std::fabs(x) >= 1.0) return 0.0;
    const double q = 1.0 - x * x;
    return c1 * std::sqrt(q) + c2 * q;
  }
};


// Gaussian profile exp(-u^2 / 2 sigma^2), unnormalised.
struct GaussianProfile {
  double inv_2sig2 = 0;

  explicit GaussianProfile(const double sigma_kms)
    : inv_2sig2(1.0 / (2.0 * sigma_kms * sigma_kms)) {}

  double operator()(const double u) const {
    return std::exp(-u * u * inv_2sig2);
  }
};


// Profile tabulated at velocities v_min + k step, linearly interpolated,
// zero outside the table.
struct TabulatedProfile {
  const std::vector<double>& table;
  double v_min = 0, step = 1;

  TabulatedProfile(const std::vector<double>& table_, const double v_min_, const double step_)
    : table(table_), v_min(v_min_), step(step_) {}

  double operator()(const double u) const {
    const double x = (u - v_min) / step;
    if (x < 0.0) return 0.0;
    const size_t k = static_cast<size_t>(x);
    if (k + 1 >= table.size()) return 0.0;
    const double t = x - static_cast<double>(k);
    return (1.0 - t) * table[k] + t * table[k + 1];
  }
};


}


#endif
