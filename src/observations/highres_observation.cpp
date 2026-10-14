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


#include "highres_observation.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <omp.h>

#include <Eigen/Dense>

#include "highres_model_transform.h"
#include "../additional/physical_const.h"
#include "../additional/exceptions.h"
#include "../CUDA_kernels/highres_loglike_kernels.h"


namespace bear {


HighResObservation::~HighResObservation()
{
  freeDeviceMemory();
}


void HighResObservation::init(const std::string& file_path)
{
  loadDataFile(file_path);

  // Determine overall wavelength range
  wavelength_min = 1e30;
  wavelength_max = 0;

  for (const auto& order : spectral_orders)
  {
    if (!order.wavelengths.empty())
    {
      wavelength_min = std::min(wavelength_min, order.wavelengths.front());
      wavelength_max = std::max(wavelength_max, order.wavelengths.back());
    }
  }

  std::cout << "High-res observation loaded: " << observation_name << "\n"
            << "  Orders: " << nb_orders
            << ", Exposures: " << nb_exposures
            << ", R = " << resolving_power << "\n"
            << "  Wavelength range: " << wavelength_min
            << " - " << wavelength_max << " nm\n";

  // Default barycentric velocities to zero if not provided in file.
  if (barycentric_velocities.empty())
    barycentric_velocities.assign(nb_exposures, 0.0);

  // Enable exposure (velocity) blurring only when both the orbital period and the
  // per-exposure integration times are available.
  exposure_blurring = (orbital_period > 0.0 && exposure_times.size() == nb_exposures);

  if (orbital_period > 0.0 && exposure_times.size() != nb_exposures
                           && !exposure_times.empty())
    std::cout << "  WARNING: #exposure_times count (" << exposure_times.size()
              << ") != nb_exposures (" << nb_exposures
              << "); exposure blurring disabled\n";

  if (exposure_blurring)
    std::cout << "  Exposure blurring enabled\n";

  loadModelScale(file_path);

  // Re-injection puts the model into detector units so that (I-P) acts on it the
  // way it acted on the data.  Without model filtering it has no meaning.
  if ((reinject_model || has_model_scale) && (!has_filtering || !filter_model))
  {
    std::string error_message =
      "#reinject_model 1 and #model_scale require a #filtering_basis section and "
      "#filter_model 1 in file: " + file_path + "\n";
    throw InvalidInput(std::string("HighResObservation::init"), error_message);
  }

  if (reinject_model && has_model_scale)
    std::cout << "  WARNING: #model_scale section ignored, #reinject_model 1 "
              << "uses S = P*raw_flux instead\n";

  if (has_filtering)
    initFiltering();

  if (has_flux_uncertainties)
    std::cout << "  Per-pixel flux uncertainties loaded\n";
}


// Precompute (I - P) projection matrices and apply filtering to observed data.
// Gibson et al. 2022, Section 3.3, Eq. 7:
//   P = U (ΛU)† Λ  where Λ = diag(1/σ̄_exp)
//   M_filtered = (I - P) M
// A column of ones is automatically appended to U to ensure mean removal.
void HighResObservation::initFiltering()
{
  const size_t ne = nb_exposures;
  const size_t nb = nb_basis_vectors;
  const size_t nb_aug = nb + 1;  // augmented: original basis + column of ones

  projection_matrices.resize(nb_orders);

  for (size_t ord = 0; ord < nb_orders; ++ord)
  {
    // Build augmented U matrix: [U_raw | 1]
    Eigen::MatrixXd U(ne, nb_aug);

    for (size_t e = 0; e < ne; ++e)
    {
      for (size_t b = 0; b < nb; ++b)
        U(e, b) = basis_vectors_raw[ord][e * nb + b];

      U(e, nb) = 1.0;  // column of ones for mean removal
    }

    // Build Λ = diag(1/σ̄) or identity if no uncertainties
    Eigen::VectorXd lambda(ne);

    if (!uncertainties[ord].empty())
    {
      for (size_t e = 0; e < ne; ++e)
        lambda(e) = 1.0 / uncertainties[ord][e];
    }
    else
    {
      lambda.setOnes();
    }

    // P = U * inv(U^T Λ² U) * U^T * Λ²
    Eigen::MatrixXd L2 = lambda.array().square().matrix().asDiagonal();
    Eigen::MatrixXd UtL2U = U.transpose() * L2 * U;
    Eigen::MatrixXd P = U * UtL2U.inverse() * U.transpose() * L2;
    Eigen::MatrixXd IminusP = Eigen::MatrixXd::Identity(ne, ne) - P;

    // Store row-major
    projection_matrices[ord].resize(ne * ne);

    for (size_t i = 0; i < ne; ++i)
      for (size_t j = 0; j < ne; ++j)
        projection_matrices[ord][i * ne + j] = IminusP(i, j);
  }

  // Offsets of each order in the flattened [order_offset * ne + exp * N + pixel]
  // layout used by the re-injection matrix.
  std::vector<size_t> ord_offsets(nb_orders);
  size_t total_pix = 0;

  for (size_t ord = 0; ord < nb_orders; ++ord)
  {
    ord_offsets[ord] = total_pix;
    total_pix += spectral_orders[ord].nb_pixels;
  }

  if (reinject_model)
    model_scale_host.assign(total_pix * ne, 0.0f);

  // Apply (I-P) to observed flux and precompute filtered statistics
  filtered_flux.resize(nb_orders);
  filtered_data_mean.resize(nb_orders * ne);
  filtered_data_sf2.resize(nb_orders * ne);

  for (size_t ord = 0; ord < nb_orders; ++ord)
  {
    const auto& order = spectral_orders[ord];
    const size_t N = order.nb_pixels;
    const auto& IminusP = projection_matrices[ord];

    // Apply (I-P) to flux: filtered_flux[ord][exp][pix] = sum_j (I-P)[exp][j] * flux[j][pix]
    filtered_flux[ord].resize(ne, std::vector<double>(N, 0.0));

    for (size_t e = 0; e < ne; ++e)
    {
      for (size_t j = 0; j < ne; ++j)
      {
        const double w = IminusP[e * ne + j];

        if (std::fabs(w) < 1e-15) continue;

        for (size_t p = 0; p < N; ++p)
          filtered_flux[ord][e][p] += w * order.flux[j][p];
      }
    }

    // CHIMERA-style re-injection background S = D - (I-P) D = P D, the part of the
    // data the filter captures.  Taken before the spectral detrending below, so
    // that (I-P) S = 0 and a model that is constant in time cancels exactly.
    if (reinject_model)
    {
      const size_t off = ord_offsets[ord];

      for (size_t e = 0; e < ne; ++e)
        for (size_t p = 0; p < N; ++p)
          model_scale_host[off * ne + e * N + p] = static_cast<float>(
            order.flux[e][p] - filtered_flux[ord][e][p]);
    }

    // For filter_model=true: quadratic-detrend each exposure to remove the
    // broad CIA/blackbody continuum curvature before per-exposure cross-correlation.
    // The filtered model is detrended the same way at every likelihood call.
    // For filter_model=false (data-only mode): no detrending.  The sums there run
    // over all exposures, and the zero temporal mean per pixel guaranteed by the
    // column-of-ones augmentation is what cancels a static continuum.
    if (filter_model)
    {
      const double p_mid = 0.5 * static_cast<double>(N - 1);

      for (size_t e = 0; e < ne; ++e)
      {
        double Sy   = 0;
        double Sxy  = 0;
        double Sx2y = 0;

        for (size_t p = 0; p < N; ++p)
        {
          const double x = static_cast<double>(p) - p_mid;
          const double v = filtered_flux[ord][e][p];
          Sy   += v;
          Sxy  += x * v;
          Sx2y += x * x * v;
        }

        double a, b, c;
        quadraticFitFromSums(static_cast<double>(N), Sy, Sxy, Sx2y, a, b, c);

        for (size_t p = 0; p < N; ++p)
        {
          const double x = static_cast<double>(p) - p_mid;
          filtered_flux[ord][e][p] -= a + b * x + c * x * x;
        }
      }
    }

    // Compute filtered data mean and sf2
    for (size_t e = 0; e < ne; ++e)
    {
      double sum = 0;
      for (size_t p = 0; p < N; ++p)
        sum += filtered_flux[ord][e][p];

      double mean = sum / static_cast<double>(N);
      filtered_data_mean[ord * ne + e] = mean;

      double sf2 = 0;
      for (size_t p = 0; p < N; ++p)
      {
        double d = filtered_flux[ord][e][p] - mean;
        sf2 += d * d;
      }

      filtered_data_sf2[ord * ne + e] = sf2 / static_cast<double>(N);
    }
  }

  if (reinject_model)
  {
    has_model_scale = true;
    std::cout << "  Re-injection enabled: model_scale = P*raw_flux (SVD background)\n";
  }

  std::cout << "  Filtering initialized: projection matrices computed for "
            << nb_orders << " orders\n";
}


double HighResObservation::kpRef() const { return kp_ref; }
double HighResObservation::vsysRef() const { return vsys_ref; }


// Diagnostic export. Runs only the interpolate+filter stage of the filtered
// likelihood path, so the result is the model exactly as the likelihood builds
// it -- including the CHIMERA re-injection scaling when the observation was
// loaded with #reinject_model 1. No likelihood is evaluated. The only state
// touched is model_filtered_dev, which is scratch workspace overwritten at the
// start of every likelihood call anyway.
std::vector<float> HighResObservation::modelMatrixGPU(
  const float* broadened_spectrum_gpu,
  const double* model_wavelengths_gpu,
  size_t nb_model_points,
  double Kp, double Vsys, double dphi,
  bool apply_projection,
  const float* stellar_spectrum_gpu) const
{
  // Same convention as computeLogLikelihoodGPU: the kernels receive total
  // velocities, not the retrieved offsets.
  Kp   += kp_ref;
  Vsys += vsys_ref;

  const float* exposure_blur_coeff_arg =
    exposure_blurring ? exposure_blur_coeff_dev : nullptr;

  // With the projection, show the model exactly as the likelihood treats it:
  // (I-P) M, or the temporal centring of the data-only mode.
  const int projection = !apply_projection ? highres_projection_none
    : (filter_model ? highres_projection_filter : highres_projection_center_time);

  launchHighResInterpFilterOnly(
    broadened_spectrum_gpu,
    model_wavelengths_gpu,
    static_cast<int>(nb_model_points),
    all_wavelengths_dev,
    order_offsets_dev,
    order_nb_pixels_dev,
    orbital_phases_dev,
    barycentric_velocities_dev,
    exposure_blur_coeff_arg,
    projection_matrices_dev,
    model_filtered_dev,
    static_cast<int>(nb_orders),
    static_cast<int>(nb_exposures),
    max_pixels_per_order,
    Kp,
    Vsys,
    dphi,
    has_model_scale ? model_scale_dev : nullptr,
    projection,
    use_phase_function,
    stellar_spectrum_gpu);

  std::vector<float> host(total_pixels * nb_exposures, 0.0f);

  //moveToHost takes T*& , which a const method cannot bind to directly;
  //the pointer itself is not modified.
  float* device_ptr = model_filtered_dev;
  moveToHost(device_ptr, host);

  return host;
}


// Precompute data-only weighted sums for Gibson Eq. 4 likelihood.
// Uses filtered flux when filtering is active, but always original uncertainties.
void HighResObservation::precomputeGibsonStatistics()
{
  const size_t ne = nb_exposures;

  gibson_S1.resize(nb_orders * ne);
  gibson_Sf.resize(nb_orders * ne);
  gibson_Sff.resize(nb_orders * ne);

  for (size_t ord = 0; ord < nb_orders; ++ord)
  {
    const auto& order = spectral_orders[ord];
    const size_t N = order.nb_pixels;

    for (size_t e = 0; e < ne; ++e)
    {
      const auto& sigma = order.flux_uncertainties[e];
      const auto& f = has_filtering ? filtered_flux[ord][e] : order.flux[e];

      double s1 = 0, sf = 0, sff = 0;

      for (size_t p = 0; p < N; ++p)
      {
        const double inv_sigma2 = 1.0 / (sigma[p] * sigma[p]);
        s1  += inv_sigma2;
        sf  += f[p] * inv_sigma2;
        sff += f[p] * f[p] * inv_sigma2;
      }

      gibson_S1[ord * ne + e] = s1;
      gibson_Sf[ord * ne + e] = sf;
      gibson_Sff[ord * ne + e] = sff;
    }
  }

  std::cout << "  Gibson Eq. 4 statistics precomputed\n";
}


void HighResObservation::applyProjection(
  size_t ord,
  std::vector<std::vector<double>>& matrix) const
{
  const size_t ne = nb_exposures;
  const size_t N = matrix[0].size();
  const auto& IminusP = projection_matrices[ord];

  std::vector<std::vector<double>> result(ne, std::vector<double>(N, 0.0));

  for (size_t e = 0; e < ne; ++e)
  {
    for (size_t j = 0; j < ne; ++j)
    {
      const double w = IminusP[e * ne + j];

      if (std::fabs(w) < 1e-15) continue;

      for (size_t p = 0; p < N; ++p)
        result[e][p] += w * matrix[j][p];
    }
  }

  matrix = std::move(result);
}


// Point interpolation of the model at the Doppler-shifted wavelengths of an
// order.  model_wavelengths is DESCENDING (µm).  For each pixel the bracket
// index is the first grid point at or below the shifted wavelength, the same
// as the binary search of the GPU kernels; pixels outside the grid stay 0.
void HighResObservation::interpolateModelOntoOrder(
  const SpectralOrder& order,
  const std::vector<double>& broadened_spectrum,
  const std::vector<double>& model_wavelengths,
  double doppler_inv,
  std::vector<double>& model_on_order) const
{
  const size_t N = order.nb_pixels;
  const size_t n_model = model_wavelengths.size();

  model_on_order.assign(N, 0.0);

  if (N == 0 || n_model < 2)
    return;

  // Seed the running bracket with a binary search at the first pixel
  const double wl_first = order.wavelengths[0] * 1e-3 * doppler_inv;
  auto it = std::lower_bound(
    model_wavelengths.begin(), model_wavelengths.end(),
    wl_first, std::greater<double>());
  size_t model_idx = std::distance(model_wavelengths.begin(), it);

  for (size_t p = 0; p < N; ++p)
  {
    const double wl_shifted = order.wavelengths[p] * 1e-3 * doppler_inv;

    // model_wavelengths[model_idx - 1] > wl_shifted >= model_wavelengths[model_idx]
    while (model_idx > 0 && model_wavelengths[model_idx - 1] <= wl_shifted)
      --model_idx;
    while (model_idx < n_model && model_wavelengths[model_idx] > wl_shifted)
      ++model_idx;

    if (model_idx < 1 || model_idx >= n_model)
      continue;

    const double w1 = model_wavelengths[model_idx - 1];
    const double w2 = model_wavelengths[model_idx];
    const double t = (wl_shifted - w1) / (w2 - w1);

    model_on_order[p] = (1.0 - t) * broadened_spectrum[model_idx - 1]
                      + t * broadened_spectrum[model_idx];
  }
}


// Exposure (boxcar) blurring via a precomputed cumulative integral.
// The mean of the model over a wavelength interval [wl_lo, wl_hi] equals
// (I(wl_hi) - I(wl_lo)) / (wl_hi - wl_lo), where I is the cumulative integral.
// Cost is O(1) per pixel, independent of the box width.  In the limit delta_v -> 0
// this reduces exactly to the linear point interpolation of interpolateModelOntoOrder.
void HighResObservation::interpolateModelOntoOrderBlurred(
  const SpectralOrder& order,
  const std::vector<double>& broadened_spectrum,
  const std::vector<double>& model_cumulative_integral,
  const std::vector<double>& model_wavelengths,
  double v_rad,
  double delta_v,
  std::vector<double>& model_on_order) const
{
  const double c_kms = constants::light_c * 1e-5;  // cm/s -> km/s
  const size_t N = order.nb_pixels;
  const size_t n_model = model_wavelengths.size();

  model_on_order.assign(N, 0.0);

  // Box edges in velocity map to Doppler factors.  A higher velocity is a larger
  // redshift (smaller rest wavelength), so inv_lo (from v_rad + half) < inv_hi.
  const double half   = 0.5 * std::abs(delta_v);
  const double inv_lo = 1.0 / (1.0 + (v_rad + half) / c_kms);
  const double inv_hi = 1.0 / (1.0 + (v_rad - half) / c_kms);

  // Cumulative integral I(wl) evaluated with a running descending bracket.
  // model_wavelengths is sorted in DESCENDING order (µm).  Returns NaN when wl
  // falls outside the model grid so the caller can skip that pixel.
  auto integral_at = [&](double wl, size_t& idx) -> double {
    if (wl >= model_wavelengths[0] || wl <= model_wavelengths[n_model - 1])
      return std::nan("");

    // model_wavelengths[idx] <= wl < model_wavelengths[idx-1]
    while (idx > 1 && model_wavelengths[idx - 1] <= wl) --idx;
    while (idx < n_model - 1 && model_wavelengths[idx] > wl) ++idx;

    const double w1 = model_wavelengths[idx - 1];  // larger
    const double w2 = model_wavelengths[idx];      // smaller
    const double f1 = broadened_spectrum[idx - 1];
    const double f2 = broadened_spectrum[idx];

    const double d    = w1 - wl;         // >= 0
    const double f_wl = f1 + (f2 - f1) * (d / (w1 - w2));

    return model_cumulative_integral[idx - 1] + 0.5 * (f1 + f_wl) * d;
  };

  // Seed the running brackets with a binary search at the first pixel's box edges.
  // The order spans only a small sub-range of the model grid, so starting from the
  // array end would force an O(n_model) linear scan on the first pixel (per order,
  // per exposure).  Subsequent pixels then advance the brackets in O(1).
  auto seed_bracket = [&](double wl) -> size_t {
    auto it = std::lower_bound(model_wavelengths.begin(), model_wavelengths.end(),
                               wl, std::greater<double>());
    size_t idx = static_cast<size_t>(std::distance(model_wavelengths.begin(), it));
    if (idx < 1) idx = 1;
    if (idx > n_model - 1) idx = n_model - 1;
    return idx;
  };

  const double wl0 = order.wavelengths[0] * 1e-3;
  size_t idx_lo = seed_bracket(wl0 * inv_lo);
  size_t idx_hi = seed_bracket(wl0 * inv_hi);

  for (size_t p = 0; p < N; ++p)
  {
    const double wl_base = order.wavelengths[p] * 1e-3;  // nm -> µm
    const double wl_lo   = wl_base * inv_lo;             // smaller rest wl
    const double wl_hi   = wl_base * inv_hi;             // larger rest wl

    const double i_lo = integral_at(wl_lo, idx_lo);
    const double i_hi = integral_at(wl_hi, idx_hi);

    if (std::isnan(i_lo) || std::isnan(i_hi))
      continue;  // box extends beyond the model grid: leave 0

    // integral_at returns I(wl) = ∫_wl^{λmax} S dλ (integrated *above* wl), so the
    // box integral ∫_{wl_lo}^{wl_hi} = I(wl_lo) - I(wl_hi) = i_lo - i_hi.
    model_on_order[p] = (i_lo - i_hi) / (wl_hi - wl_lo);
  }
}


void HighResObservation::initDeviceMemory()
{
  // The filtered GPU path keeps one pixel column of the model (all exposures) in
  // a fixed-size local array.
  if (has_filtering && nb_exposures > static_cast<size_t>(highres_max_exposures_gpu))
  {
    std::string error_message =
      "The filtered GPU likelihood supports at most "
      + std::to_string(highres_max_exposures_gpu) + " exposures per observation, but "
      + observation_name + " has " + std::to_string(nb_exposures)
      + ". Split the sequence into several observation files or run on the CPU.\n";
    throw InvalidInput(std::string("HighResObservation::initDeviceMemory"), error_message);
  }

  // Compute total pixels across all orders and max pixels per order
  total_pixels = 0;
  max_pixels_per_order = 0;

  for (const auto& order : spectral_orders)
  {
    total_pixels += order.nb_pixels;
    if (static_cast<int>(order.nb_pixels) > max_pixels_per_order)
      max_pixels_per_order = static_cast<int>(order.nb_pixels);
  }

  // Build flattened wavelength array and per-order offsets
  // Observed wavelengths kept in double on the device: as float (~500 nm) they carry
  // ~18 m/s of positional error, which shifts the Doppler-interpolation point and, on
  // sharp lines, changes the model far more than float rounding — the dominant source
  // of the CPU/GPU high-res likelihood mismatch.
  std::vector<double> all_wavelengths;
  std::vector<int> offsets(nb_orders);
  std::vector<int> nb_pixels_vec(nb_orders);

  all_wavelengths.reserve(total_pixels);
  size_t offset = 0;

  for (size_t ord = 0; ord < nb_orders; ++ord)
  {
    const auto& order = spectral_orders[ord];
    offsets[ord] = static_cast<int>(offset);
    nb_pixels_vec[ord] = static_cast<int>(order.nb_pixels);

    for (size_t p = 0; p < order.nb_pixels; ++p)
      all_wavelengths.push_back(order.wavelengths[p]);

    offset += order.nb_pixels;
  }

  moveToDevice(all_wavelengths_dev, all_wavelengths);
  moveToDevice(order_offsets_dev, offsets);
  moveToDevice(order_nb_pixels_dev, nb_pixels_vec);

  // Build flattened flux array.
  // When filtering is enabled, use filtered data; otherwise use raw data.
  // Layout: all_flux[order_offset * nb_exposures + exp * N + pixel]
  std::vector<float> all_flux(total_pixels * nb_exposures);

  for (size_t ord = 0; ord < nb_orders; ++ord)
  {
    const auto& order = spectral_orders[ord];
    const size_t N = order.nb_pixels;
    const size_t ord_offset = offsets[ord];

    for (size_t e = 0; e < nb_exposures; ++e)
      for (size_t p = 0; p < N; ++p)
      {
        const double val = has_filtering ? filtered_flux[ord][e][p]
                                         : order.flux[e][p];
        all_flux[ord_offset * nb_exposures + e * N + p] =
          static_cast<float>(val);
      }
  }

  moveToDevice(all_flux_dev, all_flux);

  // Orbital phases and barycentric velocities, in double as on the CPU
  std::vector<double> phases(orbital_phases.begin(), orbital_phases.end());
  moveToDevice(orbital_phases_dev, phases);

  std::vector<double> vbary(barycentric_velocities.begin(), barycentric_velocities.end());
  moveToDevice(barycentric_velocities_dev, vbary);

  // Upload the per-exposure exposure-blur coefficient (2*pi/P)*t_exp.  Kept null
  // when blurring is disabled so the kernels fall back to point interpolation.
  if (exposure_blurring)
  {
    std::vector<float> blur_coeff(nb_exposures);
    for (size_t e = 0; e < nb_exposures; ++e)
      blur_coeff[e] = static_cast<float>(
        (2.0 * constants::pi / orbital_period) * exposure_times[e]);
    moveToDevice(exposure_blur_coeff_dev, blur_coeff);
  }

  // Precompute per-order-per-exposure data mean and sf2 (variance) from the
  // float values the kernels actually read.  The kernels form R_xf and R_ff from
  // those values; taking sf2 from the same numbers keeps R_xf^2 <= N sf2 R_ff
  // (Cauchy-Schwarz), so the likelihood argument cannot turn negative through
  // float/double rounding when the noise is tiny compared with the signal.
  std::vector<float> data_mean_h(nb_orders * nb_exposures);
  std::vector<double> data_sf2_h(nb_orders * nb_exposures);

  for (size_t ord = 0; ord < nb_orders; ++ord)
  {
    const size_t N = spectral_orders[ord].nb_pixels;
    const size_t ord_offset = offsets[ord];

    for (size_t e = 0; e < nb_exposures; ++e)
    {
      const float* row = &all_flux[ord_offset * nb_exposures + e * N];

      double sum = 0;
      for (size_t p = 0; p < N; ++p)
        sum += static_cast<double>(row[p]);

      const float mean = static_cast<float>(sum / static_cast<double>(N));
      data_mean_h[ord * nb_exposures + e] = mean;

      double sf2 = 0;
      for (size_t p = 0; p < N; ++p)
      {
        const double d = static_cast<double>(row[p]) - static_cast<double>(mean);
        sf2 += d * d;
      }

      data_sf2_h[ord * nb_exposures + e] = sf2 / static_cast<double>(N);
    }
  }

  moveToDevice(data_mean_dev, data_mean_h);
  moveToDevice(data_sf2_dev, data_sf2_h);

  // Upload projection matrices and allocate filtered model workspace (GPU)
  if (has_filtering)
  {
    std::vector<float> proj_flat(nb_orders * nb_exposures * nb_exposures);

    for (size_t ord = 0; ord < nb_orders; ++ord)
      for (size_t i = 0; i < nb_exposures * nb_exposures; ++i)
        proj_flat[ord * nb_exposures * nb_exposures + i] =
          static_cast<float>(projection_matrices[ord][i]);

    moveToDevice(projection_matrices_dev, proj_flat);

    // Workspace for filtered model: [total_pixels * nb_exposures]
    allocateOnDevice(model_filtered_dev, total_pixels * nb_exposures);
  }

  // Upload model scale matrix for re-injection (CHIMERA-style)
  if (has_model_scale)
    moveToDevice(model_scale_dev, model_scale_host);

  // Upload flux uncertainties and Gibson statistics for Gibson likelihood
  if (likelihood_mode == HighResLikelihoodMode::gibson)
  {
    // Flatten flux uncertainties: same layout as all_flux
    std::vector<float> all_unc(total_pixels * nb_exposures);

    for (size_t ord = 0; ord < nb_orders; ++ord)
    {
      const auto& order = spectral_orders[ord];
      const size_t N = order.nb_pixels;
      const size_t ord_offset = offsets[ord];

      for (size_t e = 0; e < nb_exposures; ++e)
        for (size_t p = 0; p < N; ++p)
          all_unc[ord_offset * nb_exposures + e * N + p] =
            static_cast<float>(order.flux_uncertainties[e][p]);
    }

    // Data sums from the same float data and uncertainties the kernels read, for
    // the same reason as sf2 above: chi^2(alpha) then stays >= 0 up to rounding.
    std::vector<double> S1_h(nb_orders * nb_exposures);
    std::vector<double> Sf_h(nb_orders * nb_exposures);
    std::vector<double> Sff_h(nb_orders * nb_exposures);

    for (size_t ord = 0; ord < nb_orders; ++ord)
    {
      const size_t N = spectral_orders[ord].nb_pixels;
      const size_t ord_offset = offsets[ord];

      for (size_t e = 0; e < nb_exposures; ++e)
      {
        const float* f = &all_flux[ord_offset * nb_exposures + e * N];
        const float* s = &all_unc[ord_offset * nb_exposures + e * N];

        double s1 = 0, sf = 0, sff = 0;

        for (size_t p = 0; p < N; ++p)
        {
          const double inv_s2 = 1.0 / (static_cast<double>(s[p]) * static_cast<double>(s[p]));
          const double fp = static_cast<double>(f[p]);
          s1  += inv_s2;
          sf  += fp * inv_s2;
          sff += fp * fp * inv_s2;
        }

        S1_h[ord * nb_exposures + e]  = s1;
        Sf_h[ord * nb_exposures + e]  = sf;
        Sff_h[ord * nb_exposures + e] = sff;
      }
    }

    moveToDevice(flux_uncertainties_dev, all_unc);
    moveToDevice(gibson_S1_dev, S1_h);
    moveToDevice(gibson_Sf_dev, Sf_h);
    moveToDevice(gibson_Sff_dev, Sff_h);
  }

  std::cout << "  GPU memory initialized: " << total_pixels << " total pixels, "
            << max_pixels_per_order << " max per order\n";
}


void HighResObservation::freeDeviceMemory()
{
  if (all_wavelengths_dev != nullptr) deleteFromDevice(all_wavelengths_dev);
  if (all_flux_dev != nullptr) deleteFromDevice(all_flux_dev);
  if (order_offsets_dev != nullptr) deleteFromDevice(order_offsets_dev);
  if (order_nb_pixels_dev != nullptr) deleteFromDevice(order_nb_pixels_dev);
  if (orbital_phases_dev != nullptr) deleteFromDevice(orbital_phases_dev);
  if (barycentric_velocities_dev != nullptr) deleteFromDevice(barycentric_velocities_dev);
  if (exposure_blur_coeff_dev != nullptr) deleteFromDevice(exposure_blur_coeff_dev);
  if (data_mean_dev != nullptr) deleteFromDevice(data_mean_dev);
  if (data_sf2_dev != nullptr) deleteFromDevice(data_sf2_dev);
  if (projection_matrices_dev != nullptr) deleteFromDevice(projection_matrices_dev);
  if (model_filtered_dev != nullptr) deleteFromDevice(model_filtered_dev);
  if (model_scale_dev != nullptr) deleteFromDevice(model_scale_dev);
  if (flux_uncertainties_dev != nullptr) deleteFromDevice(flux_uncertainties_dev);
  if (gibson_S1_dev != nullptr) deleteFromDevice(gibson_S1_dev);
  if (gibson_Sf_dev != nullptr) deleteFromDevice(gibson_Sf_dev);
  if (gibson_Sff_dev != nullptr) deleteFromDevice(gibson_Sff_dev);
}


// Linear interpolation on the descending model grid, 0 outside it (the same
// bracket rule as the GPU kernels).
static double interpolateDescending(
  const std::vector<double>& grid,
  const double* values,
  const size_t n,
  const double wl)
{
  auto it = std::lower_bound(grid.begin(), grid.begin() + n, wl, std::greater<double>());
  const size_t idx = static_cast<size_t>(std::distance(grid.begin(), it));

  if (idx < 1 || idx >= n)
    return 0.0;

  const double t = (wl - grid[idx - 1]) / (grid[idx] - grid[idx - 1]);

  return (1.0 - t) * values[idx - 1] + t * values[idx];
}


// CPU log-likelihood.  Mirrors the GPU path step by step (see
// highres_loglike_kernels.cu): model value at the Doppler target (box-averaged
// over the exposure smear when enabled), emission stellar-frame correction or
// transmission depth conversion, phase weight, re-injection scale, then
//   filtered, filter_model = true: (I-P), quadratic detrend, per-exposure logL
//   filtered, filter_model = false: temporal centring, one logL per order
//   no filtering:                  mean subtraction, per-exposure logL
// with the Gibson form using inverse-variance weighted means instead.
double HighResObservation::computeLogLikelihood(
  const std::vector<double>& broadened_spectrum,
  const std::vector<double>& model_wavelengths,
  double Kp, double Vsys, double dphi, double alpha,
  const double* stellar_spectrum,
  size_t nb_stellar_points) const
{
  const double c_kms = constants::light_c * 1e-5;  // cm/s -> km/s
  const size_t ne = nb_exposures;
  const size_t n_model = model_wavelengths.size();

  const bool emission  = stellar_spectrum != nullptr && nb_stellar_points > 0;
  const size_t n_stellar = std::min(n_model, nb_stellar_points);
  const bool filtered  = has_filtering && filter_model;
  const bool data_only = has_filtering && !filter_model;
  const bool gibson    = likelihood_mode == HighResLikelihoodMode::gibson;

  // Brogi & Line forms: a negative alpha selects the maximum-likelihood alpha
  const double alpha_bl =
    (likelihood_mode == HighResLikelihoodMode::free_alpha) ? alpha : -1.0;

  // Exposure blurring: precompute the cumulative trapezoidal integral of the model
  // once (shared by all orders and exposures).  model_wavelengths is descending, so
  // (lambda[k-1] - lambda[k]) > 0 and the integral increases with index.  Kept in
  // double for precision over the full high-res grid.
  std::vector<double> model_cumint;
  if (exposure_blurring)
  {
    model_cumint.assign(n_model, 0.0);
    for (size_t k = 1; k < n_model; ++k)
      model_cumint[k] = model_cumint[k-1]
        + 0.5 * (broadened_spectrum[k-1] + broadened_spectrum[k])
              * (model_wavelengths[k-1] - model_wavelengths[k]);
  }

  // Per-exposure geometry, identical to exposureGeometry() in the GPU kernels
  const double kp_total   = kp_ref + Kp;
  const double vsys_total = vsys_ref + Vsys;

  std::vector<double> v_rad(ne), inv_dop(ne), delta_v(ne, 0.0), phase_w(ne, 1.0);

  for (size_t e = 0; e < ne; ++e)
  {
    const double phase = orbital_phases[e] + dphi;

    v_rad[e]   = kp_total * std::sin(2.0 * constants::pi * phase)
               + vsys_total + barycentric_velocities[e];
    inv_dop[e] = 1.0 / (1.0 + v_rad[e] / c_kms);

    // Full width of the intra-exposure velocity smear (km/s).  Only the planet's
    // orbital term accelerates over an exposure; Vsys and v_bary are ~constant.
    if (exposure_blurring)
      delta_v[e] = kp_total * std::cos(2.0 * constants::pi * phase)
                 * (2.0 * constants::pi / orbital_period) * exposure_times[e];

    if (use_phase_function)
      phase_w[e] = highResPhaseWeight(phase);
  }

  // Order offsets in the flattened re-injection matrix
  std::vector<size_t> ord_offsets(nb_orders, 0);
  for (size_t o = 1; o < nb_orders; ++o)
    ord_offsets[o] = ord_offsets[o - 1] + spectral_orders[o - 1].nb_pixels;

  double total_log_like = 0;

  // Process each order (parallelized over orders)
  #pragma omp parallel for reduction(+:total_log_like) schedule(dynamic, 1)
  for (size_t ord = 0; ord < nb_orders; ++ord)
  {
    const auto& order = spectral_orders[ord];
    const size_t N = order.nb_pixels;
    const double dN = static_cast<double>(N);

    // Fs at the rest wavelength of each pixel (emission), the same for all exposures
    std::vector<double> fs_rest;

    if (emission)
    {
      fs_rest.resize(N);
      for (size_t p = 0; p < N; ++p)
        fs_rest[p] = interpolateDescending(
          model_wavelengths, stellar_spectrum, n_stellar, order.wavelengths[p] * 1e-3);
    }

    // 1) Model at all exposures
    std::vector<std::vector<double>> model_matrix(ne, std::vector<double>(N, 0.0));

    for (size_t exp = 0; exp < ne; ++exp)
    {
      auto& row = model_matrix[exp];

      // Box-average at the Doppler target when the smear is non-negligible (near
      // conjunction); otherwise point-interpolate.
      if (exposure_blurring && std::abs(delta_v[exp]) > exposure_blur_min_kms)
        interpolateModelOntoOrderBlurred(order, broadened_spectrum, model_cumint,
                                         model_wavelengths, v_rad[exp], delta_v[exp], row);
      else
        interpolateModelOntoOrder(order, broadened_spectrum, model_wavelengths,
                                  inv_dop[exp], row);

      for (size_t p = 0; p < N; ++p)
      {
        double m = row[p];

        if (!emission)
        {
          // transmission: transit depth (ppm) -> relative in-transit flux change
          m *= transit_depth_to_flux;
        }
        else if (fs_rest[p] > 0.0)
        {
          // Fs(λ')/Fs(λ_rest) undoes the Doppler shift that interpolating Fp/Fs at
          // λ' gives the stellar spectrum: only Fp moves, Fs stays in the stellar
          // rest frame of the data.
          m *= interpolateDescending(model_wavelengths, stellar_spectrum, n_stellar,
                                     order.wavelengths[p] * 1e-3 * inv_dop[exp])
             / fs_rest[p];
        }

        m *= phase_w[exp];

        // Re-injection (Line et al. 2021): model in detector units, m * S
        if (has_model_scale)
          m *= static_cast<double>(model_scale_host[ord_offsets[ord] * ne + exp * N + p]);

        row[p] = m;
      }
    }

    // 2) Temporal treatment of the model
    if (filtered)
    {
      applyProjection(ord, model_matrix);
    }
    else if (data_only)
    {
      // N_PCA = 0 model filtering: remove the temporal mean of each pixel
      for (size_t p = 0; p < N; ++p)
      {
        double mean = 0.0;
        for (size_t exp = 0; exp < ne; ++exp)
          mean += model_matrix[exp][p];

        mean /= static_cast<double>(ne);

        for (size_t exp = 0; exp < ne; ++exp)
          model_matrix[exp][p] -= mean;
      }
    }

    // 3a) Data-only mode: sums over all exposures of the order, one logL per order.
    //     The (I-P)-filtered data have zero temporal mean per pixel and the model is
    //     centred in time, so a continuum that does not move with the planet cancels.
    if (data_only)
    {
      double rxf = 0, rff = 0, sff = 0;
      double sm = 0, sfm = 0, smm = 0;

      for (size_t exp = 0; exp < ne; ++exp)
      {
        const auto& d_row = filtered_flux[ord][exp];
        const auto& m_row = model_matrix[exp];

        for (size_t p = 0; p < N; ++p)
        {
          const double d = d_row[p];
          const double m = m_row[p];

          rxf += d * m;
          rff += m * m;
          sff += d * d;

          if (gibson)
          {
            const double s = order.flux_uncertainties[exp][p];
            const double inv_s2 = 1.0 / (s * s);
            sm  += m * inv_s2;
            sfm += d * m * inv_s2;
            smm += m * m * inv_s2;
          }
        }
      }

      const double dN_total = dN * static_cast<double>(ne);

      if (gibson)
      {
        double S1 = 0, Sf = 0, Sff = 0;

        for (size_t exp = 0; exp < ne; ++exp)
        {
          S1  += gibson_S1[ord * ne + exp];
          Sf  += gibson_Sf[ord * ne + exp];
          Sff += gibson_Sff[ord * ne + exp];
        }

        total_log_like += gibsonLogLike(S1, Sf, Sff, sm, sfm, smm, dN_total, alpha);
      }
      else
      {
        total_log_like += brogiLineLogLike(rxf, rff, sff / dN_total, dN_total, alpha_bl);
      }

      continue;
    }

    // 3b) Per-exposure likelihood
    const double p_mid = 0.5 * static_cast<double>(N - 1);

    for (size_t exp = 0; exp < ne; ++exp)
    {
      const auto& flux_ref = has_filtering ? filtered_flux[ord][exp] : order.flux[exp];
      auto& m_row = model_matrix[exp];

      // Model detrending: a quadratic in pixel index when the model is filtered (the
      // filtered data were detrended the same way at initialisation), otherwise the
      // mean.  The Gibson form removes weighted means through its sums instead.
      if (filtered)
      {
        double Sy = 0, Sxy = 0, Sx2y = 0;

        for (size_t p = 0; p < N; ++p)
        {
          const double x = static_cast<double>(p) - p_mid;
          Sy   += m_row[p];
          Sxy  += x * m_row[p];
          Sx2y += x * x * m_row[p];
        }

        double a, b, c;
        quadraticFitFromSums(dN, Sy, Sxy, Sx2y, a, b, c);

        for (size_t p = 0; p < N; ++p)
        {
          const double x = static_cast<double>(p) - p_mid;
          m_row[p] -= a + b * x + c * x * x;
        }
      }
      else if (!gibson)
      {
        double mean = 0;
        for (size_t p = 0; p < N; ++p)
          mean += m_row[p];

        mean /= dN;

        for (size_t p = 0; p < N; ++p)
          m_row[p] -= mean;
      }

      if (gibson)
      {
        const auto& sigma = order.flux_uncertainties[exp];

        double Sm = 0, Sfm = 0, Smm = 0;

        for (size_t p = 0; p < N; ++p)
        {
          const double m = m_row[p];
          const double inv_sigma2 = 1.0 / (sigma[p] * sigma[p]);
          Sm  += m * inv_sigma2;
          Sfm += flux_ref[p] * m * inv_sigma2;
          Smm += m * m * inv_sigma2;
        }

        const size_t idx = ord * ne + exp;

        total_log_like += gibsonLogLike(
          gibson_S1[idx], gibson_Sf[idx], gibson_Sff[idx], Sm, Sfm, Smm, dN, alpha);
      }
      else
      {
        double data_mean, sf2;

        if (has_filtering)
        {
          data_mean = filtered_data_mean[ord * ne + exp];
          sf2 = filtered_data_sf2[ord * ne + exp];
        }
        else
        {
          double sum = 0;
          for (size_t p = 0; p < N; ++p)
            sum += flux_ref[p];
          data_mean = sum / dN;

          sf2 = 0;
          for (size_t p = 0; p < N; ++p)
          {
            const double d = flux_ref[p] - data_mean;
            sf2 += d * d;
          }
          sf2 /= dN;
        }

        double rxf = 0, rff = 0;

        for (size_t p = 0; p < N; ++p)
        {
          const double d = flux_ref[p] - data_mean;
          rxf += d * m_row[p];
          rff += m_row[p] * m_row[p];
        }

        if (rff > 0.0)
          total_log_like += brogiLineLogLike(rxf, rff, sf2, dN, alpha_bl);
      }
    }
  }

  return total_log_like;
}


void HighResObservation::computeLogLikelihoodGPU(
  const float* broadened_spectrum_gpu,
  const double* model_wavelengths_gpu,
  size_t nb_model_points,
  double Kp, double Vsys, double dphi, double alpha,
  double* d_log_like_dev,
  const float* stellar_spectrum_gpu) const
{
  // Fold in reference offsets so CUDA kernels receive the total velocity
  Kp   += kp_ref;
  Vsys += vsys_ref;

  const bool gibson = likelihood_mode == HighResLikelihoodMode::gibson;
  const int likelihood_form = gibson ? highres_form_gibson
    : (likelihood_mode == HighResLikelihoodMode::free_alpha ? highres_form_free_alpha
                                                            : highres_form_marginalized_alpha);

  // Exposure blurring: the kernels box-average the model at the Doppler target using
  // exposure_blur_coeff_dev (null when disabled -> point interpolation).
  const float* exposure_blur_coeff_arg =
    exposure_blurring ? exposure_blur_coeff_dev : nullptr;

  if (has_filtering)
  {
    launchHighResLogLikeFiltered(
      broadened_spectrum_gpu,
      model_wavelengths_gpu,
      static_cast<int>(nb_model_points),
      all_wavelengths_dev,
      all_flux_dev,
      order_offsets_dev,
      order_nb_pixels_dev,
      data_mean_dev,
      data_sf2_dev,
      orbital_phases_dev,
      barycentric_velocities_dev,
      exposure_blur_coeff_arg,
      projection_matrices_dev,
      model_filtered_dev,
      static_cast<int>(nb_orders),
      static_cast<int>(nb_exposures),
      max_pixels_per_order,
      Kp, Vsys, dphi,
      alpha,
      likelihood_form,
      flux_uncertainties_dev,
      gibson_S1_dev,
      gibson_Sf_dev,
      gibson_Sff_dev,
      has_model_scale ? model_scale_dev : nullptr,
      filter_model,
      use_phase_function,
      stellar_spectrum_gpu,
      d_log_like_dev);
  }
  else if (gibson)
  {
    launchHighResLogLikeGibson(
      broadened_spectrum_gpu,
      model_wavelengths_gpu,
      static_cast<int>(nb_model_points),
      all_wavelengths_dev,
      all_flux_dev,
      order_offsets_dev,
      order_nb_pixels_dev,
      flux_uncertainties_dev,
      gibson_S1_dev,
      gibson_Sf_dev,
      gibson_Sff_dev,
      orbital_phases_dev,
      barycentric_velocities_dev,
      exposure_blur_coeff_arg,
      static_cast<int>(nb_orders),
      static_cast<int>(nb_exposures),
      max_pixels_per_order,
      Kp, Vsys, dphi,
      alpha,
      d_log_like_dev,
      stellar_spectrum_gpu,
      use_phase_function);
  }
  else
  {
    launchHighResLogLike(
      broadened_spectrum_gpu,
      model_wavelengths_gpu,
      static_cast<int>(nb_model_points),
      all_wavelengths_dev,
      all_flux_dev,
      order_offsets_dev,
      order_nb_pixels_dev,
      data_mean_dev,
      data_sf2_dev,
      orbital_phases_dev,
      barycentric_velocities_dev,
      exposure_blur_coeff_arg,
      static_cast<int>(nb_orders),
      static_cast<int>(nb_exposures),
      max_pixels_per_order,
      Kp, Vsys, dphi,
      likelihood_form == highres_form_free_alpha ? alpha : -1.0,
      d_log_like_dev,
      stellar_spectrum_gpu,
      use_phase_function);
  }
}


}
