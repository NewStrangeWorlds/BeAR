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
#include <omp.h>

#include <Eigen/Dense>

#include "../additional/physical_const.h"
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

  std::cout << "  Filtering initialized: projection matrices computed for "
            << nb_orders << " orders\n";
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

  // Binary search for starting position
  const double wl_first = order.wavelengths[0] * 1e-3 * doppler_inv;
  auto it = std::lower_bound(
    model_wavelengths.begin(), model_wavelengths.end(),
    wl_first, std::greater<double>());
  size_t model_idx = std::distance(model_wavelengths.begin(), it);

  if (model_idx == 0) model_idx = 1;

  for (size_t p = 0; p < N; ++p)
  {
    const double wl_shifted = order.wavelengths[p] * 1e-3 * doppler_inv;

    while (model_idx > 1 && model_wavelengths[model_idx - 1] < wl_shifted)
      --model_idx;

    if (model_idx < 1 || model_idx >= n_model)
      continue;

    double w1 = model_wavelengths[model_idx - 1];
    double w2 = model_wavelengths[model_idx];
    double t = (wl_shifted - w1) / (w2 - w1);
    // Convert from transit depth in ppm to normalised flux (1 - depth)
    // so that model and data are in the same units and alpha ~ +1
    model_on_order[p] = 1.0 - ((1.0 - t) * broadened_spectrum[model_idx - 1]
                              + t * broadened_spectrum[model_idx]) * 1e-6;
  }
}


void HighResObservation::initDeviceMemory()
{
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
  std::vector<float> all_wavelengths;
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
      all_wavelengths.push_back(static_cast<float>(order.wavelengths[p]));

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

  // Upload orbital phases
  std::vector<float> phases_float(orbital_phases.begin(), orbital_phases.end());
  moveToDevice(orbital_phases_dev, phases_float);

  // Precompute per-order-per-exposure data mean and sf2 (variance)
  std::vector<float> data_mean_h(nb_orders * nb_exposures);
  std::vector<double> data_sf2_h(nb_orders * nb_exposures);

  for (size_t ord = 0; ord < nb_orders; ++ord)
  {
    const auto& order = spectral_orders[ord];
    const size_t N = order.nb_pixels;

    for (size_t e = 0; e < nb_exposures; ++e)
    {
      if (has_filtering)
      {
        data_mean_h[ord * nb_exposures + e] =
          static_cast<float>(filtered_data_mean[ord * nb_exposures + e]);
        data_sf2_h[ord * nb_exposures + e] =
          filtered_data_sf2[ord * nb_exposures + e];
      }
      else
      {
        double sum = 0;
        for (size_t p = 0; p < N; ++p)
          sum += order.flux[e][p];

        double mean = sum / static_cast<double>(N);
        data_mean_h[ord * nb_exposures + e] = static_cast<float>(mean);

        double sf2 = 0;
        for (size_t p = 0; p < N; ++p)
        {
          double d = order.flux[e][p] - mean;
          sf2 += d * d;
        }
        sf2 /= static_cast<double>(N);
        data_sf2_h[ord * nb_exposures + e] = sf2;
      }
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

    moveToDevice(flux_uncertainties_dev, all_unc);
    moveToDevice(gibson_S1_dev, gibson_S1);
    moveToDevice(gibson_Sf_dev, gibson_Sf);
    moveToDevice(gibson_Sff_dev, gibson_Sff);
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
  if (data_mean_dev != nullptr) deleteFromDevice(data_mean_dev);
  if (data_sf2_dev != nullptr) deleteFromDevice(data_sf2_dev);
  if (projection_matrices_dev != nullptr) deleteFromDevice(projection_matrices_dev);
  if (model_filtered_dev != nullptr) deleteFromDevice(model_filtered_dev);
  if (flux_uncertainties_dev != nullptr) deleteFromDevice(flux_uncertainties_dev);
  if (gibson_S1_dev != nullptr) deleteFromDevice(gibson_S1_dev);
  if (gibson_Sf_dev != nullptr) deleteFromDevice(gibson_Sf_dev);
  if (gibson_Sff_dev != nullptr) deleteFromDevice(gibson_Sff_dev);
}


// CPU log-likelihood with unified order-first structure.
// When filtering is enabled, the model is interpolated at all exposures,
// then filtered via (I-P) before computing the cross-correlation likelihood.
double HighResObservation::computeLogLikelihood(
  const std::vector<double>& broadened_spectrum,
  const std::vector<double>& model_wavelengths,
  double Kp, double Vsys, double alpha) const
{
  const double c_kms = constants::light_c * 1e-5;  // cm/s -> km/s
  double total_log_like = 0;

  // Precompute inverse Doppler factors for all exposures
  std::vector<double> doppler_inv(nb_exposures);
  for (size_t exp = 0; exp < nb_exposures; ++exp)
  {
    const double phase = orbital_phases[exp];
    const double v_rad = Kp * std::sin(2.0 * constants::pi * phase) + Vsys;
    doppler_inv[exp] = 1.0 / (1.0 + v_rad / c_kms);
  }

  // Process each order (parallelized over orders)
  #pragma omp parallel for reduction(+:total_log_like) schedule(dynamic, 1)
  for (size_t ord = 0; ord < nb_orders; ++ord)
  {
    const auto& order = spectral_orders[ord];
    const size_t N = order.nb_pixels;

    // 1) Interpolate model at all exposures
    std::vector<std::vector<double>> model_matrix(
      nb_exposures, std::vector<double>(N, 0.0));

    for (size_t exp = 0; exp < nb_exposures; ++exp)
    {
      interpolateModelOntoOrder(
        order, broadened_spectrum, model_wavelengths,
        doppler_inv[exp], model_matrix[exp]);
    }

    // 2) Apply filtering if enabled: model_matrix = (I-P) @ model_matrix
    if (has_filtering)
      applyProjection(ord, model_matrix);

    // 3) Cross-correlation likelihood for each exposure
    for (size_t exp = 0; exp < nb_exposures; ++exp)
    {
      // Select appropriate data source
      const auto& flux_ref = has_filtering ? filtered_flux[ord][exp]
                                           : order.flux[exp];

      double data_mean, sf2;

      if (has_filtering)
      {
        data_mean = filtered_data_mean[ord * nb_exposures + exp];
        sf2 = filtered_data_sf2[ord * nb_exposures + exp];
      }
      else
      {
        double sum = 0;
        for (size_t p = 0; p < N; ++p)
          sum += flux_ref[p];
        data_mean = sum / static_cast<double>(N);

        sf2 = 0;
        for (size_t p = 0; p < N; ++p)
        {
          double d = flux_ref[p] - data_mean;
          sf2 += d * d;
        }
        sf2 /= static_cast<double>(N);
      }

      // Compute model mean
      double model_mean = 0;
      for (size_t p = 0; p < N; ++p)
        model_mean += model_matrix[exp][p];
      model_mean /= static_cast<double>(N);

      if (likelihood_mode == HighResLikelihoodMode::gibson)
      {
        // Gibson et al. 2022 Eq. 4: per-pixel uncertainty weighting, beta marginalized
        // chi2 = (Sff - Sf^2/S1) + alpha^2*(Smm - Sm^2/S1) - 2*alpha*(Sfm - Sf*Sm/S1)
        // ln L = -N/2 * ln(chi2 / N)
        const auto& sigma = spectral_orders[ord].flux_uncertainties[exp];
        const double S1  = gibson_S1[ord * nb_exposures + exp];
        const double Sf  = gibson_Sf[ord * nb_exposures + exp];
        const double Sff = gibson_Sff[ord * nb_exposures + exp];

        double Sm = 0, Sfm = 0, Smm = 0;

        for (size_t p = 0; p < N; ++p)
        {
          const double m = model_matrix[exp][p];
          const double inv_sigma2 = 1.0 / (sigma[p] * sigma[p]);
          Sm  += m * inv_sigma2;
          Sfm += flux_ref[p] * m * inv_sigma2;
          Smm += m * m * inv_sigma2;
        }

        const double dN = static_cast<double>(N);
        const double chi2 = (Sff - Sf * Sf / S1)
                           + alpha * alpha * (Smm - Sm * Sm / S1)
                           - 2.0 * alpha * (Sfm - Sf * Sm / S1);

        if (chi2 > 0)
          total_log_like += -0.5 * dN * std::log(chi2 / dN);
        else
          total_log_like += -1e30;
      }
      else
      {
        // Brogi & Line 2019 cross-correlation likelihood
        double R_xf = 0, R_ff = 0;

        for (size_t p = 0; p < N; ++p)
        {
          const double d = flux_ref[p] - data_mean;
          const double m = model_matrix[exp][p] - model_mean;
          R_xf += d * m;
          R_ff += m * m;
        }

        if (R_ff > 0.0)
        {
          const double dN = static_cast<double>(N);
          double arg;

          if (likelihood_mode == HighResLikelihoodMode::marginalized_alpha)
          {
            // Normalized: ln L = -N/2 * ln((sf2 - R_xf^2/(N*R_ff)) / sf2)
            // = -N/2 * ln(1 - r^2), r = cross-correlation coefficient.
            // Dividing by sf2 removes the data-scale constant -N/2*ln(sf2),
            // which is ~1e8 for PCA-filtered data and breaks nested sampling.
            arg = (sf2 - (R_xf * R_xf) / (dN * R_ff)) / sf2;
          }
          else
          {
            // Explicit alpha: normalized by sg2 for the same reason.
            const double sg2 = R_ff / dN;
            const double R = R_xf / dN;
            arg = (sf2 + alpha * alpha * sg2 - 2.0 * alpha * R) / sf2;
          }

          if (arg > 0)
            total_log_like += -0.5 * N * std::log(arg);
          else
            total_log_like += -1e30;
        }
      }
    }
  }

  return total_log_like;
}


void HighResObservation::computeLogLikelihoodGPU(
  const float* broadened_spectrum_gpu,
  const double* model_wavelengths_gpu,
  size_t nb_model_points,
  double Kp, double Vsys, double alpha,
  double* d_log_like_dev) const
{
  // alpha_gpu: positive = explicit alpha; negative = marginalize
  const float alpha_gpu = (likelihood_mode == HighResLikelihoodMode::marginalized_alpha)
    ? -1.0f : static_cast<float>(alpha);

  if (likelihood_mode == HighResLikelihoodMode::gibson)
  {
    if (has_filtering)
    {
      launchHighResLogLikeFilteredGibson(
        broadened_spectrum_gpu,
        model_wavelengths_gpu,
        static_cast<int>(nb_model_points),
        all_wavelengths_dev,
        all_flux_dev,
        order_offsets_dev,
        order_nb_pixels_dev,
        orbital_phases_dev,
        projection_matrices_dev,
        model_filtered_dev,
        flux_uncertainties_dev,
        gibson_S1_dev,
        gibson_Sf_dev,
        gibson_Sff_dev,
        static_cast<int>(nb_orders),
        static_cast<int>(nb_exposures),
        max_pixels_per_order,
        static_cast<float>(Kp),
        static_cast<float>(Vsys),
        static_cast<float>(alpha),
        d_log_like_dev);
    }
    else
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
        static_cast<int>(nb_orders),
        static_cast<int>(nb_exposures),
        max_pixels_per_order,
        static_cast<float>(Kp),
        static_cast<float>(Vsys),
        static_cast<float>(alpha),
        d_log_like_dev);
    }
  }
  else if (has_filtering)
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
      projection_matrices_dev,
      model_filtered_dev,
      static_cast<int>(nb_orders),
      static_cast<int>(nb_exposures),
      max_pixels_per_order,
      static_cast<float>(Kp),
      static_cast<float>(Vsys),
      alpha_gpu,
      d_log_like_dev);
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
      static_cast<int>(nb_orders),
      static_cast<int>(nb_exposures),
      max_pixels_per_order,
      static_cast<float>(Kp),
      static_cast<float>(Vsys),
      alpha_gpu,
      d_log_like_dev);
  }
}


}
