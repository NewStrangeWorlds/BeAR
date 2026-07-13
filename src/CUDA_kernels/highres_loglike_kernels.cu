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


#include "highres_loglike_kernels.h"
#include "reduce_kernels.h"
#include "error_check.h"
#include "../additional/physical_const.h"

#include <cmath>
#include <cstdio>


namespace bear {


// Binary search in a descending array: find first index i such that arr[i] <= val
__device__ __forceinline__
int binarySearchDescending(
  const double* __restrict__ arr,
  int n,
  double val)
{
  int lo = 0, hi = n - 1;

  while (lo < hi)
  {
    int mid = (lo + hi) / 2;

    if (arr[mid] > val)
      lo = mid + 1;
    else
      hi = mid;
  }

  return lo;
}


// Interpolate the broadened model spectrum at a Doppler-shifted wavelength.
// Returns the interpolated value, or 0 if outside the model range.
// For the non-filtered (transmission) path: converts transit depth in ppm to
// normalised flux via (1 - depth * 1e-6) so that model and data share the
// same scale and alpha ~ +1 for a fixed-alpha retrieval.
__device__ __forceinline__
float interpolateModel(
  const float*  __restrict__ broadened_spectrum,
  const double* __restrict__ model_wavelengths,
  int n_model,
  double wl_shifted)
{
  int idx = binarySearchDescending(model_wavelengths, n_model, wl_shifted);

  if (idx < 1 || idx >= n_model)
    return 0.0f;

  double w1 = model_wavelengths[idx - 1];
  double w2 = model_wavelengths[idx];
  float t = (float)((wl_shifted - w1) / (w2 - w1));

  // Convert from transit depth in ppm to normalised flux (1 - depth)
  // so that model and data are in the same units and alpha ~ +1
  return 1.0f - ((1.0f - t) * broadened_spectrum[idx - 1]
               + t * broadened_spectrum[idx]) * 1e-6f;
}


// Raw linear interpolation of the model spectrum at a Doppler-shifted wavelength.
// Used in the filtered (emission) path where the spectrum is in physical flux units
// (W m⁻² cm), not in transit-depth ppm.  The Brogi & Line Pearson-r formula is
// scale-invariant, so no unit conversion is needed; only the spectral shape matters.
__device__ __forceinline__
float interpolateModelRaw(
  const float*  __restrict__ broadened_spectrum,
  const double* __restrict__ model_wavelengths,
  int n_model,
  double wl_shifted)
{
  int idx = binarySearchDescending(model_wavelengths, n_model, wl_shifted);

  if (idx < 1 || idx >= n_model)
    return 0.0f;

  double w1 = model_wavelengths[idx - 1];
  double w2 = model_wavelengths[idx];
  float t = (float)((wl_shifted - w1) / (w2 - w1));

  return (1.0f - t) * broadened_spectrum[idx - 1] + t * broadened_spectrum[idx];
}


// Like interpolateModelRaw but also returns the bracket index via out_idx,
// so the same index can be reused to interpolate a second spectrum at the
// same wavelength without a second binary search.
__device__ __forceinline__
float interpolateModelRawIdx(
  const float*  __restrict__ spectrum,
  const double* __restrict__ model_wavelengths,
  int n_model,
  double wl,
  int& out_idx)
{
  int idx = binarySearchDescending(model_wavelengths, n_model, wl);
  out_idx = idx;

  if (idx < 1 || idx >= n_model)
    return 0.0f;

  double w1 = model_wavelengths[idx - 1];
  double w2 = model_wavelengths[idx];
  float t = (float)((wl - w1) / (w2 - w1));

  return (1.0f - t) * spectrum[idx - 1] + t * spectrum[idx];
}


// Interpolate a spectrum at a pre-computed bracket index (no binary search).
// Use when the index was already found by interpolateModelRawIdx for another
// spectrum at the same wavelength.
__device__ __forceinline__
float interpolateAtIdx(
  const float*  __restrict__ spectrum,
  const double* __restrict__ model_wavelengths,
  int n_model,
  int idx,
  double wl)
{
  if (idx < 1 || idx >= n_model)
    return 0.0f;

  double w1 = model_wavelengths[idx - 1];
  double w2 = model_wavelengths[idx];
  float t = (float)((wl - w1) / (w2 - w1));

  return (1.0f - t) * spectrum[idx - 1] + t * spectrum[idx];
}


// Exposure (boxcar) blurring: mean of the spectrum over the wavelength interval
// [wl_lo, wl_hi] (wl_lo < wl_hi), computed as a direct trapezoidal integral over the
// descending model grid divided by the interval width.  Numerically identical to the
// CPU cumulative-integral formulation.  Returns false (and leaves out_mean untouched)
// when the interval falls outside the model range, so the caller can substitute 0.
// If hint >= 1 is passed (the bracket index of a nearby wavelength, e.g. the box
// centre already located by the caller), both edge brackets are found by short
// linear scans from it and no binary search is performed at all.
__device__ __forceinline__
bool boxMeanRaw(
  const float*  __restrict__ spectrum,
  const double* __restrict__ model_wavelengths,
  int n_model,
  double wl_lo,
  double wl_hi,
  double& out_mean,
  int hint = -1)
{
  // In range requires model_wavelengths[n-1] < wl_lo and wl_hi < model_wavelengths[0].
  if (wl_lo <= model_wavelengths[n_model - 1] || wl_hi >= model_wavelengths[0])
    return false;

  // Bracket for the upper edge: model_wavelengths[i] <= wl < model_wavelengths[i-1].
  // Seed from the caller's hint (scan) when available, else binary-search.
  int i_hi;
  if (hint >= 1 && hint < n_model)
  {
    i_hi = hint;
    while (i_hi > 1        && model_wavelengths[i_hi - 1] <= wl_hi) --i_hi;
    while (i_hi < n_model  && model_wavelengths[i_hi]     >  wl_hi) ++i_hi;
  }
  else
  {
    i_hi = binarySearchDescending(model_wavelengths, n_model, wl_hi);
  }

  // wl_lo < wl_hi and the grid is descending, so the lower-edge bracket is at an
  // index >= i_hi and only a few pixels away (the box is narrow).  Find it with a
  // short linear scan instead of a second binary search.
  int i_lo = i_hi;
  while (i_lo < n_model && model_wavelengths[i_lo] > wl_lo) ++i_lo;

  if (i_hi < 1 || i_lo < 1 || i_lo >= n_model)
    return false;

  // Model value at the two edges (linear interpolation within their cells).
  const double whi1 = model_wavelengths[i_hi - 1];  // larger
  const double whi2 = model_wavelengths[i_hi];       // smaller
  const double f_hi = (double)spectrum[i_hi - 1]
    + ((double)spectrum[i_hi] - (double)spectrum[i_hi - 1])
      * ((whi1 - wl_hi) / (whi1 - whi2));

  const double wlo1 = model_wavelengths[i_lo - 1];  // larger
  const double wlo2 = model_wavelengths[i_lo];       // smaller
  const double f_lo = (double)spectrum[i_lo - 1]
    + ((double)spectrum[i_lo] - (double)spectrum[i_lo - 1])
      * ((wlo1 - wl_lo) / (wlo1 - wlo2));

  double integral;

  if (i_hi == i_lo)
  {
    // Both edges share a cell.
    integral = 0.5 * (f_lo + f_hi) * (wl_hi - wl_lo);
  }
  else
  {
    // Top partial cell: wl_hi down to grid point whi2.
    integral = 0.5 * (f_hi + (double)spectrum[i_hi]) * (wl_hi - whi2);

    // Full interior cells between whi2 and wlo1.
    for (int i = i_hi; i <= i_lo - 2; ++i)
      integral += 0.5 * ((double)spectrum[i] + (double)spectrum[i + 1])
                      * (model_wavelengths[i] - model_wavelengths[i + 1]);

    // Bottom partial cell: grid point wlo1 down to wl_lo.
    integral += 0.5 * ((double)spectrum[i_lo - 1] + f_lo) * (wlo1 - wl_lo);
  }

  out_mean = integral / (wl_hi - wl_lo);
  return true;
}


// Exposure blurring, pre-convolution stage.
// Produces a per-exposure boxcar-blurred copy of the model so that the likelihood
// kernels can point-interpolate it at full speed (no per-pixel box averaging in the
// hot loop).  One thread per (exposure, model point): the boxcar half-width in
// velocity is dV(e)/2 with dV(e) = Kp cos(2*pi*(phi_e+dphi)) * coeff[e], applied as a
// wavelength interval around model_wavelengths[k] and averaged via boxMeanRaw (with
// the exact bracket hint k, so the edge scans are only a few pixels).
// blurred_model layout: [exposure * n_model + k].
__global__
void boxBlurModelKernel(
  const float*  __restrict__ broadened_spectrum,
  const double* __restrict__ model_wavelengths,
  const int                  n_model,
  const float*  __restrict__ orbital_phases,
  const float*  __restrict__ exposure_blur_coeff,
  const int                  nb_exposures,
  const float                Kp,
  const float                dphi,
  float*        __restrict__ blurred_model)
{
  const long idx = (long)blockIdx.x * blockDim.x + threadIdx.x;
  const long total = (long)nb_exposures * n_model;
  if (idx >= total) return;

  const int e = (int)(idx / n_model);
  const int k = (int)(idx % n_model);

  const float delta_v = Kp * cosf(2.0f * (float)M_PI * (orbital_phases[e] + dphi))
                      * exposure_blur_coeff[e];

  const float c_kms = (float)(constants::light_c * 1e-5);

  // Negligible smear (< ~1e-3 km/s): copy the model point unchanged.
  if (fabsf(delta_v) <= 1.0e-3f)
  {
    blurred_model[idx] = broadened_spectrum[k];
    return;
  }

  const double half  = 0.5 * (double)fabsf(delta_v) / (double)c_kms;  // fractional
  const double wl_k  = model_wavelengths[k];
  const double wl_lo = wl_k * (1.0 - half);  // smaller wavelength
  const double wl_hi = wl_k * (1.0 + half);  // larger wavelength

  double mean;
  blurred_model[idx] =
    boxMeanRaw(broadened_spectrum, model_wavelengths, n_model, wl_lo, wl_hi, mean, k)
      ? (float)mean : broadened_spectrum[k];
}


__host__
void launchBoxBlurModel(
    const float* broadened_spectrum_dev,
    const double* model_wavelengths_dev,
    int n_model,
    const float* orbital_phases_dev,
    const float* exposure_blur_coeff_dev,
    int nb_exposures,
    float Kp, float dphi,
    float* blurred_model_dev)
{
  const int threads = 256;
  const long total = (long)nb_exposures * n_model;
  const int blocks = (int)((total + threads - 1) / threads);

  boxBlurModelKernel<<<blocks, threads>>>(
    broadened_spectrum_dev,
    model_wavelengths_dev,
    n_model,
    orbital_phases_dev,
    exposure_blur_coeff_dev,
    nb_exposures,
    Kp, dphi,
    blurred_model_dev);

  CUDA_CHECK_AFTER_KERNEL();
}


// Each block handles one (order, exposure) pair.
// No shared memory required — the model is interpolated twice (cheap ALU)
// to avoid storing per-pixel values.
//
// Pass 1: Interpolate model onto order grid, accumulate model_sum → model_mean
// Pass 2: Re-interpolate, mean-subtract, compute cross-correlation sums (R_xf, R_ff)
// Thread 0: Compute log-likelihood and atomicAdd to global result
__global__
void highResLogLikeKernel(
  const float*  __restrict__ broadened_spectrum,
  const double* __restrict__ model_wavelengths,
  const int                  n_model,
  const float*  __restrict__ order_wavelengths,
  const float*  __restrict__ order_flux,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const float*  __restrict__ data_mean,
  const double* __restrict__ data_sf2,
  const float*  __restrict__ orbital_phases,
  const float*  __restrict__ v_bary,
  const float*  __restrict__ per_exposure_model,
  const int                  nb_orders,
  const int                  nb_exposures,
  const float                Kp,
  const float                Vsys,
  const float                dphi,
  const float                alpha,
  double*       __restrict__ d_log_like)
{
  const int task = blockIdx.x;
  const int exp  = task / nb_orders;
  const int ord  = task % nb_orders;
  const int tid  = threadIdx.x;

  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];

  // Pointers to this order's data
  const float* wl_order = order_wavelengths + offset;
  const float* flux = order_flux + offset * nb_exposures + exp * N;

  // Model source: the per-exposure exposure-blurred copy when active, else the
  // single broadened model (point interpolation, no in-kernel box averaging).
  const float* model_src = (per_exposure_model != nullptr)
    ? per_exposure_model + (size_t)exp * n_model
    : broadened_spectrum;

  // Doppler factor for this exposure
  const float c_kms = (float)(constants::light_c * 1e-5);
  const float phase = orbital_phases[exp];
  const float v_rad = Kp * sinf(2.0f * (float)M_PI * (phase + dphi)) + Vsys + v_bary[exp];
  const double inv_doppler = 1.0 / (1.0 + (double)v_rad / (double)c_kms);

  // ---- Pass 1: Interpolate model, accumulate model_sum ----
  float model_sum = 0.0f;

  for (int p = tid; p < N; p += blockDim.x)
  {
    const double wl_base = (double)wl_order[p] * 1e-3;
    model_sum += interpolateModel(model_src, model_wavelengths, n_model,
                                  wl_base * inv_doppler);
  }

  model_sum = blockReduceSum(model_sum);

  __shared__ float s_model_mean;

  if (tid == 0)
    s_model_mean = model_sum / (float)N;

  __syncthreads();

  // ---- Pass 2: Re-interpolate, mean-subtract, cross-correlation sums ----
  const float mmean = s_model_mean;
  const float dmean = data_mean[ord * nb_exposures + exp];

  double local_rxf = 0.0;
  double local_rff = 0.0;

  for (int p = tid; p < N; p += blockDim.x)
  {
    const double wl_base = (double)wl_order[p] * 1e-3;
    const float model_val = interpolateModel(model_src, model_wavelengths, n_model,
                                             wl_base * inv_doppler);

    double d = (double)flux[p] - (double)dmean;
    double m = (double)model_val - (double)mmean;
    local_rxf += d * m;
    local_rff += m * m;
  }

  local_rxf = blockReduceSum(local_rxf);
  local_rff = blockReduceSum(local_rff);

  // ---- Thread 0: Compute log-likelihood ----
  if (tid == 0)
  {
    const double dN = (double)N;
    const double sf2 = data_sf2[ord * nb_exposures + exp];

    if (local_rff > 0.0)
    {
      double arg;

      if (alpha < 0.0f)
      {
        // Normalized: ln L = -N/2 * ln(1 - r^2), r = cross-correlation coefficient.
        // Dividing by sf2 removes the constant -N/2*ln(sf2) that otherwise
        // reaches ~1e8 for PCA-filtered data and breaks nested sampling convergence.
        arg = (sf2 - (local_rxf * local_rxf) / (dN * local_rff)) / sf2;
      }
      else
      {
        const double a = (double)alpha;
        const double sg2 = local_rff / dN;
        const double R = local_rxf / dN;
        arg = (sf2 + a * a * sg2 - 2.0 * a * R) / sf2;
      }

      if (arg > 0.0)
        atomicAdd(d_log_like, -0.5 * dN * log(arg));
      else
        atomicAdd(d_log_like, -1e30);
    }
  }
}


__host__
void launchHighResLogLike(
    const float* broadened_spectrum_dev,
    const double* model_wavelengths_dev,
    int n_model,
    const float* order_wavelengths_dev,
    const float* order_flux_dev,
    const int* order_offsets_dev,
    const int* order_nb_pixels_dev,
    const float* data_mean_dev,
    const double* data_sf2_dev,
    const float* orbital_phases_dev,
    const float* v_bary_dev,
    const float* per_exposure_model_dev,
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    float Kp, float Vsys, float dphi,
    float alpha,
    double* d_log_like_dev,
    const float* stellar_spectrum_dev)
{
  const int threads = 256;
  const int blocks = nb_orders * nb_exposures;

  highResLogLikeKernel<<<blocks, threads>>>(
    broadened_spectrum_dev,
    model_wavelengths_dev,
    n_model,
    order_wavelengths_dev,
    order_flux_dev,
    order_offsets_dev,
    order_nb_pixels_dev,
    data_mean_dev,
    data_sf2_dev,
    orbital_phases_dev,
    v_bary_dev,
    per_exposure_model_dev,
    nb_orders,
    nb_exposures,
    Kp, Vsys, dphi,
    alpha,
    d_log_like_dev);

  CUDA_CHECK_AFTER_KERNEL();
}


// ============================================================================
// Filtered path: Gibson et al. 2022 fast model filtering
// ============================================================================


// Kernel 1: Interpolate model at all Doppler shifts, then apply (I-P) per order.
// One block per order. Each thread processes a pixel, looping over all exposures.
// (I-P) matrix is loaded into shared memory (nb_exp * nb_exp floats).
//
// Output: model_filtered[ord_offset * nb_exp + exp * N + p]
__global__
void highResInterpFilterKernel(
  const float*  __restrict__ broadened_spectrum,
  const double* __restrict__ model_wavelengths,
  const int                  n_model,
  const float*  __restrict__ order_wavelengths,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const float*  __restrict__ orbital_phases,
  const float*  __restrict__ v_bary,
  const float*  __restrict__ per_exposure_model,
  const float*  __restrict__ projection_matrices,
  float*        __restrict__ model_filtered,
  const int                  nb_orders,
  const int                  nb_exposures,
  const float                Kp,
  const float                Vsys,
  const float                dphi,
  const float*               model_scale = nullptr,
  const bool                 apply_projection = true,
  const bool                 use_phase_function = false,
  const float*               stellar_spectrum = nullptr)
{
  const int ord = blockIdx.x;
  const int tid = threadIdx.x;

  if (ord >= nb_orders) return;

  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];
  const float* wl_order = order_wavelengths + offset;

  // Load (I-P) matrix for this order into shared memory
  extern __shared__ float s_IminusP[];

  const int mat_size = nb_exposures * nb_exposures;
  const float* proj_ptr = projection_matrices + ord * mat_size;

  for (int i = tid; i < mat_size; i += blockDim.x)
    s_IminusP[i] = proj_ptr[i];

  __syncthreads();

  // Precompute Doppler factors
  const float c_kms = (float)(constants::light_c * 1e-5);

  // Process pixels in stride
  for (int p = tid; p < N; p += blockDim.x)
  {
    const double wl_um_rest = (double)wl_order[p] * 1e-3;

    // Interpolate model at all exposures for this pixel
    // Store raw (unfiltered) model values in local array.
    // For nb_exposures up to ~64, this fits in registers/local memory.
    float raw_model[128];  // max exposures supported

    // Fs(λ_rest) is constant across all exposures for this pixel — compute it once.
    // Correction factor Fs(λ·v_dop)/Fs(λ_rest) requires only one search per exposure
    // (for the Doppler-shifted wavelength), whose index is also reused for broadened_spectrum.
    int   idx_rest = 0;
    float fs_rest  = 0.0f;
    if (stellar_spectrum != nullptr)
      fs_rest = interpolateModelRawIdx(stellar_spectrum, model_wavelengths, n_model,
                                       wl_um_rest, idx_rest);

    for (int exp = 0; exp < nb_exposures; ++exp)
    {
      const float phase   = orbital_phases[exp];
      const float v_rad   = Kp * sinf(2.0f * (float)M_PI * (phase + dphi))
                          + Vsys + v_bary[exp];
      const double inv_dop   = 1.0 / (1.0 + (double)v_rad / (double)c_kms);
      const double wl_um_dop = wl_um_rest * inv_dop;

      // Model source: the per-exposure exposure-blurred copy when active, else the
      // shared broadened model.  Only the planet model is blurred; the stellar
      // template (below) always uses the original spectrum at the exposure centre.
      const float* model_src = (per_exposure_model != nullptr)
        ? per_exposure_model + (size_t)exp * n_model
        : broadened_spectrum;

      // One binary search for the Doppler-shifted wavelength; index reused for
      // both the model and the stellar spectrum at the same position.
      int idx_dop;
      raw_model[exp] = interpolateModelRawIdx(model_src, model_wavelengths, n_model,
                                              wl_um_dop, idx_dop);

      // Per-pixel correction: Fs(λ·v_dop)/Fs(λ_rest) undoes the erroneous Doppler
      // shift of the stellar template. fs_rest already computed outside the exp loop.
      if (stellar_spectrum != nullptr && fs_rest > 0.0f)
      {
        const float fs_dop = interpolateAtIdx(stellar_spectrum, model_wavelengths, n_model,
                                              idx_dop, wl_um_dop);
        raw_model[exp] *= fs_dop / fs_rest;
      }

      // Lambertian dayside phase function: 0.5*(1+cos(2*pi*phase-pi))^2
      // (Pelletier et al. 2025 §3.3 step 4, Herman et al. 2022)
      if (use_phase_function)
      {
        const float pf = 1.0f + cosf(2.0f * (float)M_PI * phase - (float)M_PI);
        raw_model[exp] *= 0.5f * pf * pf;
      }

      // Re-injection: multiply by per-exposure scale factor if provided.
      // This implements CHIMERA's "model × data_scale" approach (Line et al. 2021):
      // the model is embedded into the systematic matrix before SVD projection,
      // placing it in the same subspace as the PCA-cleaned data.
      if (model_scale != nullptr)
        raw_model[exp] *= model_scale[offset * nb_exposures + exp * N + p];
    }

    // Apply (I-P) projection or temporal-mean subtraction, depending on apply_projection.
    //
    // apply_projection=true:  standard path — apply (I-P) to put model and data
    //   in the same filtered temporal subspace.
    //
    // apply_projection=false: data-only filtering path — subtract only the temporal
    //   mean per pixel from the model (equivalent to N_PCA=0 model filtering).
    //   This zeroes out a velocity-independent flat model (CIA/blackbody barely
    //   changes with the small Doppler shifts of a hot Jupiter) while preserving
    //   the Doppler trail of molecular lines.  Kernel 2 then detrends spectrally
    //   and correlates with the fully PCA-cleaned, spectrally detrended data.
    if (apply_projection)
    {
      for (int exp = 0; exp < nb_exposures; ++exp)
      {
        float val = 0.0f;

        for (int j = 0; j < nb_exposures; ++j)
          val += s_IminusP[exp * nb_exposures + j] * raw_model[j];

        model_filtered[offset * nb_exposures + exp * N + p] = val;
      }
    }
    else
    {
      // filter_model=false: store raw model flux (no projection, no normalization).
      // CIA cancels in the total-CCF kernel because CIA is spectrally smooth:
      // a ±Kp Doppler shift moves CIA by <1–2 pixels at R=45000, so
      // Fp_CIA[e,p] ≈ C[p] for all e, and ∑_e D_pca[e,p] × C[p] = C[p] × 0 = 0
      // (per-pixel temporal mean of D_pca is exactly zero by construction).
      // Molecular lines are narrow: they shift by many pixels across exposures,
      // creating the Doppler trail that gives non-zero R_xf at the correct velocity.
      for (int exp = 0; exp < nb_exposures; ++exp)
        model_filtered[offset * nb_exposures + exp * N + p] = raw_model[exp];
    }
  }
}


// Kernel 2: Compute log-likelihood from precomputed filtered model.
// One block per (order, exposure) pair — same grid as the unfiltered kernel.
//
// The CHIMERA PCA-cleaned data has its spectral continuum (mean + slope) removed per
// exposure.  The filtered BeAR model, however, retains a spectral linear trend because
// Doppler-shifting the blackbody continuum by ±95 km/s shifts the spectral slope,
// producing a per-pixel residual proportional to (dB/dλ) * v/c after temporal-mean
// removal.  This slope contributes heavily to rff (model power at all ~1848 pixels)
// without contributing to rxf (data has no slope), diluting r = rxf/√(rff·sf2) and
// burying the molecular signal.
//
// Fix: linearly detrend the model per exposure (fit and subtract mean + slope across
// pixels) so that rff only captures molecular-feature power (~10% of pixels).  The
// data is already spectrally detrended by CHIMERA, so no data detrending is needed.
__global__
void highResLogLikeFromFilteredKernel(
  const float*  __restrict__ model_filtered,
  const float*  __restrict__ order_flux,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const float*  __restrict__ data_mean,
  const double* __restrict__ data_sf2,
  const int                  nb_orders,
  const int                  nb_exposures,
  const float                alpha,
  double*       __restrict__ d_log_like)
{
  const int task = blockIdx.x;
  const int exp  = task / nb_orders;
  const int ord  = task % nb_orders;
  const int tid  = threadIdx.x;

  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];

  // Pointers to this order/exposure data
  const float* flux = order_flux + offset * nb_exposures + exp * N;
  const float* model = model_filtered + offset * nb_exposures + exp * N;

  // --- Pass 1: Compute model quadratic spectral fit coefficients ---
  // Fit: model[p] = a + b*(p-p_mid) + c*(p-p_mid)^2 + residual
  // Normal equations (symmetric x = p - p_mid, so odd cross terms vanish):
  //   a = (sigma4*Sy - sigma2*Sx2y) / det
  //   b = Sxy / sigma2
  //   c = (N*Sx2y - sigma2*Sy) / det
  // sigma2 = N(N^2-1)/12,  sigma4 = N(N^2-1)(3N^2-7)/240,  det = N*sigma4 - sigma2^2
  const float p_mid = 0.5f * (float)(N - 1);

  float model_sum  = 0.0f;
  float model_xsum = 0.0f;
  float model_x2sum = 0.0f;

  for (int p = tid; p < N; p += blockDim.x)
  {
    float x = (float)p - p_mid;
    float m = model[p];
    model_sum   += m;
    model_xsum  += x * m;
    model_x2sum += x * x * m;
  }

  model_sum   = blockReduceSum(model_sum);
  model_xsum  = blockReduceSum(model_xsum);
  model_x2sum = blockReduceSum(model_x2sum);

  __shared__ float s_model_a;
  __shared__ float s_model_b;
  __shared__ float s_model_c;

  if (tid == 0)
  {
    const float dN     = (float)N;
    const float sigma2 = dN * ((float)(N + 1) * (float)(N - 1)) / 12.0f;
    const float sigma4 = dN * ((float)(N + 1) * (float)(N - 1))
                            * (3.0f * dN * dN - 7.0f) / 240.0f;
    const float det    = dN * sigma4 - sigma2 * sigma2;

    s_model_a = (det > 0.0f)    ? (sigma4 * model_sum  - sigma2 * model_x2sum) / det
                                 : model_sum / dN;
    s_model_b = (sigma2 > 0.0f) ? model_xsum / sigma2 : 0.0f;
    s_model_c = (det > 0.0f)    ? (dN * model_x2sum - sigma2 * model_sum) / det : 0.0f;
  }

  __syncthreads();

  // --- Pass 2: Cross-correlation sums ---
  const float ma     = s_model_a;
  const float mb     = s_model_b;
  const float mc     = s_model_c;
  const float dmean  = data_mean[ord * nb_exposures + exp];

  double local_rxf = 0.0;
  double local_rff = 0.0;

  for (int p = tid; p < N; p += blockDim.x)
  {
    float x = (float)p - p_mid;
    double d = (double)flux[p] - (double)dmean;
    // Subtract quadratic spectral fit from model per exposure.
    // This removes the broad CIA/blackbody continuum curvature, leaving only
    // narrow-band molecular features for cross-correlation with the filtered data.
    double m = (double)model[p] - (double)(ma + mb * x + mc * x * x);
    local_rxf += d * m;
    local_rff += m * m;
  }

  local_rxf = blockReduceSum(local_rxf);
  local_rff = blockReduceSum(local_rff);

  // --- Thread 0: Compute log-likelihood ---
  if (tid == 0)
  {
    const double dN = (double)N;
    const double sf2 = data_sf2[ord * nb_exposures + exp];

    if (local_rff > 0.0)
    {
      double arg;

      if (alpha < 0.0f)
      {
        arg = (sf2 - (local_rxf * local_rxf) / (dN * local_rff)) / sf2;
      }
      else
      {
        const double a = (double)alpha;
        const double sg2 = local_rff / dN;
        const double R = local_rxf / dN;
        arg = (sf2 + a * a * sg2 - 2.0 * a * R) / sf2;
      }

      if (arg > 0.0)
        atomicAdd(d_log_like, -0.5 * dN * log(arg));
      else
        atomicAdd(d_log_like, -1e30);
    }
  }
}


// Total-CCF kernel for the filter_model=false path.
// One block per order.  Accumulates R_xf, R_ff, sf2 over ALL exposures × pixels
// and then computes a single log-likelihood contribution per order.
//
// The CIA/blackbody continuum cancels in the total CCF because the (I-P)-filtered
// data has zero temporal mean per pixel: sum_e data[e,p] = 0 for all p.
// Hence sum_{e,p} data[e,p]*CIA[p] = sum_p CIA[p]*(sum_e data[e,p]) = 0.
// This cancellation only holds when the sum is taken over ALL exposures; a
// per-exposure correlation would have a non-zero CIA contribution for each
// individual exposure.
__global__
void highResLogLikeTotalCCFKernel(
  const float*  __restrict__ model_filtered,   // raw Fp, [ord_offset*nb_exp + exp*N + p]
  const float*  __restrict__ order_flux,       // (I-P)-filtered data, same layout
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const int                  nb_orders,
  const int                  nb_exposures,
  double*       __restrict__ d_log_like)
{
  const int ord = blockIdx.x;
  const int tid = threadIdx.x;

  if (ord >= nb_orders) return;

  const int N      = order_nb_pixels[ord];
  const int offset = order_offsets[ord];

  const int total = N * nb_exposures;

  double local_rxf = 0.0;
  double local_rff = 0.0;
  double local_sf2 = 0.0;

  for (int idx = tid; idx < total; idx += blockDim.x)
  {
    const int exp = idx / N;
    const int p   = idx % N;
    const int base = offset * nb_exposures + exp * N + p;

    const double d = (double)order_flux[base];
    const double m = (double)model_filtered[base];

    local_rxf += d * m;
    local_rff += m * m;
    local_sf2 += d * d;
  }

  local_rxf = blockReduceSum(local_rxf);
  local_rff = blockReduceSum(local_rff);
  local_sf2 = blockReduceSum(local_sf2);

  if (tid == 0 && local_rff > 0.0 && local_sf2 > 0.0)
  {
    const double dN = (double)(N * nb_exposures);
    const double r2 = (local_rxf * local_rxf) / (local_sf2 * local_rff);

    if (r2 < 1.0)
      atomicAdd(d_log_like, -0.5 * dN * log(1.0 - r2));
    else
      atomicAdd(d_log_like, -1e30);
  }
}


__host__
void launchHighResLogLikeFiltered(
    const float* broadened_spectrum_dev,
    const double* model_wavelengths_dev,
    int n_model,
    const float* order_wavelengths_dev,
    const float* order_flux_dev,
    const int* order_offsets_dev,
    const int* order_nb_pixels_dev,
    const float* data_mean_dev,
    const double* data_sf2_dev,
    const float* orbital_phases_dev,
    const float* v_bary_dev,
    const float* per_exposure_model_dev,
    const float* projection_matrices_dev,
    float* model_filtered_dev,
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    float Kp, float Vsys, float dphi,
    float alpha,
    double* d_log_like_dev,
    const float* model_scale_dev,
    bool apply_model_projection,
    bool use_phase_function,
    const float* stellar_spectrum_dev)
{
  // Kernel 1: Interpolate + filter, one block per order
  {
    const int threads = 256;
    const int blocks = nb_orders;
    const size_t shared_mem = nb_exposures * nb_exposures * sizeof(float);

    highResInterpFilterKernel<<<blocks, threads, shared_mem>>>(
      broadened_spectrum_dev,
      model_wavelengths_dev,
      n_model,
      order_wavelengths_dev,
      order_offsets_dev,
      order_nb_pixels_dev,
      orbital_phases_dev,
      v_bary_dev,
      per_exposure_model_dev,
      projection_matrices_dev,
      model_filtered_dev,
      nb_orders,
      nb_exposures,
      Kp, Vsys, dphi,
      model_scale_dev,
      apply_model_projection,
      use_phase_function,
      stellar_spectrum_dev);

    CUDA_CHECK_AFTER_KERNEL();
  }

  // Kernel 2: Likelihood from model
  if (apply_model_projection)
  {
    // filter_model=true: per-(order,exposure) logL with quadratic-detrended model
    const int threads = 256;
    const int blocks = nb_orders * nb_exposures;

    highResLogLikeFromFilteredKernel<<<blocks, threads>>>(
      model_filtered_dev,
      order_flux_dev,
      order_offsets_dev,
      order_nb_pixels_dev,
      data_mean_dev,
      data_sf2_dev,
      nb_orders,
      nb_exposures,
      alpha,
      d_log_like_dev);

    CUDA_CHECK_AFTER_KERNEL();
  }
  else
  {
    // filter_model=false: total CCF per order (sum over all exposures).
    // CIA cancels because data has zero temporal mean per pixel.
    const int threads = 256;
    const int blocks = nb_orders;

    highResLogLikeTotalCCFKernel<<<blocks, threads>>>(
      model_filtered_dev,
      order_flux_dev,
      order_offsets_dev,
      order_nb_pixels_dev,
      nb_orders,
      nb_exposures,
      d_log_like_dev);

    CUDA_CHECK_AFTER_KERNEL();
  }
}


// ============================================================================
// Gibson et al. 2022 Eq. 4: per-pixel uncertainty weighting
// ============================================================================


// Unfiltered Gibson kernel: one block per (order, exposure).
// Single pass: interpolate model, accumulate weighted sums Sm, Sfm, Smm.
// Thread 0: combine with precomputed S1, Sf, Sff → chi2 → log-likelihood.
__global__
void highResLogLikeGibsonKernel(
  const float*  __restrict__ broadened_spectrum,
  const double* __restrict__ model_wavelengths,
  const int                  n_model,
  const float*  __restrict__ order_wavelengths,
  const float*  __restrict__ order_flux,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const float*  __restrict__ flux_uncertainties,
  const double* __restrict__ gibson_S1,
  const double* __restrict__ gibson_Sf,
  const double* __restrict__ gibson_Sff,
  const float*  __restrict__ orbital_phases,
  const float*  __restrict__ v_bary,
  const float*  __restrict__ per_exposure_model,
  const int                  nb_orders,
  const int                  nb_exposures,
  const float                Kp,
  const float                Vsys,
  const float                dphi,
  const float                alpha,
  double*       __restrict__ d_log_like)
{
  const int task = blockIdx.x;
  const int exp  = task / nb_orders;
  const int ord  = task % nb_orders;
  const int tid  = threadIdx.x;

  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];

  const float* wl_order = order_wavelengths + offset;
  const float* flux = order_flux + offset * nb_exposures + exp * N;
  const float* sigma = flux_uncertainties + offset * nb_exposures + exp * N;

  // Model source: per-exposure exposure-blurred copy when active, else the shared
  // broadened model (point interpolation, no in-kernel box averaging).
  const float* model_src = (per_exposure_model != nullptr)
    ? per_exposure_model + (size_t)exp * n_model
    : broadened_spectrum;

  // Doppler factor for this exposure
  const float c_kms = (float)(constants::light_c * 1e-5);
  const float phase = orbital_phases[exp];
  const float v_rad = Kp * sinf(2.0f * (float)M_PI * (phase + dphi)) + Vsys + v_bary[exp];
  const double inv_doppler = 1.0 / (1.0 + (double)v_rad / (double)c_kms);

  // Single pass: interpolate model, accumulate weighted sums
  double local_Sm = 0.0, local_Sfm = 0.0, local_Smm = 0.0;

  for (int p = tid; p < N; p += blockDim.x)
  {
    const double wl_base = (double)wl_order[p] * 1e-3;
    const float model_val = interpolateModel(model_src, model_wavelengths, n_model,
                                             wl_base * inv_doppler);

    const double m = (double)model_val;
    const double f = (double)flux[p];
    const double s = (double)sigma[p];
    const double inv_s2 = 1.0 / (s * s);

    local_Sm  += m * inv_s2;
    local_Sfm += f * m * inv_s2;
    local_Smm += m * m * inv_s2;
  }

  local_Sm  = blockReduceSum(local_Sm);
  local_Sfm = blockReduceSum(local_Sfm);
  local_Smm = blockReduceSum(local_Smm);

  if (tid == 0)
  {
    const int idx = ord * nb_exposures + exp;
    const double S1  = gibson_S1[idx];
    const double Sf  = gibson_Sf[idx];
    const double Sff = gibson_Sff[idx];
    const double a = (double)alpha;
    const double dN = (double)N;

    const double chi2_data = Sff - Sf * Sf / S1;
    const double chi2 = chi2_data
                       + a * a * (local_Smm - local_Sm * local_Sm / S1)
                       - 2.0 * a * (local_Sfm - Sf * local_Sm / S1);

    if (chi2 > 0.0 && chi2_data > 0.0)
      atomicAdd(d_log_like, -0.5 * dN * log(chi2 / chi2_data));
    else
      atomicAdd(d_log_like, -1e30);
  }
}


__host__
void launchHighResLogLikeGibson(
    const float* broadened_spectrum_dev,
    const double* model_wavelengths_dev,
    int n_model,
    const float* order_wavelengths_dev,
    const float* order_flux_dev,
    const int* order_offsets_dev,
    const int* order_nb_pixels_dev,
    const float* flux_uncertainties_dev,
    const double* gibson_S1_dev,
    const double* gibson_Sf_dev,
    const double* gibson_Sff_dev,
    const float* orbital_phases_dev,
    const float* v_bary_dev,
    const float* per_exposure_model_dev,
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    float Kp, float Vsys, float dphi,
    float alpha,
    double* d_log_like_dev,
    const float* stellar_spectrum_dev)
{
  const int threads = 256;
  const int blocks = nb_orders * nb_exposures;

  highResLogLikeGibsonKernel<<<blocks, threads>>>(
    broadened_spectrum_dev,
    model_wavelengths_dev,
    n_model,
    order_wavelengths_dev,
    order_flux_dev,
    order_offsets_dev,
    order_nb_pixels_dev,
    flux_uncertainties_dev,
    gibson_S1_dev,
    gibson_Sf_dev,
    gibson_Sff_dev,
    orbital_phases_dev,
    v_bary_dev,
    per_exposure_model_dev,
    nb_orders,
    nb_exposures,
    Kp, Vsys, dphi,
    alpha,
    d_log_like_dev);

  CUDA_CHECK_AFTER_KERNEL();
}


// Kernel 2 for filtered Gibson: read precomputed filtered model and filtered flux,
// apply per-pixel σ_i weighting. One block per (order, exposure).
__global__
void highResLogLikeFromFilteredGibsonKernel(
  const float*  __restrict__ model_filtered,
  const float*  __restrict__ order_flux,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const float*  __restrict__ flux_uncertainties,
  const double* __restrict__ gibson_S1,
  const double* __restrict__ gibson_Sf,
  const double* __restrict__ gibson_Sff,
  const int                  nb_orders,
  const int                  nb_exposures,
  const float                alpha,
  double*       __restrict__ d_log_like)
{
  const int task = blockIdx.x;
  const int exp  = task / nb_orders;
  const int ord  = task % nb_orders;
  const int tid  = threadIdx.x;

  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];

  const float* model = model_filtered + offset * nb_exposures + exp * N;
  const float* flux  = order_flux + offset * nb_exposures + exp * N;
  const float* sigma = flux_uncertainties + offset * nb_exposures + exp * N;

  double local_Sm = 0.0, local_Sfm = 0.0, local_Smm = 0.0;

  for (int p = tid; p < N; p += blockDim.x)
  {
    const double m = (double)model[p];
    const double f = (double)flux[p];
    const double s = (double)sigma[p];
    const double inv_s2 = 1.0 / (s * s);

    local_Sm  += m * inv_s2;
    local_Sfm += f * m * inv_s2;
    local_Smm += m * m * inv_s2;
  }

  local_Sm  = blockReduceSum(local_Sm);
  local_Sfm = blockReduceSum(local_Sfm);
  local_Smm = blockReduceSum(local_Smm);

  if (tid == 0)
  {
    const int idx = ord * nb_exposures + exp;
    const double S1  = gibson_S1[idx];
    const double Sf  = gibson_Sf[idx];
    const double Sff = gibson_Sff[idx];
    const double a = (double)alpha;
    const double dN = (double)N;

    const double chi2_data = Sff - Sf * Sf / S1;
    const double chi2 = chi2_data
                       + a * a * (local_Smm - local_Sm * local_Sm / S1)
                       - 2.0 * a * (local_Sfm - Sf * local_Sm / S1);

    if (chi2 > 0.0 && chi2_data > 0.0)
      atomicAdd(d_log_like, -0.5 * dN * log(chi2 / chi2_data));
    else
      atomicAdd(d_log_like, -1e30);
  }
}


__host__
void launchHighResLogLikeFilteredGibson(
    const float* broadened_spectrum_dev,
    const double* model_wavelengths_dev,
    int n_model,
    const float* order_wavelengths_dev,
    const float* order_flux_dev,
    const int* order_offsets_dev,
    const int* order_nb_pixels_dev,
    const float* orbital_phases_dev,
    const float* v_bary_dev,
    const float* per_exposure_model_dev,
    const float* projection_matrices_dev,
    float* model_filtered_dev,
    const float* flux_uncertainties_dev,
    const double* gibson_S1_dev,
    const double* gibson_Sf_dev,
    const double* gibson_Sff_dev,
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    float Kp, float Vsys, float dphi,
    float alpha,
    double* d_log_like_dev,
    const float* stellar_spectrum_dev)
{
  // Kernel 1: Interpolate + filter (reuse existing kernel)
  {
    const int threads = 256;
    const int blocks = nb_orders;
    const size_t shared_mem = nb_exposures * nb_exposures * sizeof(float);

    highResInterpFilterKernel<<<blocks, threads, shared_mem>>>(
      broadened_spectrum_dev,
      model_wavelengths_dev,
      n_model,
      order_wavelengths_dev,
      order_offsets_dev,
      order_nb_pixels_dev,
      orbital_phases_dev,
      v_bary_dev,
      per_exposure_model_dev,
      projection_matrices_dev,
      model_filtered_dev,
      nb_orders,
      nb_exposures,
      Kp, Vsys, dphi,
      nullptr,   // model_scale
      true,      // apply_projection
      false,     // use_phase_function
      stellar_spectrum_dev);

    CUDA_CHECK_AFTER_KERNEL();
  }

  // Kernel 2: Gibson likelihood from filtered model
  {
    const int threads = 256;
    const int blocks = nb_orders * nb_exposures;

    highResLogLikeFromFilteredGibsonKernel<<<blocks, threads>>>(
      model_filtered_dev,
      order_flux_dev,
      order_offsets_dev,
      order_nb_pixels_dev,
      flux_uncertainties_dev,
      gibson_S1_dev,
      gibson_Sf_dev,
      gibson_Sff_dev,
      nb_orders,
      nb_exposures,
      alpha,
      d_log_like_dev);

    CUDA_CHECK_AFTER_KERNEL();
  }
}


} // namespace bear
