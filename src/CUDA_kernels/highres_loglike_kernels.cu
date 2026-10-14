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
#include "../observations/highres_model_transform.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>


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


// Same bracket as binarySearchDescending, but found by a short linear scan
// from a nearby hint index. Used when consecutive lookups target nearly
// identical wavelengths (e.g. the same pixel at consecutive exposures, whose
// Doppler shifts differ by only a few km/s ~ a few grid points).
__device__ __forceinline__
int bracketFromHint(
  const double* __restrict__ arr,
  int n,
  double val,
  int hint)
{
  int i = min(max(hint, 0), n - 1);

  while (i > 0 && arr[i - 1] <= val) --i;
  while (i < n - 1 && arr[i] > val) ++i;

  return i;
}


// Linear interpolation of a spectrum on the (descending) model grid, also
// returning the bracket index so that a second spectrum can be interpolated at
// the same wavelength without another binary search.
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


// Per-exposure geometry, shared by all kernels (and mirrored on the CPU):
// Doppler factor, full intra-exposure velocity smear, and phase weight.
// Kp and Vsys are the total velocities (reference + retrieved offset).
__device__ __forceinline__
void exposureGeometry(
  const double phase_obs,
  const double v_bary,
  const float* __restrict__ exposure_blur_coeff,
  const int exp,
  const double Kp,
  const double Vsys,
  const double dphi,
  const bool use_phase_function,
  double& inv_dop,
  double& delta_v,
  double& phase_w)
{
  const double c_kms = constants::light_c * 1e-5;
  const double phase = phase_obs + dphi;
  const double v_rad = Kp * sin(2.0 * M_PI * phase) + Vsys + v_bary;

  inv_dop = 1.0 / (1.0 + v_rad / c_kms);

  // Only the orbital term accelerates during an exposure; Vsys and v_bary are
  // constant over it.  Zero when blurring is disabled.
  delta_v = (exposure_blur_coeff != nullptr)
    ? Kp * cos(2.0 * M_PI * phase) * (double)exposure_blur_coeff[exp]
    : 0.0;

  phase_w = use_phase_function ? highResPhaseWeight(phase) : 1.0;
}


// Model value at observed wavelength wl_rest (µm) for one exposure:
//   * box-averaged over the intra-exposure smear delta_v, or point-interpolated
//     when the smear is negligible, at the Doppler target wl_rest * inv_dop;
//   * emission (stellar_spectrum != nullptr): multiplied by Fs(λ')/Fs(λ_rest), so
//     that only Fp is Doppler shifted and Fs stays in the stellar rest frame;
//     fs_rest = Fs(λ_rest) is computed once per pixel by the caller;
//   * transmission (stellar_spectrum == nullptr): converted from transit depth
//     in ppm to the relative in-transit flux change, -1e-6 depth.
// Returns 0 outside the model grid.
// idx_hint: bracket index of a previous lookup at a nearby wavelength (-1 if
// none); the bracket is then found by a short scan instead of a binary search.
// Updated to the bracket found here.
__device__ __forceinline__
float modelValueAtTarget(
  const float*  __restrict__ model_src,
  const double* __restrict__ model_wavelengths,
  const int n_model,
  const double wl_rest,
  const double inv_dop,
  const double delta_v,
  const float*  __restrict__ stellar_spectrum,
  const float fs_rest,
  int& idx_hint)
{
  const double wl_dop = wl_rest * inv_dop;
  const int idx = (idx_hint < 0)
    ? binarySearchDescending(model_wavelengths, n_model, wl_dop)
    : bracketFromHint(model_wavelengths, n_model, wl_dop, idx_hint);

  idx_hint = idx;

  if (idx < 1 || idx >= n_model || model_wavelengths[idx] > wl_dop)
    return 0.0f;

  float m;

  if (fabs(delta_v) > exposure_blur_min_kms)
  {
    // Box-at-target: edges at v_rad ± dV/2 (1/inv_dop = 1 + v_rad/c).
    const double c_kms  = constants::light_c * 1e-5;
    const double half   = 0.5 * fabs(delta_v);
    const double inv_lo = 1.0 / (1.0 / inv_dop + half / c_kms);  // larger v -> smaller inv
    const double inv_hi = 1.0 / (1.0 / inv_dop - half / c_kms);
    double mean;

    if (!boxMeanRaw(model_src, model_wavelengths, n_model,
                    wl_rest * inv_lo, wl_rest * inv_hi, mean, idx))
      return 0.0f;

    m = (float)mean;
  }
  else
  {
    m = interpolateAtIdx(model_src, model_wavelengths, n_model, idx, wl_dop);
  }

  if (stellar_spectrum == nullptr)
    return (float)(transit_depth_to_flux * (double)m);

  if (fs_rest > 0.0f)
    m *= interpolateAtIdx(stellar_spectrum, model_wavelengths, n_model, idx, wl_dop) / fs_rest;

  return m;
}


__device__ __forceinline__
float modelValueAtTarget(
  const float*  __restrict__ model_src,
  const double* __restrict__ model_wavelengths,
  const int n_model,
  const double wl_rest,
  const double inv_dop,
  const double delta_v,
  const float*  __restrict__ stellar_spectrum,
  const float fs_rest)
{
  int idx_hint = -1;

  return modelValueAtTarget(model_src, model_wavelengths, n_model, wl_rest,
                            inv_dop, delta_v, stellar_spectrum, fs_rest, idx_hint);
}


__device__ __forceinline__
float stellarRestValue(
  const float*  __restrict__ stellar_spectrum,
  const double* __restrict__ model_wavelengths,
  const int n_model,
  const double wl_rest)
{
  if (stellar_spectrum == nullptr)
    return 0.0f;

  int idx;
  const float fs = interpolateModelRawIdx(stellar_spectrum, model_wavelengths, n_model, wl_rest, idx);

  return (idx >= 1 && idx < n_model && model_wavelengths[idx] <= wl_rest) ? fs : 0.0f;
}


// Model value at the first pixel of an order, used as a reference that is
// subtracted from every model value before the sums are accumulated.  Both
// unfiltered likelihood forms are exactly invariant to a constant shift of the
// model; subtracting one keeps float rounding of the data (e.g. ~0.97 for a
// transmission spectrum) from being multiplied by a large model offset (e.g. the
// -1e-6 x transit depth continuum).  Must be called by all threads of the block.
__device__
double blockModelReference(
  const float*  __restrict__ model_src,
  const double* __restrict__ model_wavelengths,
  const int n_model,
  const double* __restrict__ wl_order,
  const double inv_dop,
  const double delta_v,
  const double phase_w,
  const float*  __restrict__ stellar_spectrum)
{
  __shared__ double s_ref;

  if (threadIdx.x == 0)
  {
    const double wl_rest = wl_order[0] * 1e-3;
    const float fs_rest = stellarRestValue(stellar_spectrum, model_wavelengths, n_model, wl_rest);
    s_ref = phase_w * (double)modelValueAtTarget(
      model_src, model_wavelengths, n_model, wl_rest, inv_dop, delta_v, stellar_spectrum, fs_rest);
  }

  __syncthreads();

  return s_ref;
}


// Least-squares quadratic a + b x + c x^2 in the centred pixel index
// x = p - (N-1)/2 of y[0..N-1].  Must be called by all threads of the block;
// every thread receives the coefficients.
__device__
void blockQuadraticFit(
  const float* __restrict__ y,
  const int N,
  double& a,
  double& b,
  double& c)
{
  const double p_mid = 0.5 * (double)(N - 1);

  double sy = 0.0, sxy = 0.0, sx2y = 0.0;

  for (int p = threadIdx.x; p < N; p += blockDim.x)
  {
    const double x = (double)p - p_mid;
    const double v = (double)y[p];
    sy   += v;
    sxy  += x * v;
    sx2y += x * x * v;
  }

  sy   = blockReduceSum(sy);
  sxy  = blockReduceSum(sxy);
  sx2y = blockReduceSum(sx2y);

  __shared__ double s_coeff[3];

  if (threadIdx.x == 0)
    quadraticFitFromSums((double)N, sy, sxy, sx2y, s_coeff[0], s_coeff[1], s_coeff[2]);

  __syncthreads();

  a = s_coeff[0];
  b = s_coeff[1];
  c = s_coeff[2];
}


// ============================================================================
// Unfiltered path
// ============================================================================


// One block per (order, exposure).  Single pass: interpolate the model and
// accumulate raw sums; the mean-subtracted sums follow algebraically:
//   R_xf = sum(d*m) - dmean*sum(m) - mmean*sum(d) + N*dmean*mmean
//   R_ff = sum(m^2) - N*mmean^2
__global__
void highResLogLikeKernel(
  const float*  __restrict__ broadened_spectrum,
  const double* __restrict__ model_wavelengths,
  const int                  n_model,
  const double* __restrict__ order_wavelengths,
  const float*  __restrict__ order_flux,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const float*  __restrict__ data_mean,
  const double* __restrict__ data_sf2,
  const double* __restrict__ orbital_phases,
  const double* __restrict__ v_bary,
  const float*  __restrict__ exposure_blur_coeff,
  const int                  nb_orders,
  const int                  nb_exposures,
  const double               Kp,
  const double               Vsys,
  const double               dphi,
  const double               alpha,
  const bool                 use_phase_function,
  const float*  __restrict__ stellar_spectrum,
  double*       __restrict__ d_log_like)
{
  const int task = blockIdx.x;
  const int exp  = task / nb_orders;
  const int ord  = task % nb_orders;
  const int tid  = threadIdx.x;

  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];

  const double* wl_order = order_wavelengths + offset;
  const float* flux = order_flux + (size_t)offset * nb_exposures + (size_t)exp * N;

  double inv_dop, delta_v, phase_w;
  exposureGeometry(orbital_phases[exp], v_bary[exp], exposure_blur_coeff, exp,
                   Kp, Vsys, dphi, use_phase_function, inv_dop, delta_v, phase_w);

  const double m_ref = blockModelReference(
    broadened_spectrum, model_wavelengths, n_model, wl_order, inv_dop, delta_v,
    phase_w, stellar_spectrum);

  double sum_m = 0.0, sum_m2 = 0.0, sum_dm = 0.0, sum_d = 0.0;

  for (int p = tid; p < N; p += blockDim.x)
  {
    const double wl_rest = wl_order[p] * 1e-3;
    const float fs_rest = stellarRestValue(stellar_spectrum, model_wavelengths, n_model, wl_rest);
    const double m = phase_w * (double)modelValueAtTarget(
      broadened_spectrum, model_wavelengths, n_model, wl_rest, inv_dop, delta_v,
      stellar_spectrum, fs_rest) - m_ref;
    const double d = (double)flux[p];

    sum_m  += m;
    sum_m2 += m * m;
    sum_dm += d * m;
    sum_d  += d;
  }

  sum_m  = blockReduceSum(sum_m);
  sum_m2 = blockReduceSum(sum_m2);
  sum_dm = blockReduceSum(sum_dm);
  sum_d  = blockReduceSum(sum_d);

  if (tid == 0)
  {
    const double dN = (double)N;
    const double sf2 = data_sf2[ord * nb_exposures + exp];
    const double mmean = sum_m / dN;
    const double dmean = (double)data_mean[ord * nb_exposures + exp];

    const double rxf = sum_dm - dmean * sum_m - mmean * sum_d + dN * dmean * mmean;
    const double rff = sum_m2 - dN * mmean * mmean;

    if (rff > 0.0)
      atomicAdd(d_log_like, brogiLineLogLike(rxf, rff, sf2, dN, alpha));
  }
}


__host__
void launchHighResLogLike(
    const float* broadened_spectrum_dev,
    const double* model_wavelengths_dev,
    int n_model,
    const double* order_wavelengths_dev,
    const float* order_flux_dev,
    const int* order_offsets_dev,
    const int* order_nb_pixels_dev,
    const float* data_mean_dev,
    const double* data_sf2_dev,
    const double* orbital_phases_dev,
    const double* v_bary_dev,
    const float* exposure_blur_coeff_dev,
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    double Kp, double Vsys, double dphi,
    double alpha,
    double* d_log_like_dev,
    const float* stellar_spectrum_dev,
    bool use_phase_function)
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
    exposure_blur_coeff_dev,
    nb_orders,
    nb_exposures,
    Kp, Vsys, dphi,
    alpha,
    use_phase_function,
    stellar_spectrum_dev,
    d_log_like_dev);

  CUDA_CHECK_AFTER_KERNEL();
}


// Unfiltered Gibson kernel: one block per (order, exposure).
// Single pass: interpolate the model, accumulate weighted sums Sm, Sfm, Smm;
// thread 0 combines them with the precomputed S1, Sf, Sff.
__global__
void highResLogLikeGibsonKernel(
  const float*  __restrict__ broadened_spectrum,
  const double* __restrict__ model_wavelengths,
  const int                  n_model,
  const double* __restrict__ order_wavelengths,
  const float*  __restrict__ order_flux,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const float*  __restrict__ flux_uncertainties,
  const double* __restrict__ gibson_S1,
  const double* __restrict__ gibson_Sf,
  const double* __restrict__ gibson_Sff,
  const double* __restrict__ orbital_phases,
  const double* __restrict__ v_bary,
  const float*  __restrict__ exposure_blur_coeff,
  const int                  nb_orders,
  const int                  nb_exposures,
  const double               Kp,
  const double               Vsys,
  const double               dphi,
  const double               alpha,
  const bool                 use_phase_function,
  const float*  __restrict__ stellar_spectrum,
  double*       __restrict__ d_log_like)
{
  const int task = blockIdx.x;
  const int exp  = task / nb_orders;
  const int ord  = task % nb_orders;
  const int tid  = threadIdx.x;

  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];

  const double* wl_order = order_wavelengths + offset;
  const float* flux  = order_flux + (size_t)offset * nb_exposures + (size_t)exp * N;
  const float* sigma = flux_uncertainties + (size_t)offset * nb_exposures + (size_t)exp * N;

  double inv_dop, delta_v, phase_w;
  exposureGeometry(orbital_phases[exp], v_bary[exp], exposure_blur_coeff, exp,
                   Kp, Vsys, dphi, use_phase_function, inv_dop, delta_v, phase_w);

  const double m_ref = blockModelReference(
    broadened_spectrum, model_wavelengths, n_model, wl_order, inv_dop, delta_v,
    phase_w, stellar_spectrum);

  double local_Sm = 0.0, local_Sfm = 0.0, local_Smm = 0.0;

  for (int p = tid; p < N; p += blockDim.x)
  {
    const double wl_rest = wl_order[p] * 1e-3;
    const float fs_rest = stellarRestValue(stellar_spectrum, model_wavelengths, n_model, wl_rest);
    const double m = phase_w * (double)modelValueAtTarget(
      broadened_spectrum, model_wavelengths, n_model, wl_rest, inv_dop, delta_v,
      stellar_spectrum, fs_rest) - m_ref;
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

    atomicAdd(d_log_like, gibsonLogLike(
      gibson_S1[idx], gibson_Sf[idx], gibson_Sff[idx],
      local_Sm, local_Sfm, local_Smm, (double)N, alpha));
  }
}


__host__
void launchHighResLogLikeGibson(
    const float* broadened_spectrum_dev,
    const double* model_wavelengths_dev,
    int n_model,
    const double* order_wavelengths_dev,
    const float* order_flux_dev,
    const int* order_offsets_dev,
    const int* order_nb_pixels_dev,
    const float* flux_uncertainties_dev,
    const double* gibson_S1_dev,
    const double* gibson_Sf_dev,
    const double* gibson_Sff_dev,
    const double* orbital_phases_dev,
    const double* v_bary_dev,
    const float* exposure_blur_coeff_dev,
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    double Kp, double Vsys, double dphi,
    double alpha,
    double* d_log_like_dev,
    const float* stellar_spectrum_dev,
    bool use_phase_function)
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
    exposure_blur_coeff_dev,
    nb_orders,
    nb_exposures,
    Kp, Vsys, dphi,
    alpha,
    use_phase_function,
    stellar_spectrum_dev,
    d_log_like_dev);

  CUDA_CHECK_AFTER_KERNEL();
}


// ============================================================================
// Filtered path: Gibson et al. 2022 fast model filtering, and data-only mode
// ============================================================================


// Kernel 1: interpolate the model at all Doppler shifts, then treat each pixel
// column of the order according to `projection`.
// 2D grid: blockIdx.x = order, blockIdx.y = pixel chunk, so the GPU is saturated
// even for a small number of orders.  Each thread processes one pixel column and
// loops over all exposures; the column lives in the local array raw[MAX_K].
//
// Shared memory: the per-exposure Doppler factors, smears and phase weights
// (3 K doubles), followed by the (I - P) matrix (K^2 floats) when it fits.
// Otherwise the matrix is read from global memory (it stays cached in L1/L2).
//
// Output: model_out[ord_offset * K + exp * N + p]
template <int MAX_K>
__global__
void highResInterpFilterKernel(
  const float*  __restrict__ broadened_spectrum,
  const double* __restrict__ model_wavelengths,
  const int                  n_model,
  const double* __restrict__ order_wavelengths,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const double* __restrict__ orbital_phases,
  const double* __restrict__ v_bary,
  const float*  __restrict__ exposure_blur_coeff,
  const float*  __restrict__ projection_matrices,
  float*        __restrict__ model_out,
  const int                  nb_orders,
  const int                  nb_exposures,
  const double               Kp,
  const double               Vsys,
  const double               dphi,
  const float*  __restrict__ model_scale,
  const int                  projection,
  const bool                 use_phase_function,
  const float*  __restrict__ stellar_spectrum,
  const bool                 matrix_in_shared)
{
  const int ord = blockIdx.x;
  const int tid = threadIdx.x;

  if (ord >= nb_orders) return;

  const int K = nb_exposures;
  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];
  const double* wl_order = order_wavelengths + offset;

  extern __shared__ double s_shared_raw[];
  double* s_inv_dop = s_shared_raw;
  double* s_delta_v = s_shared_raw + K;
  double* s_phase_w = s_shared_raw + 2 * K;
  float*  s_matrix  = (float*)(s_shared_raw + 3 * K);

  const float* proj_ptr = projection_matrices + (size_t)ord * K * K;
  const float* IminusP = matrix_in_shared ? s_matrix : proj_ptr;

  if (projection == highres_projection_filter && matrix_in_shared)
    for (int i = tid; i < K * K; i += blockDim.x)
      s_matrix[i] = proj_ptr[i];

  for (int e = tid; e < K; e += blockDim.x)
    exposureGeometry(orbital_phases[e], v_bary[e], exposure_blur_coeff, e,
                     Kp, Vsys, dphi, use_phase_function,
                     s_inv_dop[e], s_delta_v[e], s_phase_w[e]);

  __syncthreads();

  for (int p = blockIdx.y * blockDim.x + tid; p < N; p += blockDim.x * gridDim.y)
  {
    const double wl_rest = wl_order[p] * 1e-3;

    // Fs(λ_rest) is the same for all exposures of this pixel.
    const float fs_rest = stellarRestValue(stellar_spectrum, model_wavelengths, n_model, wl_rest);

    float raw[MAX_K];

    // Bracket index of the previous exposure; consecutive exposures shift the
    // wavelength by only a few grid points, so after the first (binary) search
    // each subsequent bracket is found by a short scan from the previous one.
    int idx_prev = -1;

    for (int e = 0; e < K; ++e)
    {
      float m = (float)(s_phase_w[e] * (double)modelValueAtTarget(
        broadened_spectrum, model_wavelengths, n_model, wl_rest,
        s_inv_dop[e], s_delta_v[e], stellar_spectrum, fs_rest, idx_prev));

      // Re-injection (Line et al. 2021): put the model in detector units by
      // multiplying with the background the filter captures, S = P D.
      if (model_scale != nullptr)
        m *= model_scale[(size_t)offset * K + (size_t)e * N + p];

      raw[e] = m;
    }

    float* out = model_out + (size_t)offset * K + p;

    if (projection == highres_projection_filter)
    {
      for (int e = 0; e < K; ++e)
      {
        const float* row = IminusP + e * K;
        float val = 0.0f;

        for (int j = 0; j < K; ++j)
          val += row[j] * raw[j];

        out[(size_t)e * N] = val;
      }
    }
    else if (projection == highres_projection_center_time)
    {
      // Data-only mode: remove only the temporal mean of each pixel (N_PCA = 0
      // model filtering).  A continuum that does not move with the planet cancels,
      // the Doppler trail of the lines is kept.
      double mean = 0.0;

      for (int e = 0; e < K; ++e)
        mean += raw[e];

      mean /= (double)K;

      for (int e = 0; e < K; ++e)
        out[(size_t)e * N] = (float)((double)raw[e] - mean);
    }
    else
    {
      for (int e = 0; e < K; ++e)
        out[(size_t)e * N] = raw[e];
    }
  }
}


// Largest dynamic shared memory a block may use on this device after opting in.
static int maxDynamicSharedMemoryBytes()
{
  static int value = -1;

  if (value < 0)
  {
    int device = 0;
    gpuErrchk(cudaGetDevice(&device));
    gpuErrchk(cudaDeviceGetAttribute(&value, cudaDevAttrMaxSharedMemoryPerBlockOptin, device));
  }

  return value;
}


template <int MAX_K, typename... Args>
static void launchInterpFilterInstance(
  const dim3 blocks, const int threads, const size_t shared_bytes, Args... args)
{
  // Above 48 KB, dynamic shared memory needs an explicit opt-in per kernel.
  static size_t opted_in_bytes = 48 * 1024;

  if (shared_bytes > opted_in_bytes)
  {
    gpuErrchk(cudaFuncSetAttribute(
      highResInterpFilterKernel<MAX_K>,
      cudaFuncAttributeMaxDynamicSharedMemorySize,
      static_cast<int>(shared_bytes)));
    opted_in_bytes = shared_bytes;
  }

  highResInterpFilterKernel<MAX_K><<<blocks, threads, shared_bytes>>>(args...);
}


static void launchInterpFilter(
    const float* broadened_spectrum_dev,
    const double* model_wavelengths_dev,
    const int n_model,
    const double* order_wavelengths_dev,
    const int* order_offsets_dev,
    const int* order_nb_pixels_dev,
    const double* orbital_phases_dev,
    const double* v_bary_dev,
    const float* exposure_blur_coeff_dev,
    const float* projection_matrices_dev,
    float* model_out_dev,
    const int nb_orders,
    const int nb_exposures,
    const int max_pixels_per_order,
    const double Kp, const double Vsys, const double dphi,
    const float* model_scale_dev,
    const int projection,
    const bool use_phase_function,
    const float* stellar_spectrum_dev)
{
  const int threads = 256;
  const dim3 blocks(nb_orders, (max_pixels_per_order + threads - 1) / threads);

  const size_t K = static_cast<size_t>(nb_exposures);
  const size_t base_bytes = 3 * K * sizeof(double);
  const size_t matrix_bytes = K * K * sizeof(float);

  const bool matrix_in_shared = projection == highres_projection_filter
    && base_bytes + matrix_bytes <= static_cast<size_t>(maxDynamicSharedMemoryBytes());
  const size_t shared_bytes = base_bytes + (matrix_in_shared ? matrix_bytes : 0);

  #define BEAR_INTERP_FILTER_ARGS \
    broadened_spectrum_dev, model_wavelengths_dev, n_model, order_wavelengths_dev, \
    order_offsets_dev, order_nb_pixels_dev, orbital_phases_dev, v_bary_dev, \
    exposure_blur_coeff_dev, projection_matrices_dev, model_out_dev, nb_orders, \
    nb_exposures, Kp, Vsys, dphi, model_scale_dev, projection, use_phase_function, \
    stellar_spectrum_dev, matrix_in_shared

  if (nb_exposures <= 128)
    launchInterpFilterInstance<128>(blocks, threads, shared_bytes, BEAR_INTERP_FILTER_ARGS);
  else if (nb_exposures <= 256)
    launchInterpFilterInstance<256>(blocks, threads, shared_bytes, BEAR_INTERP_FILTER_ARGS);
  else if (nb_exposures <= highres_max_exposures_gpu)
    launchInterpFilterInstance<highres_max_exposures_gpu>(
      blocks, threads, shared_bytes, BEAR_INTERP_FILTER_ARGS);
  else
    throw std::runtime_error(
      "launchInterpFilter: " + std::to_string(nb_exposures) + " exposures exceed the GPU limit of "
      + std::to_string(highres_max_exposures_gpu) + "\n");

  #undef BEAR_INTERP_FILTER_ARGS

  CUDA_CHECK_AFTER_KERNEL();
}


// Kernel 2 (filter_model = true, Brogi & Line forms): one block per (order,
// exposure).  A quadratic in pixel index is removed from the filtered model,
// which takes out the continuum curvature that survives (I - P), e.g. the
// blackbody slope shifted by ±Kp.  The filtered data were detrended the same
// way once at initialisation.
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
  const double               alpha,
  double*       __restrict__ d_log_like)
{
  const int task = blockIdx.x;
  const int exp  = task / nb_orders;
  const int ord  = task % nb_orders;
  const int tid  = threadIdx.x;

  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];

  const float* flux  = order_flux + (size_t)offset * nb_exposures + (size_t)exp * N;
  const float* model = model_filtered + (size_t)offset * nb_exposures + (size_t)exp * N;

  double ma, mb, mc;
  blockQuadraticFit(model, N, ma, mb, mc);

  const double p_mid = 0.5 * (double)(N - 1);
  const double dmean = (double)data_mean[ord * nb_exposures + exp];

  double local_rxf = 0.0;
  double local_rff = 0.0;

  for (int p = tid; p < N; p += blockDim.x)
  {
    const double x = (double)p - p_mid;
    const double d = (double)flux[p] - dmean;
    const double m = (double)model[p] - (ma + mb * x + mc * x * x);
    local_rxf += d * m;
    local_rff += m * m;
  }

  local_rxf = blockReduceSum(local_rxf);
  local_rff = blockReduceSum(local_rff);

  if (tid == 0)
  {
    const double sf2 = data_sf2[ord * nb_exposures + exp];

    if (local_rff > 0.0)
      atomicAdd(d_log_like, brogiLineLogLike(local_rxf, local_rff, sf2, (double)N, alpha));
  }
}


// Kernel 2 (filter_model = true, Gibson form): one block per (order, exposure).
// The model gets the same quadratic detrending as in the Brogi & Line kernel, so
// that data and model are treated alike; the inverse-variance weighted means are
// then removed through the S_f/S_1 and S_m/S_1 terms.
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
  const double               alpha,
  double*       __restrict__ d_log_like)
{
  const int task = blockIdx.x;
  const int exp  = task / nb_orders;
  const int ord  = task % nb_orders;
  const int tid  = threadIdx.x;

  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];

  const float* model = model_filtered + (size_t)offset * nb_exposures + (size_t)exp * N;
  const float* flux  = order_flux + (size_t)offset * nb_exposures + (size_t)exp * N;
  const float* sigma = flux_uncertainties + (size_t)offset * nb_exposures + (size_t)exp * N;

  double ma, mb, mc;
  blockQuadraticFit(model, N, ma, mb, mc);

  const double p_mid = 0.5 * (double)(N - 1);

  double local_Sm = 0.0, local_Sfm = 0.0, local_Smm = 0.0;

  for (int p = tid; p < N; p += blockDim.x)
  {
    const double x = (double)p - p_mid;
    const double m = (double)model[p] - (ma + mb * x + mc * x * x);
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

    atomicAdd(d_log_like, gibsonLogLike(
      gibson_S1[idx], gibson_Sf[idx], gibson_Sff[idx],
      local_Sm, local_Sfm, local_Smm, (double)N, alpha));
  }
}


// Data-only mode (filter_model = false), stage 1.  One block per (order,
// exposure); each block reduces its exposure's sums and atomically adds them to
// the per-order partials [ord * 6 + k]:
//   k = 0..2: R_xf = sum d m, R_ff = sum m^2, sum d^2          (Brogi & Line forms)
//   k = 3..5: S_m, S_fm, S_mm with 1/sigma^2 weights          (Gibson form only)
// The model is centred in time per pixel (kernel 1) and the (I - P)-filtered data
// have zero temporal mean per pixel, so a continuum that does not move with the
// planet cancels from every sum.
__global__
void highResTotalCCFSumsKernel(
  const float*  __restrict__ model_centred,
  const float*  __restrict__ order_flux,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const float*  __restrict__ flux_uncertainties,
  const int                  nb_orders,
  const int                  nb_exposures,
  double*       __restrict__ order_partials)
{
  const int task = blockIdx.x;
  const int exp  = task / nb_orders;
  const int ord  = task % nb_orders;
  const int tid  = threadIdx.x;

  const int N      = order_nb_pixels[ord];
  const int offset = order_offsets[ord];

  const float* flux  = order_flux + (size_t)offset * nb_exposures + (size_t)exp * N;
  const float* model = model_centred + (size_t)offset * nb_exposures + (size_t)exp * N;
  const float* sigma = flux_uncertainties != nullptr
    ? flux_uncertainties + (size_t)offset * nb_exposures + (size_t)exp * N
    : nullptr;

  double rxf = 0.0, rff = 0.0, sff = 0.0;
  double sm = 0.0, sfm = 0.0, smm = 0.0;

  for (int p = tid; p < N; p += blockDim.x)
  {
    const double d = (double)flux[p];
    const double m = (double)model[p];

    rxf += d * m;
    rff += m * m;
    sff += d * d;

    if (sigma != nullptr)
    {
      const double inv_s2 = 1.0 / ((double)sigma[p] * (double)sigma[p]);
      sm  += m * inv_s2;
      sfm += d * m * inv_s2;
      smm += m * m * inv_s2;
    }
  }

  rxf = blockReduceSum(rxf);
  rff = blockReduceSum(rff);
  sff = blockReduceSum(sff);

  if (sigma != nullptr)
  {
    sm  = blockReduceSum(sm);
    sfm = blockReduceSum(sfm);
    smm = blockReduceSum(smm);
  }

  if (tid == 0)
  {
    double* partials = order_partials + ord * 6;
    atomicAdd(&partials[0], rxf);
    atomicAdd(&partials[1], rff);
    atomicAdd(&partials[2], sff);

    if (sigma != nullptr)
    {
      atomicAdd(&partials[3], sm);
      atomicAdd(&partials[4], sfm);
      atomicAdd(&partials[5], smm);
    }
  }
}


// Data-only mode, stage 2: one log-likelihood per order from the partials, with
// N = pixels x exposures of the order.
__global__
void highResTotalCCFLogLikeKernel(
  const double* __restrict__ order_partials,
  const int*    __restrict__ order_nb_pixels,
  const double* __restrict__ gibson_S1,
  const double* __restrict__ gibson_Sf,
  const double* __restrict__ gibson_Sff,
  const int                  nb_orders,
  const int                  nb_exposures,
  const int                  likelihood_form,
  const double               alpha,
  double*       __restrict__ d_log_like)
{
  for (int ord = blockIdx.x * blockDim.x + threadIdx.x; ord < nb_orders; ord += blockDim.x * gridDim.x)
  {
    const double* partials = order_partials + ord * 6;
    const double dN = (double)order_nb_pixels[ord] * (double)nb_exposures;

    double log_like = 0.0;

    if (likelihood_form == highres_form_gibson)
    {
      double S1 = 0.0, Sf = 0.0, Sff = 0.0;

      for (int e = 0; e < nb_exposures; ++e)
      {
        S1  += gibson_S1[ord * nb_exposures + e];
        Sf  += gibson_Sf[ord * nb_exposures + e];
        Sff += gibson_Sff[ord * nb_exposures + e];
      }

      log_like = gibsonLogLike(S1, Sf, Sff, partials[3], partials[4], partials[5], dN, alpha);
    }
    else
    {
      const double a = (likelihood_form == highres_form_free_alpha) ? alpha : -1.0;
      log_like = brogiLineLogLike(partials[0], partials[1], partials[2] / dN, dN, a);
    }

    atomicAdd(d_log_like, log_like);
  }
}


__host__
void launchHighResLogLikeFiltered(
    const float* broadened_spectrum_dev,
    const double* model_wavelengths_dev,
    int n_model,
    const double* order_wavelengths_dev,
    const float* order_flux_dev,
    const int* order_offsets_dev,
    const int* order_nb_pixels_dev,
    const float* data_mean_dev,
    const double* data_sf2_dev,
    const double* orbital_phases_dev,
    const double* v_bary_dev,
    const float* exposure_blur_coeff_dev,
    const float* projection_matrices_dev,
    float* model_filtered_dev,
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    double Kp, double Vsys, double dphi,
    double alpha,
    int likelihood_form,
    const float* flux_uncertainties_dev,
    const double* gibson_S1_dev,
    const double* gibson_Sf_dev,
    const double* gibson_Sff_dev,
    const float* model_scale_dev,
    bool filter_model,
    bool use_phase_function,
    const float* stellar_spectrum_dev,
    double* d_log_like_dev)
{
  // Kernel 1: interpolate, weight, scale, and filter or centre the model
  launchInterpFilter(
    broadened_spectrum_dev,
    model_wavelengths_dev,
    n_model,
    order_wavelengths_dev,
    order_offsets_dev,
    order_nb_pixels_dev,
    orbital_phases_dev,
    v_bary_dev,
    exposure_blur_coeff_dev,
    projection_matrices_dev,
    model_filtered_dev,
    nb_orders,
    nb_exposures,
    max_pixels_per_order,
    Kp, Vsys, dphi,
    model_scale_dev,
    filter_model ? highres_projection_filter : highres_projection_center_time,
    use_phase_function,
    stellar_spectrum_dev);

  const int threads = 256;
  const bool gibson = likelihood_form == highres_form_gibson;

  if (filter_model)
  {
    // Kernel 2: per (order, exposure) likelihood
    const int blocks = nb_orders * nb_exposures;

    if (gibson)
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
    else
      highResLogLikeFromFilteredKernel<<<blocks, threads>>>(
        model_filtered_dev,
        order_flux_dev,
        order_offsets_dev,
        order_nb_pixels_dev,
        data_mean_dev,
        data_sf2_dev,
        nb_orders,
        nb_exposures,
        likelihood_form == highres_form_free_alpha ? alpha : -1.0,
        d_log_like_dev);

    CUDA_CHECK_AFTER_KERNEL();
  }
  else
  {
    // Data-only mode: sums over all exposures of an order, one logarithm per order.
    // The small per-order partials buffer is cached between calls.
    static double* order_partials_dev = nullptr;
    static int order_partials_capacity = 0;

    if (order_partials_capacity < nb_orders)
    {
      if (order_partials_dev != nullptr)
        gpuErrchk(cudaFree(order_partials_dev));

      gpuErrchk(cudaMalloc(&order_partials_dev, nb_orders * 6 * sizeof(double)));
      order_partials_capacity = nb_orders;
    }

    gpuErrchk(cudaMemset(order_partials_dev, 0, nb_orders * 6 * sizeof(double)));

    highResTotalCCFSumsKernel<<<nb_orders * nb_exposures, threads>>>(
      model_filtered_dev,
      order_flux_dev,
      order_offsets_dev,
      order_nb_pixels_dev,
      gibson ? flux_uncertainties_dev : nullptr,
      nb_orders,
      nb_exposures,
      order_partials_dev);

    CUDA_CHECK_AFTER_KERNEL();

    highResTotalCCFLogLikeKernel<<<1, 128>>>(
      order_partials_dev,
      order_nb_pixels_dev,
      gibson_S1_dev,
      gibson_Sf_dev,
      gibson_Sff_dev,
      nb_orders,
      nb_exposures,
      likelihood_form,
      alpha,
      d_log_like_dev);

    CUDA_CHECK_AFTER_KERNEL();
  }
}


__host__
void launchHighResInterpFilterOnly(
    const float* broadened_spectrum_dev,
    const double* model_wavelengths_dev,
    int n_model,
    const double* order_wavelengths_dev,
    const int* order_offsets_dev,
    const int* order_nb_pixels_dev,
    const double* orbital_phases_dev,
    const double* v_bary_dev,
    const float* exposure_blur_coeff_dev,
    const float* projection_matrices_dev,
    float* model_filtered_dev,
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    double Kp, double Vsys, double dphi,
    const float* model_scale_dev,
    int projection,
    bool use_phase_function,
    const float* stellar_spectrum_dev)
{
  launchInterpFilter(
    broadened_spectrum_dev,
    model_wavelengths_dev,
    n_model,
    order_wavelengths_dev,
    order_offsets_dev,
    order_nb_pixels_dev,
    orbital_phases_dev,
    v_bary_dev,
    exposure_blur_coeff_dev,
    projection_matrices_dev,
    model_filtered_dev,
    nb_orders,
    nb_exposures,
    max_pixels_per_order,
    Kp, Vsys, dphi,
    model_scale_dev,
    projection,
    use_phase_function,
    stellar_spectrum_dev);
}


} // namespace bear
