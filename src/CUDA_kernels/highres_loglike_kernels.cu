// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-FileCopyrightText: 2026 Esther Wang
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "highres_loglike_kernels.h"
#include "reduce_kernels.h"
#include "error_check.h"
#include "../additional/physical_const.h"
#include "../observations/highres_model_transform.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>


namespace bear {


// ============================================================================
// Model lookup on the precomputed per-order slices (see HighResModelGridGPU)
// ============================================================================


// Slice of the model grid owned by one order, with all pointers resolved.
struct OrderSlice {
  const float* ln;        // ln(lambda) - ln_ref, descending (hi part)
  const float* ln_lo;     // low part
  const float* inv_dln;   // 1 / cell width in ln(lambda)
  const float* dln;       // cell width in ln(lambda)
  int start;              // model index of ln[0]
  int len;
  float inv_mean_dln;
};


__device__ __forceinline__
OrderSlice orderSlice(const HighResModelGridGPU& g, const int ord)
{
  OrderSlice s;

  const int off = g.slice_offset[ord];
  s.ln = g.slice_ln + off;
  s.ln_lo = g.slice_ln_lo + off;
  s.inv_dln = g.slice_inv_dln + off;
  s.dln = g.slice_dln + off;
  s.start = g.slice_start[ord];
  s.len = g.slice_len[ord];
  s.inv_mean_dln = g.order_inv_mean_dln[ord];

  return s;
}


// A position in ln(lambda) as a double-float pair: hi carries the float value,
// lo the rounding remainder, so that hi + lo holds ~16 significant digits.
struct Pos {
  float hi;
  float lo;
};


// Exact sum of two positions (two-sum of the hi parts, the error joins the lo parts).
__device__ __forceinline__
Pos addPos(const Pos a, const Pos b)
{
  const float s = a.hi + b.hi;
  const float bb = s - a.hi;
  const float err = (a.hi - (s - bb)) + (b.hi - bb);

  Pos r;
  r.hi = s;
  r.lo = err + (a.lo + b.lo);

  return r;
}


// Signed distance ln[j] - u of a slice point from a position, > 0 when the
// grid point lies above u.  Only meaningful when the two are within a grid cell
// or so, where the hi parts subtract exactly (Sterbenz) and lo carries the rest.
__device__ __forceinline__
float distPos(const OrderSlice& sl, const int j, const Pos u)
{
  return (sl.ln[j] - u.hi) + (sl.ln_lo[j] - u.lo);
}


// Per-exposure Doppler geometry in the units of the slices: the shift of
// ln(lambda) at the exposure's radial velocity, the shifts of the two edges of
// the exposure-blur box, the inverse box width and the phase weight.
struct ExposureShift {
  Pos shift;         // ln(1 / (1 + v_rad/c))
  Pos shift_lo;      // same for the box edge at v_rad + |dV|/2 (smaller shift)
  Pos shift_hi;      // same for the box edge at v_rad - |dV|/2 (larger shift)
  float inv_box;     // 1 / (shift_hi - shift_lo)
  float phase_w;
  int blur;          // 1: box-average over the exposure smear, 0: point value
};


// Bracket of u on the descending slice: the first index j with ln[j] <= u,
// found by a short scan of the hi parts from a nearby hint, then settled with
// the exact distances.  The bracket is valid when 1 <= j < len and u < ln[0];
// the callers test the latter themselves.
__device__ __forceinline__
int bracketScan(const OrderSlice& sl, const Pos u, const int hint)
{
  int j = min(max(hint, 1), sl.len - 1);

  while (j > 1 && sl.ln[j - 1] <= u.hi) --j;
  while (j < sl.len && sl.ln[j] > u.hi) ++j;

  // the hi-part scan can be off by one where u is within a float ulp of a grid point
  if (j > 1 && distPos(sl, j - 1, u) <= 0.0f)
    --j;
  else if (j < sl.len && distPos(sl, j, u) > 0.0f)
    ++j;

  return j;
}


// Linear interpolation of spectrum within cell j of the slice at position u.
__device__ __forceinline__
float lerpInCell(
  const float* __restrict__ spectrum,
  const OrderSlice& sl,
  const int j,
  const Pos u)
{
  const float t = distPos(sl, j - 1, u) * sl.inv_dln[j];
  const float* s = spectrum + sl.start;

  return (1.0f - t) * s[j - 1] + t * s[j];
}


// Mean of the spectrum over the box [u_lo, u_hi] in ln(lambda), as a trapezoidal
// integral over the slice cells.  j is the bracket of the box centre, from which
// both edges are found by short scans.  Over a box a few cells wide, integrating
// in ln(lambda) instead of lambda changes the mean by O(dV/c)^2, far below float
// resolution.  Returns false when an edge leaves the slice.
__device__ __forceinline__
bool boxMeanOnSlice(
  const float* __restrict__ spectrum,
  const OrderSlice& sl,
  const int j,
  const Pos u_lo,
  const Pos u_hi,
  const float inv_box,
  float& out_mean)
{
  if (distPos(sl, 0, u_hi) <= 0.0f) return false;

  const int j_hi = bracketScan(sl, u_hi, j);
  const int j_lo = bracketScan(sl, u_lo, j);

  if (j_hi < 1 || j_lo >= sl.len) return false;

  const float* s = spectrum + sl.start;

  const float t_hi = distPos(sl, j_hi - 1, u_hi) * sl.inv_dln[j_hi];
  const float f_hi = (1.0f - t_hi) * s[j_hi - 1] + t_hi * s[j_hi];

  const float t_lo = distPos(sl, j_lo - 1, u_lo) * sl.inv_dln[j_lo];
  const float f_lo = (1.0f - t_lo) * s[j_lo - 1] + t_lo * s[j_lo];

  if (j_hi == j_lo)
  {
    // both edges share a cell
    out_mean = 0.5f * (f_lo + f_hi);
    return true;
  }

  // top partial cell: u_hi down to grid point j_hi
  float integral = 0.5f * (f_hi + s[j_hi]) * (1.0f - t_hi) * sl.dln[j_hi];

  // full interior cells
  for (int i = j_hi; i <= j_lo - 2; ++i)
    integral += 0.5f * (s[i] + s[i + 1]) * sl.dln[i + 1];

  // bottom partial cell: grid point j_lo - 1 down to u_lo
  integral += 0.5f * (s[j_lo - 1] + f_lo) * t_lo * sl.dln[j_lo];

  out_mean = integral * inv_box;

  return true;
}


// Model value of one pixel at one exposure:
//   * point-interpolated at the Doppler-shifted position, or box-averaged over
//     the intra-exposure smear;
//   * emission (stellar != nullptr): multiplied by Fs(lambda')/Fs(lambda_rest), so
//     that only Fp is Doppler shifted and Fs stays in the stellar rest frame;
//     inv_fs_rest = 1/Fs(lambda_rest), computed once per pixel by the caller, or
//     a negative value when Fs(lambda_rest) is not positive (no correction);
//   * transmission (stellar == nullptr): transit depth in ppm converted to the
//     relative in-transit flux change.
// Returns 0 outside the slice.  j_hint is the bracket of a previous lookup at a
// nearby position (or a prediction) and is updated to the bracket found here.
__device__ __forceinline__
float modelValueOnSlice(
  const float* __restrict__ spectrum,
  const float* __restrict__ stellar,
  const OrderSlice& sl,
  const Pos u_rest,
  const ExposureShift& ex,
  const float inv_fs_rest,
  int& j_hint)
{
  const Pos u = addPos(u_rest, ex.shift);

  const int j = bracketScan(sl, u, j_hint);
  j_hint = j;

  if (j < 1 || j >= sl.len || distPos(sl, 0, u) <= 0.0f)
    return 0.0f;

  float m;

  if (ex.blur)
  {
    if (!boxMeanOnSlice(spectrum, sl, j, addPos(u_rest, ex.shift_lo), addPos(u_rest, ex.shift_hi),
                        ex.inv_box, m))
      return 0.0f;
  }
  else
  {
    m = lerpInCell(spectrum, sl, j, u);
  }

  if (stellar == nullptr)
    return static_cast<float>(transit_depth_to_flux) * m;

  if (inv_fs_rest > 0.0f)
    m *= lerpInCell(stellar, sl, j, u) * inv_fs_rest;

  return m;
}


// 1/Fs at the rest wavelength of a pixel (emission), or -1 when there is no
// stellar spectrum, the pixel lies outside the slice or Fs is not positive.
__device__ __forceinline__
float inverseStellarRest(
  const float* __restrict__ stellar,
  const OrderSlice& sl,
  const int j0,
  const float t0)
{
  if (stellar == nullptr || j0 < 1 || j0 >= sl.len)
    return -1.0f;

  const float* s = stellar + sl.start;
  const float fs = (1.0f - t0) * s[j0 - 1] + t0 * s[j0];

  return fs > 0.0f ? 1.0f / fs : -1.0f;
}


// Rest position of a pixel as a double-float pair.
__device__ __forceinline__
Pos pixelPosition(const HighResModelGridGPU& g, const int p)
{
  Pos u;
  u.hi = g.pix_ln[p];
  u.lo = g.pix_ln_lo[p];
  return u;
}


// Bracket prediction for a pixel at a given exposure shift from its bracket at
// zero velocity: the grid is descending, so a positive shift moves to lower indices.
__device__ __forceinline__
int predictBracket(const int j0, const float shift, const OrderSlice& sl)
{
  const int j = (j0 >= 0) ? j0 : sl.len / 2;

  return j - __float2int_rn(shift * sl.inv_mean_dln);
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


__device__ __forceinline__
Pos splitPosition(const double x)
{
  Pos p;
  p.hi = (float)x;
  p.lo = (float)(x - (double)p.hi);
  return p;
}


// The same geometry converted to the ln(lambda) shifts of the slices.  The box
// edges sit at v_rad -/+ |dV|/2 (1/inv_dop = 1 + v_rad/c).
__device__ __forceinline__
ExposureShift exposureShift(
  const double phase_obs,
  const double v_bary,
  const float* __restrict__ exposure_blur_coeff,
  const int exp,
  const double Kp,
  const double Vsys,
  const double dphi,
  const bool use_phase_function)
{
  double inv_dop, delta_v, phase_w;
  exposureGeometry(phase_obs, v_bary, exposure_blur_coeff, exp,
                   Kp, Vsys, dphi, use_phase_function, inv_dop, delta_v, phase_w);

  ExposureShift ex;
  ex.shift = splitPosition(log(inv_dop));
  ex.phase_w = (float)phase_w;
  ex.blur = fabs(delta_v) > exposure_blur_min_kms ? 1 : 0;
  ex.shift_lo = ex.shift;
  ex.shift_hi = ex.shift;
  ex.inv_box = 0.0f;

  if (ex.blur)
  {
    const double c_kms = constants::light_c * 1e-5;
    const double half = 0.5 * fabs(delta_v);
    const double ln_lo = -log(1.0 / inv_dop + half / c_kms);
    const double ln_hi = -log(1.0 / inv_dop - half / c_kms);

    ex.shift_lo = splitPosition(ln_lo);
    ex.shift_hi = splitPosition(ln_hi);
    ex.inv_box = (float)(1.0 / (ln_hi - ln_lo));
  }

  return ex;
}


// Model value at the first pixel of an order, used as a reference that is
// subtracted from every model value before the sums are accumulated.  Both
// unfiltered likelihood forms are exactly invariant to a constant shift of the
// model; subtracting one keeps float rounding of the data (e.g. ~0.97 for a
// transmission spectrum) from being multiplied by a large model offset (e.g. the
// -1e-6 x transit depth continuum).  Must be called by all threads of the block.
__device__
float blockModelReference(
  const float*  __restrict__ model_src,
  const float*  __restrict__ stellar_spectrum,
  const HighResModelGridGPU& grid,
  const OrderSlice& sl,
  const int offset,
  const ExposureShift& ex)
{
  __shared__ float s_ref;

  if (threadIdx.x == 0)
  {
    const int j0 = grid.pix_j0[offset];
    const float inv_fs = inverseStellarRest(stellar_spectrum, sl, j0, grid.pix_t0[offset]);
    int j = predictBracket(j0, ex.shift.hi, sl);

    s_ref = ex.phase_w * modelValueOnSlice(
      model_src, stellar_spectrum, sl, pixelPosition(grid, offset), ex, inv_fs, j);
  }

  __syncthreads();

  return s_ref;
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
  const HighResModelGridGPU  grid,
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

  const float* flux = order_flux + (size_t)offset * nb_exposures + (size_t)exp * N;

  const OrderSlice sl = orderSlice(grid, ord);
  const ExposureShift ex = exposureShift(orbital_phases[exp], v_bary[exp], exposure_blur_coeff, exp,
                                         Kp, Vsys, dphi, use_phase_function);

  const float m_ref = blockModelReference(broadened_spectrum, stellar_spectrum, grid, sl, offset, ex);

  double sum_m = 0.0, sum_m2 = 0.0, sum_dm = 0.0, sum_d = 0.0;

  for (int p = tid; p < N; p += blockDim.x)
  {
    const int j0 = grid.pix_j0[offset + p];
    const float inv_fs = inverseStellarRest(stellar_spectrum, sl, j0, grid.pix_t0[offset + p]);
    int j = predictBracket(j0, ex.shift.hi, sl);

    const double m = (double)(ex.phase_w * modelValueOnSlice(
      broadened_spectrum, stellar_spectrum, sl, pixelPosition(grid, offset + p), ex, inv_fs, j) - m_ref);
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
    const HighResModelGridGPU& model_grid,
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
    model_grid,
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
  const HighResModelGridGPU  grid,
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

  const float* flux  = order_flux + (size_t)offset * nb_exposures + (size_t)exp * N;
  const float* sigma = flux_uncertainties + (size_t)offset * nb_exposures + (size_t)exp * N;

  const OrderSlice sl = orderSlice(grid, ord);
  const ExposureShift ex = exposureShift(orbital_phases[exp], v_bary[exp], exposure_blur_coeff, exp,
                                         Kp, Vsys, dphi, use_phase_function);

  const float m_ref = blockModelReference(broadened_spectrum, stellar_spectrum, grid, sl, offset, ex);

  double local_Sm = 0.0, local_Sfm = 0.0, local_Smm = 0.0;

  for (int p = tid; p < N; p += blockDim.x)
  {
    const int j0 = grid.pix_j0[offset + p];
    const float inv_fs = inverseStellarRest(stellar_spectrum, sl, j0, grid.pix_t0[offset + p]);
    int j = predictBracket(j0, ex.shift.hi, sl);

    const double m = (double)(ex.phase_w * modelValueOnSlice(
      broadened_spectrum, stellar_spectrum, sl, pixelPosition(grid, offset + p), ex, inv_fs, j) - m_ref);
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
    const HighResModelGridGPU& model_grid,
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
    model_grid,
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
//
// Kernel 1a interpolates the raw model matrix of every order.  In filter mode
// the pixel moments of the raw exposures then give, through (I - P), the
// quadratic that detrends each filtered exposure.  Kernel 1b applies (I - P)
// (or removes the temporal mean in data-only mode) tile by tile, subtracts the
// quadratic in registers and accumulates the likelihood sums per tile; the
// filtered model itself is only written out for the diagnostic export.
// Kernel 2 combines the tiles and evaluates the logarithms.


// Kernel 1a: interpolate the model at all Doppler shifts, apply the phase weight
// and the optional re-injection scale, and write the raw model matrix of the order.
// 2D grid: blockIdx.x = order, blockIdx.y = pixel chunk, so the GPU is saturated
// even for a small number of orders.  Each thread processes one pixel and loops
// over all exposures, chaining the bracket from one exposure to the next.
//
// Shared memory: the per-exposure shifts (K ExposureShift).
//
// Output: model_out[ord_offset * K + exp * N + p]
__global__
void highResInterpKernel(
  const float*  __restrict__ broadened_spectrum,
  const HighResModelGridGPU  grid,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const double* __restrict__ orbital_phases,
  const double* __restrict__ v_bary,
  const float*  __restrict__ exposure_blur_coeff,
  float*        __restrict__ model_out,
  const int                  nb_orders,
  const int                  nb_exposures,
  const double               Kp,
  const double               Vsys,
  const double               dphi,
  const float*  __restrict__ model_scale,
  const bool                 use_phase_function,
  const float*  __restrict__ stellar_spectrum)
{
  const int ord = blockIdx.x;
  const int tid = threadIdx.x;

  if (ord >= nb_orders) return;

  const int K = nb_exposures;
  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];

  extern __shared__ double s_shared_raw[];
  ExposureShift* s_ex = (ExposureShift*)s_shared_raw;

  for (int e = tid; e < K; e += blockDim.x)
    s_ex[e] = exposureShift(orbital_phases[e], v_bary[e], exposure_blur_coeff, e,
                            Kp, Vsys, dphi, use_phase_function);

  __syncthreads();

  const OrderSlice sl = orderSlice(grid, ord);

  for (int p = blockIdx.y * blockDim.x + tid; p < N; p += blockDim.x * gridDim.y)
  {
    const Pos u_rest = pixelPosition(grid, offset + p);
    const int j0 = grid.pix_j0[offset + p];

    // Fs(λ_rest) is the same for all exposures of this pixel.
    const float inv_fs = inverseStellarRest(stellar_spectrum, sl, j0, grid.pix_t0[offset + p]);

    float* out = model_out + (size_t)offset * K + p;

    // Bracket of the previous exposure; consecutive exposures shift the position
    // by only a few grid points, so after the first prediction each subsequent
    // bracket is found by a short scan from the previous one.
    int j = predictBracket(j0, s_ex[0].shift.hi, sl);

    for (int e = 0; e < K; ++e)
    {
      float m = s_ex[e].phase_w * modelValueOnSlice(
        broadened_spectrum, stellar_spectrum, sl, u_rest, s_ex[e], inv_fs, j);

      // Re-injection (Line et al. 2021): put the model in detector units by
      // multiplying with the background the filter captures, S = P D.
      if (model_scale != nullptr)
        m *= model_scale[(size_t)offset * K + (size_t)e * N + p];

      out[(size_t)e * N] = m;
    }
  }
}


// Compensated float accumulation of products: the product is split exactly with
// an FMA and the running sum with a two-sum, so that the total reaches double
// accuracy at float cost.
struct CompensatedSum {
  float sum = 0.0f;
  float err = 0.0f;

  __device__ __forceinline__ void add(const float a, const float b)
  {
    const float p = a * b;
    const float p_err = fmaf(a, b, -p);
    const float t = sum + p;
    const float bb = t - sum;

    err += ((sum - (t - bb)) + (p - bb)) + p_err;
    sum = t;
  }

  __device__ __forceinline__ double value() const
  {
    return (double)sum + (double)err;
  }
};


// Kernel 1b-pre, filter mode: pixel moments sum m, sum x m, sum x^2 m of every raw
// exposure (x the centred pixel index), one warp per (order, exposure).  The
// quadratic detrending of the filtered model needs these moments of the
// filtered exposures, which follow from the raw ones through the same (I - P)
// (both are linear in the model); the fit is then known before the filter runs
// and the residual can be formed in registers, see highResProjectKernel.
// The filter removes the static part of the model, and with it all but ~1e-4 of
// the raw moments, so these are accumulated with compensated arithmetic: plain
// float sums would leave the fit coefficients with errors of ~1e-3.  Rounding
// that is the same for every exposure, like that of x^2 for very long orders,
// passes through (I - P) unamplified and is harmless.
__global__
void highResRawMomentsKernel(
  const float*  __restrict__ model_raw,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const int                  nb_exposures,
  double*       __restrict__ raw_moments)
{
  const int ord = blockIdx.x;
  const int warp = threadIdx.x >> 5;
  const int lane = threadIdx.x & 31;
  const int e = blockIdx.y * (blockDim.x >> 5) + warp;

  if (e >= nb_exposures) return;

  const int N = order_nb_pixels[ord];
  const float* column = model_raw + (size_t)order_offsets[ord] * nb_exposures + (size_t)e * N;
  const float p_mid = 0.5f * (float)(N - 1);

  CompensatedSum s0, s1, s2;

  for (int p = lane; p < N; p += 32)
  {
    const float m = column[p];
    const float x = (float)p - p_mid;

    s0.add(m, 1.0f);
    s1.add(x, m);
    s2.add(x * x, m);
  }

  const double m0 = warpReduceSum(s0.value());
  const double m1 = warpReduceSum(s1.value());
  const double m2 = warpReduceSum(s2.value());

  if (lane == 0)
  {
    double* out = raw_moments + ((size_t)ord * nb_exposures + e) * 3;
    out[0] = m0;
    out[1] = m1;
    out[2] = m2;
  }
}


// Coefficients a, b, c of the least-squares quadratic of every filtered
// exposure, from the raw moments and the transposed (I - P); one thread per
// (order, exposure).  Errors in the coefficients act on the likelihood sums only
// to second order (the residual is orthogonal to the quadratics).
__global__
void highResFitCoefficientsKernel(
  const double* __restrict__ raw_moments,
  const float*  __restrict__ projection_matrices_t,
  const int*    __restrict__ order_nb_pixels,
  const int                  nb_orders,
  const int                  nb_exposures,
  const int                  nb_exposures_padded,
  float4*       __restrict__ fit_coeff)
{
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int ord = idx / nb_exposures;
  const int e = idx - ord * nb_exposures;

  if (ord >= nb_orders) return;

  const float* matrix_t =
    projection_matrices_t + (size_t)ord * nb_exposures * nb_exposures_padded + e;
  const double* moments = raw_moments + (size_t)ord * nb_exposures * 3;

  double s0 = 0.0, s1 = 0.0, s2 = 0.0;

  #pragma unroll 8
  for (int j = 0; j < nb_exposures; ++j)
  {
    const double w = __ldg(matrix_t + (size_t)j * nb_exposures_padded);

    s0 += w * moments[j * 3];
    s1 += w * moments[j * 3 + 1];
    s2 += w * moments[j * 3 + 2];
  }

  double a, b, c;
  quadraticFitFromSums((double)order_nb_pixels[ord], s0, s1, s2, a, b, c);

  fit_coeff[idx] = make_float4((float)a, (float)b, (float)c, 0.0f);
}


// Kernel 1b: temporal treatment of the raw model matrix of an order and the
// likelihood sums.  One block per (order, tile of TILE pixels); the tile's K x TILE
// raw values are staged in shared memory, laid out [exposure][pixel] so that the
// reads are conflict-free.  The block has TILE x PROJECT_GROUPS threads: the groups
// split the exposures, and each thread accumulates four exposures per shared-memory
// read from a float4 of the transposed, padded projection matrix (layout: see
// initDeviceMemory).
//   projection_filter:      (I - P) raw, then minus the quadratic of the exposure
//   projection_center_time: minus the temporal mean of each pixel
// LIKELIHOOD_SUMS = false writes the treated model back in place (diagnostics).
// LIKELIHOOD_SUMS = true keeps it in registers and accumulates the sums the
// likelihood needs over the tile, in float within a warp and in double across the
// warps, per (order, exposure) in filter mode:
//   Brogi & Line forms:  sum d r, sum r^2                       (r = detrended model)
//   Gibson form:         sum r/s^2, sum f r/s^2, sum r^2/s^2
// and per order in data-only mode, over all its exposures:
//   sum d r, sum r^2, sum d^2, and the three Gibson sums.
// The data d carry the per-exposure data mean (Brogi & Line) already removed.
constexpr int highres_project_groups = 4;
constexpr int highres_nb_sums = 3;          // partial sums per (order, exposure), filter mode
constexpr int highres_nb_order_sums = 6;    // partial sums per order, data-only mode


// Bytes of dynamic shared memory of the projection kernel: the model tile, the
// quadratic coefficients and data means of the exposures, the warp partial sums.
static size_t projectSharedBytes(const int tile, const int nb_exposures)
{
  const int warps_per_group = tile / 32;
  const int nb_warps = warps_per_group * highres_project_groups;
  const int partial_floats = std::max(
    nb_exposures * highres_nb_sums * warps_per_group, nb_warps * highres_nb_order_sums);

  return static_cast<size_t>(nb_exposures) * tile * sizeof(float)
       + static_cast<size_t>(nb_exposures) * sizeof(float4)
       + static_cast<size_t>(partial_floats) * sizeof(float);
}


template <int TILE, bool LIKELIHOOD_SUMS>
__global__
void highResProjectKernel(
  const float*  __restrict__ projection_matrices_t,
  float*        __restrict__ model_inout,
  const float*  __restrict__ order_flux,
  const float*  __restrict__ flux_uncertainties,
  const float*  __restrict__ data_mean,
  const float4* __restrict__ fit_coeff,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  double*       __restrict__ partials,
  const int                  nb_orders,
  const int                  nb_exposures,
  const int                  nb_exposures_padded,
  const int                  nb_tiles,
  const int                  projection)
{
  constexpr int warps_per_group = TILE / 32;
  constexpr int nb_warps = warps_per_group * highres_project_groups;
  constexpr int nb_threads = nb_warps * 32;

  const int ord = blockIdx.x;
  const int tile = blockIdx.y;
  const int tid = threadIdx.x;
  const int group = threadIdx.y;
  const int flat = group * TILE + tid;
  const int lane = tid & 31;
  const int warp = tid >> 5;   // within the group

  if (ord >= nb_orders) return;

  const int K = nb_exposures;
  const int Kp = nb_exposures_padded;
  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];
  const int p = tile * TILE + tid;

  extern __shared__ float s_shared[];
  float* s_tile = s_shared;                                                    // [K][TILE]
  float4* s_coeff = reinterpret_cast<float4*>(s_tile + (size_t)K * TILE);     // [K]
  float* s_part = reinterpret_cast<float*>(s_coeff + K);                      // warp partials

  const size_t column_offset = (size_t)offset * K + p;
  float* column = model_inout + column_offset;
  const float* flux = order_flux + column_offset;
  const float* sigma = (flux_uncertainties != nullptr) ? flux_uncertainties + column_offset : nullptr;

  // threads beyond the order's last pixel only take part in the barriers and
  // reductions, with zero contributions
  const bool active = p < N;

  // the model, data and uncertainties are streamed once (evict-first); the
  // projection matrix is read K times per tile
  if (active)
    for (int e = group; e < K; e += highres_project_groups)
      s_tile[e * TILE + tid] = __ldcs(column + (size_t)e * N);

  if (LIKELIHOOD_SUMS && projection == highres_projection_filter)
    for (int e = flat; e < K; e += nb_threads)
    {
      float4 c = fit_coeff[ord * K + e];
      c.w = (data_mean != nullptr) ? data_mean[ord * K + e] : 0.0f;
      s_coeff[e] = c;
    }

  __syncthreads();

  if (projection == highres_projection_filter)
  {
    const float* matrix_t = projection_matrices_t + (size_t)ord * K * Kp;
    const float x = (float)p - 0.5f * (float)(N - 1);

    // exposures of this group, in blocks of four (the padding makes the float4
    // loads aligned and the extra accumulators harmless)
    const int chunk = ((Kp / 4 + highres_project_groups - 1) / highres_project_groups) * 4;
    const int e_begin = group * chunk;
    const int e_end = min(e_begin + chunk, Kp);

    for (int e0 = e_begin; e0 < e_end; e0 += 4)
    {
      float4 acc = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
      const float4* mt = reinterpret_cast<const float4*>(matrix_t + e0);
      const int stride4 = Kp / 4;

      for (int jj = 0; jj < K; ++jj)
      {
        const float v = s_tile[jj * TILE + tid];
        const float4 m = __ldg(mt + jj * stride4);

        acc.x += m.x * v;
        acc.y += m.y * v;
        acc.z += m.z * v;
        acc.w += m.w * v;
      }

      if (!LIKELIHOOD_SUMS)
      {
        if (active)
        {
          column[(size_t)e0 * N] = acc.x;
          if (e0 + 1 < K) column[(size_t)(e0 + 1) * N] = acc.y;
          if (e0 + 2 < K) column[(size_t)(e0 + 2) * N] = acc.z;
          if (e0 + 3 < K) column[(size_t)(e0 + 3) * N] = acc.w;
        }

        continue;
      }

      const float filtered[4] = {acc.x, acc.y, acc.z, acc.w};

      #pragma unroll
      for (int i = 0; i < 4; ++i)
      {
        const int e = e0 + i;

        if (e >= K) break;   // uniform across the block

        float v0 = 0.0f, v1 = 0.0f, v2 = 0.0f;

        if (active)
        {
          const float4 c = s_coeff[e];
          const float r = filtered[i] - (c.x + x * (c.y + x * c.z));
          const float d = __ldcs(flux + (size_t)e * N) - c.w;

          if (sigma != nullptr)
          {
            const float s = __ldcs(sigma + (size_t)e * N);
            const float w = 1.0f / (s * s);

            v0 = r * w;
            v1 = d * r * w;
            v2 = r * r * w;
          }
          else
          {
            v0 = d * r;
            v1 = r * r;
          }
        }

        v0 = warpReduceSum(v0);
        v1 = warpReduceSum(v1);
        v2 = warpReduceSum(v2);

        if (lane == 0)
        {
          float* part = s_part + (e * highres_nb_sums) * warps_per_group + warp;

          part[0] = v0;
          part[warps_per_group] = v1;
          part[2 * warps_per_group] = v2;
        }
      }
    }

    if (!LIKELIHOOD_SUMS) return;

    __syncthreads();

    for (int e = flat; e < K; e += nb_threads)
    {
      double* out = partials + (((size_t)ord * K + e) * nb_tiles + tile) * highres_nb_sums;

      #pragma unroll
      for (int k = 0; k < highres_nb_sums; ++k)
      {
        const float* part = s_part + (e * highres_nb_sums + k) * warps_per_group;
        double s = 0.0;

        for (int w = 0; w < warps_per_group; ++w)
          s += part[w];

        out[k] = s;
      }
    }
  }
  else
  {
    // Data-only mode: remove only the temporal mean of each pixel (N_PCA = 0
    // model filtering).  A continuum that does not move with the planet cancels,
    // the Doppler trail of the lines is kept.
    float mean = 0.0f;

    for (int e = 0; e < K; ++e)
      mean += s_tile[e * TILE + tid];

    mean /= (float)K;

    if (!LIKELIHOOD_SUMS)
    {
      if (active)
        for (int e = group; e < K; e += highres_project_groups)
          column[(size_t)e * N] = s_tile[e * TILE + tid] - mean;

      return;
    }

    float v[highres_nb_order_sums] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

    if (active)
      for (int e = group; e < K; e += highres_project_groups)
      {
        const float r = s_tile[e * TILE + tid] - mean;
        const float d = __ldcs(flux + (size_t)e * N);

        v[0] += d * r;
        v[1] += r * r;
        v[2] += d * d;

        if (sigma != nullptr)
        {
          const float s = __ldcs(sigma + (size_t)e * N);
          const float w = 1.0f / (s * s);

          v[3] += r * w;
          v[4] += d * r * w;
          v[5] += r * r * w;
        }
      }

    const int warp_flat = flat >> 5;

    #pragma unroll
    for (int k = 0; k < highres_nb_order_sums; ++k)
    {
      v[k] = warpReduceSum(v[k]);

      if (lane == 0) s_part[k * nb_warps + warp_flat] = v[k];
    }

    __syncthreads();

    if (flat < highres_nb_order_sums)
    {
      double s = 0.0;

      for (int w = 0; w < nb_warps; ++w)
        s += s_part[flat * nb_warps + w];

      partials[((size_t)ord * nb_tiles + tile) * highres_nb_order_sums + flat] = s;
    }
  }
}


// Kernel 2, filter mode: one log-likelihood per (order, exposure) from the tile
// partials.
__global__
void highResFilteredLogLikeKernel(
  const double* __restrict__ partials,
  const int*    __restrict__ order_nb_pixels,
  const double* __restrict__ data_sf2,
  const double* __restrict__ gibson_S1,
  const double* __restrict__ gibson_Sf,
  const double* __restrict__ gibson_Sff,
  const int                  nb_orders,
  const int                  nb_exposures,
  const int                  nb_tiles,
  const int                  likelihood_form,
  const double               alpha,
  double*       __restrict__ d_log_like)
{
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int ord = idx / nb_exposures;

  if (ord >= nb_orders) return;

  const double* part = partials + (size_t)idx * nb_tiles * highres_nb_sums;

  double s0 = 0.0, s1 = 0.0, s2 = 0.0;

  for (int t = 0; t < nb_tiles; ++t)
  {
    s0 += part[t * highres_nb_sums];
    s1 += part[t * highres_nb_sums + 1];
    s2 += part[t * highres_nb_sums + 2];
  }

  const double dN = (double)order_nb_pixels[ord];
  double log_like;

  if (likelihood_form == highres_form_gibson)
  {
    log_like = gibsonLogLike(gibson_S1[idx], gibson_Sf[idx], gibson_Sff[idx], s0, s1, s2, dN, alpha);
  }
  else
  {
    if (!(s1 > 0.0)) return;

    const double a = (likelihood_form == highres_form_free_alpha) ? alpha : -1.0;
    log_like = brogiLineLogLike(s0, s1, data_sf2[idx], dN, a);
  }

  atomicAdd(d_log_like, log_like);
}


// Kernel 2, data-only mode: one log-likelihood per order from the tile partials,
// with N = pixels x exposures of the order.
__global__
void highResTotalCCFLogLikeKernel(
  const double* __restrict__ partials,
  const int*    __restrict__ order_nb_pixels,
  const double* __restrict__ gibson_S1,
  const double* __restrict__ gibson_Sf,
  const double* __restrict__ gibson_Sff,
  const int                  nb_orders,
  const int                  nb_exposures,
  const int                  nb_tiles,
  const int                  likelihood_form,
  const double               alpha,
  double*       __restrict__ d_log_like)
{
  for (int ord = blockIdx.x * blockDim.x + threadIdx.x; ord < nb_orders; ord += blockDim.x * gridDim.x)
  {
    const double* part = partials + (size_t)ord * nb_tiles * highres_nb_order_sums;
    double sums[highres_nb_order_sums] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (int t = 0; t < nb_tiles; ++t)
      for (int k = 0; k < highres_nb_order_sums; ++k)
        sums[k] += part[t * highres_nb_order_sums + k];

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

      log_like = gibsonLogLike(S1, Sf, Sff, sums[3], sums[4], sums[5], dN, alpha);
    }
    else
    {
      const double a = (likelihood_form == highres_form_free_alpha) ? alpha : -1.0;
      log_like = brogiLineLogLike(sums[0], sums[1], sums[2] / dN, dN, a);
    }

    atomicAdd(d_log_like, log_like);
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


// Opts a kernel into dynamic shared memory above 48 KB, once per kernel.
template <typename Kernel>
static void optInSharedMemory(Kernel kernel, const size_t shared_bytes, size_t& opted_in_bytes)
{
  if (shared_bytes > opted_in_bytes)
  {
    gpuErrchk(cudaFuncSetAttribute(
      kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(shared_bytes)));
    opted_in_bytes = shared_bytes;
  }
}


// Kernel 1a
static void launchInterp(
    const float* broadened_spectrum_dev,
    const HighResModelGridGPU& model_grid,
    const int* order_offsets_dev,
    const int* order_nb_pixels_dev,
    const double* orbital_phases_dev,
    const double* v_bary_dev,
    const float* exposure_blur_coeff_dev,
    float* model_out_dev,
    const int nb_orders,
    const int nb_exposures,
    const int max_pixels_per_order,
    const double Kp, const double Vsys, const double dphi,
    const float* model_scale_dev,
    const bool use_phase_function,
    const float* stellar_spectrum_dev)
{
  if (nb_exposures > highres_max_exposures_gpu)
    throw std::runtime_error(
      "launchInterp: " + std::to_string(nb_exposures) + " exposures exceed the GPU limit of "
      + std::to_string(highres_max_exposures_gpu) + "\n");

  static size_t opted_in_bytes = 48 * 1024;

  const int threads = 256;
  const dim3 blocks(nb_orders, (max_pixels_per_order + threads - 1) / threads);
  const size_t shared_bytes = static_cast<size_t>(nb_exposures) * sizeof(ExposureShift);

  optInSharedMemory(highResInterpKernel, shared_bytes, opted_in_bytes);

  highResInterpKernel<<<blocks, threads, shared_bytes>>>(
    broadened_spectrum_dev, model_grid, order_offsets_dev, order_nb_pixels_dev,
    orbital_phases_dev, v_bary_dev, exposure_blur_coeff_dev, model_out_dev,
    nb_orders, nb_exposures, Kp, Vsys, dphi, model_scale_dev,
    use_phase_function, stellar_spectrum_dev);

  CUDA_CHECK_AFTER_KERNEL();
}


template <int TILE, bool LIKELIHOOD_SUMS>
static int launchProjectInstance(
  const float* projection_matrices_dev,
  float* model_dev,
  const float* order_flux_dev,
  const float* flux_uncertainties_dev,
  const float* data_mean_dev,
  const float4* fit_coeff_dev,
  const int* order_offsets_dev,
  const int* order_nb_pixels_dev,
  double* partials_dev,
  const int nb_orders,
  const int nb_exposures,
  const int max_pixels_per_order,
  const int projection)
{
  static size_t opted_in_bytes = 48 * 1024;

  const size_t shared_bytes = projectSharedBytes(TILE, nb_exposures);
  optInSharedMemory(highResProjectKernel<TILE, LIKELIHOOD_SUMS>, shared_bytes, opted_in_bytes);

  const int nb_tiles = (max_pixels_per_order + TILE - 1) / TILE;
  const dim3 blocks(nb_orders, nb_tiles);
  const dim3 threads(TILE, highres_project_groups);
  const int nb_exposures_padded = ((nb_exposures + 3) / 4) * 4;

  highResProjectKernel<TILE, LIKELIHOOD_SUMS><<<blocks, threads, shared_bytes>>>(
    projection_matrices_dev, model_dev, order_flux_dev, flux_uncertainties_dev,
    data_mean_dev, fit_coeff_dev, order_offsets_dev, order_nb_pixels_dev, partials_dev,
    nb_orders, nb_exposures, nb_exposures_padded, nb_tiles, projection);

  CUDA_CHECK_AFTER_KERNEL();

  return nb_tiles;
}


// Kernel 1b with the widest tile that fits into the dynamic shared memory;
// returns the number of pixel tiles per order the partial sums are laid out for.
template <bool LIKELIHOOD_SUMS>
static int launchProject(
  const float* projection_matrices_dev,
  float* model_dev,
  const float* order_flux_dev,
  const float* flux_uncertainties_dev,
  const float* data_mean_dev,
  const float4* fit_coeff_dev,
  const int* order_offsets_dev,
  const int* order_nb_pixels_dev,
  double* partials_dev,
  const int nb_orders,
  const int nb_exposures,
  const int max_pixels_per_order,
  const int projection)
{
  const size_t max_shared = static_cast<size_t>(maxDynamicSharedMemoryBytes());

  if (projectSharedBytes(128, nb_exposures) <= max_shared)
    return launchProjectInstance<128, LIKELIHOOD_SUMS>(
      projection_matrices_dev, model_dev, order_flux_dev, flux_uncertainties_dev, data_mean_dev,
      fit_coeff_dev, order_offsets_dev, order_nb_pixels_dev, partials_dev,
      nb_orders, nb_exposures, max_pixels_per_order, projection);

  if (projectSharedBytes(64, nb_exposures) <= max_shared)
    return launchProjectInstance<64, LIKELIHOOD_SUMS>(
      projection_matrices_dev, model_dev, order_flux_dev, flux_uncertainties_dev, data_mean_dev,
      fit_coeff_dev, order_offsets_dev, order_nb_pixels_dev, partials_dev,
      nb_orders, nb_exposures, max_pixels_per_order, projection);

  if (projectSharedBytes(highres_min_tile, nb_exposures) <= max_shared)
    return launchProjectInstance<highres_min_tile, LIKELIHOOD_SUMS>(
      projection_matrices_dev, model_dev, order_flux_dev, flux_uncertainties_dev, data_mean_dev,
      fit_coeff_dev, order_offsets_dev, order_nb_pixels_dev, partials_dev,
      nb_orders, nb_exposures, max_pixels_per_order, projection);

  throw std::runtime_error(
    "launchProject: " + std::to_string(nb_exposures)
    + " exposures do not fit into the shared memory of the projection kernel\n");
}


__host__
void launchHighResLogLikeFiltered(
    const float* broadened_spectrum_dev,
    const HighResModelGridGPU& model_grid,
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
    const HighResFilteredScratchGPU& scratch,
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
  // Kernel 1a: the raw model matrix
  launchInterp(
    broadened_spectrum_dev, model_grid, order_offsets_dev, order_nb_pixels_dev,
    orbital_phases_dev, v_bary_dev, exposure_blur_coeff_dev, model_filtered_dev,
    nb_orders, nb_exposures, max_pixels_per_order, Kp, Vsys, dphi, model_scale_dev,
    use_phase_function, stellar_spectrum_dev);

  const bool gibson = likelihood_form == highres_form_gibson;
  const float* uncertainties = gibson ? flux_uncertainties_dev : nullptr;
  float4* fit_coeff_dev = reinterpret_cast<float4*>(scratch.fit_coeff);

  if (filter_model)
  {
    // the quadratic of every filtered exposure, from the raw moments
    const int warps = 8;
    const dim3 moment_blocks(nb_orders, (nb_exposures + warps - 1) / warps);

    highResRawMomentsKernel<<<moment_blocks, 32 * warps>>>(
      model_filtered_dev, order_offsets_dev, order_nb_pixels_dev, nb_exposures, scratch.raw_moments);

    CUDA_CHECK_AFTER_KERNEL();

    const int nb_tasks = nb_orders * nb_exposures;
    const int threads = 128;
    const int nb_exposures_padded = ((nb_exposures + 3) / 4) * 4;

    highResFitCoefficientsKernel<<<(nb_tasks + threads - 1) / threads, threads>>>(
      scratch.raw_moments, projection_matrices_dev, order_nb_pixels_dev,
      nb_orders, nb_exposures, nb_exposures_padded, fit_coeff_dev);

    CUDA_CHECK_AFTER_KERNEL();

    // Kernel 1b: (I - P), detrending and the sums per (order, exposure)
    const int nb_tiles = launchProject<true>(
      projection_matrices_dev, model_filtered_dev, order_flux_dev, uncertainties,
      gibson ? nullptr : data_mean_dev, fit_coeff_dev, order_offsets_dev, order_nb_pixels_dev,
      scratch.partials, nb_orders, nb_exposures, max_pixels_per_order, highres_projection_filter);

    // Kernel 2: the logarithms
    highResFilteredLogLikeKernel<<<(nb_tasks + threads - 1) / threads, threads>>>(
      scratch.partials, order_nb_pixels_dev, data_sf2_dev,
      gibson_S1_dev, gibson_Sf_dev, gibson_Sff_dev,
      nb_orders, nb_exposures, nb_tiles, likelihood_form, alpha, d_log_like_dev);

    CUDA_CHECK_AFTER_KERNEL();
  }
  else
  {
    // Data-only mode: sums over all exposures of an order, one logarithm per order
    const int nb_tiles = launchProject<true>(
      projection_matrices_dev, model_filtered_dev, order_flux_dev, uncertainties,
      nullptr, nullptr, order_offsets_dev, order_nb_pixels_dev,
      scratch.partials, nb_orders, nb_exposures, max_pixels_per_order,
      highres_projection_center_time);

    highResTotalCCFLogLikeKernel<<<1, 128>>>(
      scratch.partials, order_nb_pixels_dev,
      gibson_S1_dev, gibson_Sf_dev, gibson_Sff_dev,
      nb_orders, nb_exposures, nb_tiles, likelihood_form, alpha, d_log_like_dev);

    CUDA_CHECK_AFTER_KERNEL();
  }
}


__host__
void launchHighResInterpFilterOnly(
    const float* broadened_spectrum_dev,
    const HighResModelGridGPU& model_grid,
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
  launchInterp(
    broadened_spectrum_dev, model_grid, order_offsets_dev, order_nb_pixels_dev,
    orbital_phases_dev, v_bary_dev, exposure_blur_coeff_dev, model_filtered_dev,
    nb_orders, nb_exposures, max_pixels_per_order, Kp, Vsys, dphi, model_scale_dev,
    use_phase_function, stellar_spectrum_dev);

  if (projection == highres_projection_none)
    return;

  // the treated model, written back in place; the quadratic detrending of the
  // likelihood is not applied here
  launchProject<false>(
    projection_matrices_dev, model_filtered_dev, nullptr, nullptr, nullptr, nullptr,
    order_offsets_dev, order_nb_pixels_dev, nullptr,
    nb_orders, nb_exposures, max_pixels_per_order, projection);
}


} // namespace bear
