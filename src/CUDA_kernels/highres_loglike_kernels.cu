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


// Kernel 1b: temporal treatment of the raw model matrix of an order, in place.
// One block per (order, tile of TILE pixels); the tile's K x TILE values are staged
// in shared memory, laid out [exposure][pixel] so that the reads are conflict-free.
// The block has TILE x PROJECT_GROUPS threads: the groups split the exposures, and
// each thread accumulates four exposures per shared-memory read from a float4 of
// the transposed, padded projection matrix (layout: see initDeviceMemory).
//   projection_filter:      out = (I - P) raw
//   projection_center_time: subtract the temporal mean of each pixel
constexpr int highres_project_groups = 4;

template <int TILE>
__global__
void highResProjectKernel(
  const float*  __restrict__ projection_matrices_t,
  float*        __restrict__ model_inout,
  const int*    __restrict__ order_offsets,
  const int*    __restrict__ order_nb_pixels,
  const int                  nb_orders,
  const int                  nb_exposures,
  const int                  nb_exposures_padded,
  const int                  projection)
{
  const int ord = blockIdx.x;
  const int tid = threadIdx.x;
  const int group = threadIdx.y;

  if (ord >= nb_orders) return;

  const int K = nb_exposures;
  const int Kp = nb_exposures_padded;
  const int N = order_nb_pixels[ord];
  const int offset = order_offsets[ord];
  const int p = blockIdx.y * TILE + tid;

  extern __shared__ float s_tile[];   // [K][TILE]

  float* column = model_inout + (size_t)offset * K + p;

  // all threads take part in the barrier; those beyond the order's last pixel
  // only skip the loads and stores
  const bool active = p < N;

  if (active)
    for (int e = group; e < K; e += highres_project_groups)
      s_tile[e * TILE + tid] = column[(size_t)e * N];

  __syncthreads();

  if (!active) return;

  if (projection == highres_projection_filter)
  {
    const float* matrix_t = projection_matrices_t + (size_t)ord * K * Kp;

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
        const float x = s_tile[jj * TILE + tid];
        const float4 m = __ldg(mt + jj * stride4);

        acc.x += m.x * x;
        acc.y += m.y * x;
        acc.z += m.z * x;
        acc.w += m.w * x;
      }

      column[(size_t)e0 * N] = acc.x;
      if (e0 + 1 < K) column[(size_t)(e0 + 1) * N] = acc.y;
      if (e0 + 2 < K) column[(size_t)(e0 + 2) * N] = acc.z;
      if (e0 + 3 < K) column[(size_t)(e0 + 3) * N] = acc.w;
    }
  }
  else
  {
    // Data-only mode: remove only the temporal mean of each pixel (N_PCA = 0
    // model filtering).  A continuum that does not move with the planet cancels,
    // the Doppler trail of the lines is kept.
    double mean = 0.0;

    for (int e = 0; e < K; ++e)
      mean += s_tile[e * TILE + tid];

    mean /= (double)K;

    for (int e = group; e < K; e += highres_project_groups)
      column[(size_t)e * N] = (float)((double)s_tile[e * TILE + tid] - mean);
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


template <int TILE>
static void launchProjectInstance(
  const float* projection_matrices_dev,
  float* model_dev,
  const int* order_offsets_dev,
  const int* order_nb_pixels_dev,
  const int nb_orders,
  const int nb_exposures,
  const int max_pixels_per_order,
  const int projection)
{
  static size_t opted_in_bytes = 48 * 1024;

  const size_t shared_bytes = static_cast<size_t>(nb_exposures) * TILE * sizeof(float);
  optInSharedMemory(highResProjectKernel<TILE>, shared_bytes, opted_in_bytes);

  const dim3 blocks(nb_orders, (max_pixels_per_order + TILE - 1) / TILE);
  const dim3 threads(TILE, highres_project_groups);
  const int nb_exposures_padded = ((nb_exposures + 3) / 4) * 4;

  highResProjectKernel<TILE><<<blocks, threads, shared_bytes>>>(
    projection_matrices_dev, model_dev, order_offsets_dev, order_nb_pixels_dev,
    nb_orders, nb_exposures, nb_exposures_padded, projection);
}


static void launchInterpFilter(
    const float* broadened_spectrum_dev,
    const HighResModelGridGPU& model_grid,
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
  if (nb_exposures > highres_max_exposures_gpu)
    throw std::runtime_error(
      "launchInterpFilter: " + std::to_string(nb_exposures) + " exposures exceed the GPU limit of "
      + std::to_string(highres_max_exposures_gpu) + "\n");

  // Kernel 1a: the raw model matrix
  {
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

  if (projection == highres_projection_none)
    return;

  // Kernel 1b: the temporal treatment, with the widest tile whose K x TILE
  // floats fit into the dynamic shared memory
  const size_t max_shared = static_cast<size_t>(maxDynamicSharedMemoryBytes());
  const size_t bytes_per_pixel = static_cast<size_t>(nb_exposures) * sizeof(float);

  if (128 * bytes_per_pixel <= max_shared)
    launchProjectInstance<128>(projection_matrices_dev, model_out_dev, order_offsets_dev,
      order_nb_pixels_dev, nb_orders, nb_exposures, max_pixels_per_order, projection);
  else if (64 * bytes_per_pixel <= max_shared)
    launchProjectInstance<64>(projection_matrices_dev, model_out_dev, order_offsets_dev,
      order_nb_pixels_dev, nb_orders, nb_exposures, max_pixels_per_order, projection);
  else if (32 * bytes_per_pixel <= max_shared)
    launchProjectInstance<32>(projection_matrices_dev, model_out_dev, order_offsets_dev,
      order_nb_pixels_dev, nb_orders, nb_exposures, max_pixels_per_order, projection);
  else
    throw std::runtime_error(
      "launchInterpFilter: " + std::to_string(nb_exposures)
      + " exposures do not fit into the shared memory of the projection kernel\n");

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
    model_grid,
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
  launchInterpFilter(
    broadened_spectrum_dev,
    model_grid,
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
