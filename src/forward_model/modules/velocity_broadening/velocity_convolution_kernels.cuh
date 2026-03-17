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


/*
 * velocity_convolution_kernels.cuh
 *
 * This is a HEADER file for GPU code (.cuh = CUDA header).
 * It defines the MATHEMATICAL FORMULAS for the two broadening kernel shapes
 * used by convolution_kernels_highres.cu.
 *
 * WHY A SEPARATE HEADER FILE?
 * ────────────────────────────
 * The formulas here are marked __device__ __forceinline__, meaning the
 * compiler pastes them directly into any GPU kernel that calls them (like
 * an inline function in Python).  They cannot live in a .cu file because
 * GPU functions can only be called from other GPU code in the same
 * compilation unit.  Putting them in a .cuh header lets multiple .cu files
 * include and use them without duplicating code.
 *
 * WHAT IS A "KERNEL" IN THIS CONTEXT?
 * ─────────────────────────────────────
 * In signal processing, a "kernel" (also called a "filter") is the weight
 * function used in a convolution:
 *
 *   broadened_spectrum[i] = SUM_j  kernel(j - i) * spectrum[j]
 *
 * Think of it as: for each output pixel i, we take a weighted average of
 * nearby input pixels j, where the weight depends on how far j is from i.
 * A Gaussian kernel gives more weight to nearby pixels and less to distant
 * ones.  The rotational kernel has a different shape set by the physics of
 * a rotating planet.
 *
 * NOTE ON UNITS: All functions below work in PIXEL units on the log-lambda
 * grid.  Physical velocities (km/s) are converted to pixels by the caller:
 *
 *   sigma_pixels = sigma_kms / delta_v_kms
 *   vsini_pixels = vsini_kms / delta_v_kms
 *
 * where delta_v_kms = c / R_grid is the velocity spacing per pixel.
 * For HARPS observations (R_grid = 110,000): delta_v ≈ 2.73 km/s/pixel.
 */


#ifndef _velocity_convolution_kernels_cuh
#define _velocity_convolution_kernels_cuh

#pragma once

#include <cmath>


namespace bear {


/* Number of GPU threads per block used by both broadening kernels.
 *
 * 128 threads is a good balance:
 *   - Large enough to keep the GPU busy during the parallel reduction step.
 *   - Small enough that blocks can be scheduled efficiently.
 *   - Must be a multiple of 32 (the GPU's "warp" size).
 * This matches BeAR's existing low-res kernel (convolution_kernels.cu).
 */
static constexpr int HIGHRES_BLOCK_SIZE = 128;


/* ---------------------------------------------------------------------------
 * Gaussian kernel value  (instrumental broadening)
 *
 * Returns the weight that input pixel at distance dx from the centre should
 * receive when computing the Gaussian-broadened output.
 *
 * The formula is the standard Gaussian (bell curve):
 *
 *   G(dx) = (1 / (sigma * sqrt(2*pi))) * exp( -dx^2 / (2 * sigma^2) )
 *
 * For efficiency, the caller pre-computes the two constants that don't
 * change between pixels:
 *   inv_2sig2 = 1 / (2 * sigma^2)   — used in the exponent
 *   norm      = 1 / (sigma * sqrt(2*pi))  — the overall scale factor
 *
 * Physical interpretation:
 *   A spectrograph with resolving power R produces a Gaussian instrumental
 *   profile with sigma = c / (R * 2.355) in km/s.  In pixel units this is
 *   sigma_pixels = sigma_kms / delta_v_kms.
 *   For HARPS (R=110,000): sigma ≈ 1.16 km/s ≈ 0.42 pixels.
 * ---------------------------------------------------------------------------
 */
__device__ __forceinline__
float gaussianKernelHR(float dx, float inv_2sig2, float norm)
{
  return norm * expf(-dx * dx * inv_2sig2);
}


/* ---------------------------------------------------------------------------
 * Rotational broadening kernel  (Gray 2005, Eq 17.12)
 *
 * Returns the weight for a pixel at distance dx from the kernel centre,
 * for a planet rotating at vsini_pixels (in pixel units).
 *
 * PHYSICAL BACKGROUND:
 * When we observe a transiting planet at the limb, different parts of the
 * planet's terminator are moving at different line-of-sight velocities due
 * to rotation.  The flux from a point at fractional velocity x = v/vsini
 * is weighted by how much area of the disk contributes at that velocity,
 * modified by limb darkening (the centre of the disk appears brighter than
 * the edges in the stellar transit geometry).
 *
 * THE FORMULA (Gray 2005, Eq 17.12):
 *
 *   G(dx) = [ c1 * sqrt(1 - x^2) + c2 * (1 - x^2) ] / denom
 *
 * where:
 *   x     = dx / vsini_pixels            (normalised distance, -1 to +1)
 *   c1    = 2 * (1 - epsilon)            (geometric/elliptic term weight)
 *   c2    = (pi/2) * epsilon             (limb-darkening term weight)
 *   denom = pi * vsini_pixels * (1 - epsilon/3)   (normalisation)
 *   epsilon = limb-darkening coefficient (0 to 1)
 *
 * The two terms inside the brackets represent:
 *   sqrt(1 - x^2) — the geometric projection of a rotating disk (ellipse)
 *   (1 - x^2)     — the additional suppression from limb darkening
 *
 * EPSILON VALUES:
 *   epsilon = 0.0 → no limb darkening; purely elliptic profile (broader,
 *                   flatter peak).  This was the original (incorrect) GPU value.
 *   epsilon = 0.6 → standard solar-type limb darkening; narrower, more peaked
 *                   profile.  Matches model_processing.py (our default).
 *   Using the wrong epsilon changes the recovered abundances by ~1.3 dex
 *   and temperature by ~200 K in the retrieval — a significant systematic.
 *
 * COMPACT SUPPORT:
 *   The kernel is exactly zero for |dx| >= vsini_pixels (when |x| >= 1).
 *   Unlike the Gaussian which technically extends to infinity, this kernel
 *   has a hard cutoff.  So the convolution window is only ±ceil(vsini_pixels)
 *   pixels wide — for WASP-121b + HARPS that's just ±3 pixels.
 *
 * Parameters:
 *   dx            — pixel distance from the kernel centre (j - i)
 *   vsini_pixels  — vsini in pixel units (= vsini_kms / delta_v_kms)
 *   epsilon       — limb-darkening coefficient, default 0.6
 * ---------------------------------------------------------------------------
 */
__device__ __forceinline__
float rotationalKernelHR(float dx, float vsini_pixels, float epsilon)
{
  float x_norm = dx / vsini_pixels;

  if (fabsf(x_norm) >= 1.0f)
    return 0.0f;

  float one_minus_x2 = 1.0f - x_norm * x_norm;

  // Gray (2005) Eq 17.12
  float c1    = 2.0f * (1.0f - epsilon);
  float c2    = 0.5f * (float)M_PI * epsilon;
  float denom = (float)M_PI * vsini_pixels * (1.0f - epsilon / 3.0f);

  return (c1 * sqrtf(one_minus_x2) + c2 * one_minus_x2) / denom;
}


} // namespace bear


#endif
