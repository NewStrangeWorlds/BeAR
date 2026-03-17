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


#ifndef HIGHRES_LOGLIKE_KERNELS_H
#define HIGHRES_LOGLIKE_KERNELS_H


namespace bear {


// Launch the GPU kernel for Brogi & Line 2019 high-res log-likelihood.
// Computes Doppler-shifted interpolation of the broadened model onto each
// spectral order's wavelength grid, then evaluates the B&L cross-correlation
// likelihood for all orders and exposures in a single kernel launch.
//
// Result is atomically accumulated into d_log_like_dev (must be zeroed before call).
// alpha parameter: if negative, alpha is analytically marginalized (default).
// If non-negative, it is used as an explicit scaling factor in the likelihood.
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
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    float Kp, float Vsys,
    float alpha,
    double* d_log_like_dev);


// Filtered variant: Gibson et al. 2022 fast model filtering.
// Two-kernel pipeline:
//   Kernel 1 (interpFilter): Interpolate model at all Doppler shifts, apply (I-P)
//   Kernel 2 (logLikeFromFiltered): Cross-correlation likelihood from precomputed model
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
    const float* projection_matrices_dev,
    float* model_filtered_dev,
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    float Kp, float Vsys,
    float alpha,
    double* d_log_like_dev);


// Gibson et al. 2022 Eq. 4: per-pixel uncertainty weighting, beta marginalized.
// Unfiltered variant: single kernel, one block per (order, exposure).
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
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    float Kp, float Vsys,
    float alpha,
    double* d_log_like_dev);


// Gibson Eq. 4 filtered variant: uses existing interp+filter kernel (Kernel 1),
// then a Gibson-specific likelihood kernel (Kernel 2).
void launchHighResLogLikeFilteredGibson(
    const float* broadened_spectrum_dev,
    const double* model_wavelengths_dev,
    int n_model,
    const float* order_wavelengths_dev,
    const float* order_flux_dev,
    const int* order_offsets_dev,
    const int* order_nb_pixels_dev,
    const float* orbital_phases_dev,
    const float* projection_matrices_dev,
    float* model_filtered_dev,
    const float* flux_uncertainties_dev,
    const double* gibson_S1_dev,
    const double* gibson_Sf_dev,
    const double* gibson_Sff_dev,
    int nb_orders,
    int nb_exposures,
    int max_pixels_per_order,
    float Kp, float Vsys,
    float alpha,
    double* d_log_like_dev);


}


#endif
