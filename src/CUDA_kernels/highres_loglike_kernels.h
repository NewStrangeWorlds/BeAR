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


// Largest number of exposures the filtered GPU path supports.  The per-pixel
// model column is held in a fixed-size local array; see highResInterpFilterKernel.
constexpr int highres_max_exposures_gpu = 512;


// Likelihood form, mirrored from HighResLikelihoodMode so that this header does
// not depend on the observation class.
constexpr int highres_form_marginalized_alpha = 0;
constexpr int highres_form_free_alpha = 1;
constexpr int highres_form_gibson = 2;


// What happens to the interpolated model matrix of an order before the sums.
constexpr int highres_projection_none = 0;         // raw model (diagnostics only)
constexpr int highres_projection_filter = 1;       // (I - P) M
constexpr int highres_projection_center_time = 2;  // subtract the temporal mean per pixel


// Unfiltered path (no #filtering_basis): Brogi & Line 2019 likelihood with the
// per-exposure means subtracted from data and model.  One block per (order,
// exposure); the result is atomically added to d_log_like_dev.
// alpha < 0: alpha maximised analytically (marginalized-alpha mode);
// alpha >= 0: explicit scale factor (free-alpha mode).
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
    const float* stellar_spectrum_dev = nullptr,
    bool use_phase_function = false);


// Unfiltered path, Gibson et al. 2022 Eq. 4 with per-pixel uncertainties
// (inverse-variance weighted means subtracted per exposure).
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
    const float* stellar_spectrum_dev = nullptr,
    bool use_phase_function = false);


// Filtered path (a #filtering_basis is present), for all three likelihood forms.
//
// Kernel 1 interpolates the model at every exposure, applies the optional
// phase weight and re-injection scale, and then either
//   filter_model = true:  applies (I - P)                    (Gibson et al. 2022), or
//   filter_model = false: subtracts the temporal mean per pixel (data-only mode).
// With filter_model = true, the likelihood is evaluated per (order, exposure)
// after removing a quadratic in pixel index from the model (the data were
// detrended the same way at initialisation).  With filter_model = false, the
// sums are accumulated over all exposures of an order before the logarithm.
//
// model_scale_dev: optional re-injection matrix (same layout as order_flux_dev),
// only valid together with filter_model = true.
// The Gibson inputs (flux_uncertainties_dev, gibson_*_dev) are only read for
// likelihood_form == highres_form_gibson.
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
    double* d_log_like_dev);


// Diagnostic export (not used by any likelihood): run only kernel 1 of the
// filtered path and leave the result in model_filtered_dev.
// projection = highres_projection_filter gives exactly the model the likelihood
// sees after (I - P); highres_projection_none gives the same model before
// filtering.  No likelihood is evaluated.
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
    const float* stellar_spectrum_dev);


}


#endif
