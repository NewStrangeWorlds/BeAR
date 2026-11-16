// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-FileCopyrightText: 2026 Esther Wang
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef HIGHRES_LOGLIKE_KERNELS_H
#define HIGHRES_LOGLIKE_KERNELS_H


namespace bear {


// Largest number of exposures the filtered GPU path supports.  The projection
// kernel stages a tile of K x 32 model values in shared memory, see
// highResProjectKernel; this bound keeps the smallest tile within 64 KB.
constexpr int highres_max_exposures_gpu = 512;


// Precomputed geometry of the model grid, built once per observation by
// HighResObservation::initModelGridGPU.  All kernels locate the Doppler-shifted
// model values in float arithmetic on these tables instead of searching the
// double wavelength grid:
//  * every order owns a contiguous slice of the (descending) model grid that
//    covers its wavelength range plus a generous Doppler margin; the slice stores
//    ln(lambda) relative to a per-order reference, together with the cell widths
//    and their inverses;
//  * every pixel stores its rest position in the same units and its bracket at
//    zero velocity, which seeds the short scans at the shifted positions.
// Positions are kept as double-float pairs (hi + lo, with lo the rounding
// remainder of hi), so that the lookups reach ~1e-7 of a grid cell with FP32
// arithmetic only.
struct HighResModelGridGPU {
  const int*   slice_start = nullptr;        // [nb_orders] model index of the first slice point
  const int*   slice_len = nullptr;          // [nb_orders] points in the slice
  const int*   slice_offset = nullptr;       // [nb_orders] offset into the slice arrays
  const float* slice_ln = nullptr;           // ln(lambda) - ln_ref, descending within a slice (hi part)
  const float* slice_ln_lo = nullptr;        // low part of slice_ln
  const float* slice_inv_dln = nullptr;      // 1 / (slice_ln[j-1] - slice_ln[j]), 0 at j = 0
  const float* slice_dln = nullptr;          // slice_ln[j-1] - slice_ln[j], 0 at j = 0
  const float* order_inv_mean_dln = nullptr; // [nb_orders] 1 / mean cell width of the slice
  const float* pix_ln = nullptr;             // [total_pixels] ln(lambda_rest) - ln_ref of the order (hi part)
  const float* pix_ln_lo = nullptr;          // low part of pix_ln
  const int*   pix_j0 = nullptr;             // [total_pixels] bracket in the slice at zero velocity, -1 if outside
  const float* pix_t0 = nullptr;             // [total_pixels] weight of point j0 at zero velocity
};


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
    const float* stellar_spectrum_dev = nullptr,
    bool use_phase_function = false);


// Unfiltered path, Gibson et al. 2022 Eq. 4 with per-pixel uncertainties
// (inverse-variance weighted means subtracted per exposure).
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
    const float* stellar_spectrum_dev = nullptr,
    bool use_phase_function = false);


// Layout of the projection matrices on the device: per order the TRANSPOSED
// (I - P) with the exposure index padded to a multiple of four,
//   projection_matrices_dev[ord * K * Kp + j * Kp + e] = (I - P)[e][j],  Kp = 4 ceil(K/4),
// so that the projection kernel reads four exposures of a column as one float4.
constexpr int highres_projection_padding = 4;


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
    double* d_log_like_dev);


// Diagnostic export (not used by any likelihood): run only kernel 1 of the
// filtered path and leave the result in model_filtered_dev.
// projection = highres_projection_filter gives exactly the model the likelihood
// sees after (I - P); highres_projection_none gives the same model before
// filtering.  No likelihood is evaluated.
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
    const float* stellar_spectrum_dev);


}


#endif
