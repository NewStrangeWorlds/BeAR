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


#ifndef HIGHRES_OBSERVATION_H
#define HIGHRES_OBSERVATION_H

#include <vector>
#include <string>

#include "../CUDA_kernels/data_management_kernels.h"


namespace bear {


// Likelihood mode for high-resolution cross-correlation spectroscopy.
// marginalized_alpha: -N/2 ln(1 - r^2) (Zucker 2003); sigma marginalized and alpha
//                     set to its maximum-likelihood value in each order and exposure
// free_alpha:         Brogi & Line 2019 with a sampled alpha, sigma marginalized
// gibson:             Gibson et al. 2022 Eq. 4, per-pixel uncertainties, beta marginalized
enum class HighResLikelihoodMode { marginalized_alpha, free_alpha, gibson };


class HighResObservation {
  public:
    struct SpectralOrder {
      std::vector<double> wavelengths;
      std::vector<std::vector<double>> flux;  // flux[exposure][pixel]
      std::vector<std::vector<double>> flux_uncertainties;  // [exposure][pixel], empty if not present
      size_t nb_pixels = 0;

      // GPU pointers (per-order, used by CPU fallback path)
      float* wavelengths_gpu = nullptr;
      float* flux_gpu = nullptr;  // flattened [nb_exposures * nb_pixels]
    };

    HighResObservation() {}
    ~HighResObservation();

    void init(const std::string& file_path);
    void initDeviceMemory();
    void freeDeviceMemory();

    const std::string& observationName() const { return observation_name; }
    double resolvingPower() const { return resolving_power; }
    size_t nbExposures() const { return nb_exposures; }
    size_t nbOrders() const { return nb_orders; }
    const std::vector<double>& orbitalPhases() const { return orbital_phases; }
    const std::vector<SpectralOrder>& orders() const { return spectral_orders; }

    // Overall wavelength range for spectral grid construction
    double wavelengthMin() const { return wavelength_min; }
    double wavelengthMax() const { return wavelength_max; }

    bool hasFiltering() const { return has_filtering; }
    bool hasFluxUncertainties() const { return has_flux_uncertainties; }

    // Diagnostic export (not part of any likelihood). The interpolated model
    // for every order and exposure, either before or after the (I-P) temporal
    // filter, computed with the same GPU kernel the likelihood uses.
    // Layout matches the flux buffer: [order_offset*nb_exp + exp*N + pixel].
    std::vector<float> modelMatrixGPU(
      const float* broadened_spectrum_gpu,
      const double* model_wavelengths_gpu,
      size_t nb_model_points,
      double Kp, double Vsys, double dphi,
      bool apply_projection,
      const float* stellar_spectrum_gpu = nullptr) const;

    size_t totalPixels() const { return total_pixels; }
    double kpRef() const;
    double vsysRef() const;
    bool filterModel() const { return filter_model; }
    bool reinjectModel() const { return reinject_model; }

    // (I-P)-filtered observed data, [order][exposure][pixel]
    const std::vector<std::vector<std::vector<double>>>& filteredFlux() const {
      return filtered_flux; }

    void precomputeGibsonStatistics();

    // Likelihood mode: marginalized_alpha (default), free_alpha, or gibson
    HighResLikelihoodMode likelihood_mode = HighResLikelihoodMode::marginalized_alpha;

    // Compute high-res log-likelihood (CPU)
    // In marginalized_alpha mode, the alpha parameter is ignored.
    // stellar_spectrum / nb_stellar_points: optional cached Fs (same wavelength grid as
    // model_wavelengths, units µm).  When non-null, a per-pixel correction factor
    // Fs(λ·v_dop)/Fs(λ_rest) is applied to undo the erroneous Fs Doppler shift.
    double computeLogLikelihood(
      const std::vector<double>& broadened_spectrum,
      const std::vector<double>& model_wavelengths,
      double Kp, double Vsys, double dphi, double alpha = 1.0,
      const double* stellar_spectrum = nullptr,
      size_t nb_stellar_points = 0) const;

    // Compute high-res log-likelihood (GPU)
    // Accumulates result into d_log_like_dev via atomicAdd.
    // stellar_spectrum_gpu: optional cached Fs on device (same grid as model), µm.
    void computeLogLikelihoodGPU(
      const float* broadened_spectrum_gpu,
      const double* model_wavelengths_gpu,
      size_t nb_model_points,
      double Kp, double Vsys, double dphi, double alpha,
      double* d_log_like_dev,
      const float* stellar_spectrum_gpu = nullptr) const;

  private:
    std::string observation_name;
    double resolving_power = 0;
    std::vector<double> orbital_phases;
    std::vector<SpectralOrder> spectral_orders;
    size_t nb_exposures = 0;
    size_t nb_orders = 0;

    // Reference velocities: retrieved Kp/Vsys are offsets from these values.
    // When absent from the data file both default to 0 (direct retrieval, old behaviour).
    double kp_ref = 0.0;
    double vsys_ref = 0.0;

    // Per-exposure barycentric velocity corrections (km/s).
    // If absent from the data file, all values default to 0 (i.e. user applied
    // barycentric correction in pre-processing).
    std::vector<double> barycentric_velocities;

    // --- Exposure (velocity) blurring ---
    // During a finite exposure the planet's line-of-sight velocity drifts, so the
    // observed spectrum is the time-average of the instantaneous Doppler-shifted
    // spectrum: a top-hat (boxcar) in velocity of full width
    //   dV(phi) = (kp_ref+Kp) * cos(2*pi*(phi+dphi)) * (2*pi/P) * t_exp   [km/s]
    // maximal at conjunction, zero at quadrature.  Applied as a running-mean of the
    // model via a precomputed cumulative integral (O(1) per pixel, width-independent).
    // Enabled when both orbital_period and exposure_times are present in the data file.
    bool exposure_blurring = false;
    double orbital_period = 0.0;             // orbital period in seconds
    std::vector<double> exposure_times;      // per-exposure integration time in seconds

    double wavelength_min = 0;
    double wavelength_max = 0;

    void loadDataFile(const std::string& file_path);

    // --- Per-pixel uncertainties (Gibson Eq. 4) ---
    bool has_flux_uncertainties = false;

    // Precomputed data-only weighted sums for Gibson likelihood [nb_orders * nb_exposures]
    std::vector<double> gibson_S1;   // sum 1/sigma_i^2
    std::vector<double> gibson_Sf;   // sum f_i/sigma_i^2
    std::vector<double> gibson_Sff;  // sum f_i^2/sigma_i^2

    // --- Model filtering (Gibson et al. 2022) ---
    bool has_filtering = false;
    // Data-only mode when false (#filter_model 0): (I-P) is applied to the data only;
    // the model is centred in time per pixel instead, and the cross-correlation sums
    // are accumulated over all exposures of an order before the logarithm.  The model
    // then does not undergo the distortion the filter imposes on the data, so this
    // mode is a diagnostic rather than an unbiased likelihood.
    bool filter_model = true;
    // CHIMERA-style re-injection (Line et al. 2021): if true, model_scale is computed
    // in initFiltering() as S = P*raw_flux (the background the filter captures).  The
    // model is multiplied by S before the (I-P) projection, which puts the planet
    // signal in detector units.  Requires filtering with filter_model = true, and is
    // only informative with a free alpha.
    bool reinject_model = false;
    // Dayside phase function [0.5*(1+cos(2*pi*phase-pi))]^2 (Pelletier et al. 2025)
    // applied per exposure to the model before filtering.  Enabled by
    // #phase_function 1 in the observation data file.
    bool use_phase_function = false;
    size_t nb_basis_vectors = 0;

    // Per-order (I - P) projection matrix, row-major [nb_exp x nb_exp]
    std::vector<std::vector<double>> projection_matrices;  // [nb_orders][nb_exp*nb_exp]

    // Per-order basis vectors U_raw (nb_exp x nb_basis), loaded from file
    std::vector<std::vector<double>> basis_vectors_raw;  // [nb_orders][nb_exp*nb_basis]

    // Per-order mean uncertainties for weighting (nb_exposures per order)
    std::vector<std::vector<double>> uncertainties;  // [nb_orders][nb_exp]

    // Filtered observed data and precomputed statistics
    std::vector<std::vector<std::vector<double>>> filtered_flux;  // [ord][exp][pix]
    std::vector<double> filtered_data_mean;   // [nb_orders * nb_exposures]
    std::vector<double> filtered_data_sf2;    // [nb_orders * nb_exposures]

    void loadFilteringBasis(const std::string& file_path);
    void loadModelScale(const std::string& file_path);
    void initFiltering();
    void applyProjection(
      size_t ord,
      std::vector<std::vector<double>>& matrix) const;

    // Interpolate model onto an order's wavelength grid
    void interpolateModelOntoOrder(
      const SpectralOrder& order,
      const std::vector<double>& broadened_spectrum,
      const std::vector<double>& model_wavelengths,
      double doppler_inv,
      std::vector<double>& model_on_order) const;

    // Interpolate model onto an order's grid with exposure (boxcar) blurring.
    // model_cumulative_integral is the cumulative trapezoidal integral of
    // broadened_spectrum over model_wavelengths (built once per likelihood call).
    // The value at each pixel is the mean of the model over the velocity interval
    // [v_rad - |delta_v|/2, v_rad + |delta_v|/2], evaluated as the difference of the
    // cumulative integral at the two box edges divided by their wavelength span.
    void interpolateModelOntoOrderBlurred(
      const SpectralOrder& order,
      const std::vector<double>& broadened_spectrum,
      const std::vector<double>& model_cumulative_integral,
      const std::vector<double>& model_wavelengths,
      double v_rad,
      double delta_v,
      std::vector<double>& model_on_order) const;

    // --- GPU buffers ---
    // Flattened GPU buffers for the unfiltered kernel
    double* all_wavelengths_dev = nullptr;   // flattened order wavelengths (nm)
    float* all_flux_dev = nullptr;          // flattened flux [per order: nb_exp * nb_pix]
    int* order_offsets_dev = nullptr;       // pixel offset per order
    int* order_nb_pixels_dev = nullptr;     // pixels per order
    double* orbital_phases_dev = nullptr;
    double* barycentric_velocities_dev = nullptr;
    // Per-exposure exposure-blur coefficient (2*pi/P) * t_exp [dimensionless];
    // delta_v = Kp * cos(2*pi*(phi+dphi)) * coeff.  Null when blurring is disabled.
    float* exposure_blur_coeff_dev = nullptr;
    float* data_mean_dev = nullptr;         // [nb_orders * nb_exposures]
    double* data_sf2_dev = nullptr;         // [nb_orders * nb_exposures]
    int max_pixels_per_order = 0;

    // Additional GPU buffers for filtered path
    float* projection_matrices_dev = nullptr;  // [nb_orders * nb_exp * nb_exp]
    float* model_filtered_dev = nullptr;       // workspace [total_pixels * nb_exposures]
    size_t total_pixels = 0;

    // Optional per-pixel, per-exposure scale matrix for model re-injection
    // (CHIMERA-style: the model is multiplied by this matrix before the (I-P)
    // projection).  Either S = P*raw_flux (#reinject_model 1) or read from a
    // #model_scale section.  Layout: same as all_flux_dev
    // [order_offset * nb_exposures + exp * N + pixel].  Null when absent.
    bool has_model_scale = false;
    std::vector<float> model_scale_host;  // flattened host copy
    float* model_scale_dev = nullptr;

    // GPU buffers for Gibson likelihood (per-pixel uncertainties)
    float* flux_uncertainties_dev = nullptr;   // same layout as all_flux_dev
    double* gibson_S1_dev = nullptr;           // [nb_orders * nb_exposures]
    double* gibson_Sf_dev = nullptr;
    double* gibson_Sff_dev = nullptr;
};


}

#endif
