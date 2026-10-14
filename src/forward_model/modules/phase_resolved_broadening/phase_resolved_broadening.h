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


#ifndef _phase_resolved_broadening_h
#define _phase_resolved_broadening_h

#include <vector>
#include <string>
#include <memory>

#include "../module.h"

#include "../../../spectral_grid/spectral_grid.h"
#include "../../../CUDA_kernels/data_management_kernels.h"

namespace bear {


// Phase-resolved broadening module based on Brogi et al. (2016) Section 4.2.
// Constructs a 2D model of the planet's atmospheric ring transiting a
// limb-darkened stellar disc.  Each ring pixel contributes a Gaussian
// (instrumental) weighted by the stellar intensity at that position.
// The resulting 1D broadening kernel is then convolved with the spectrum.
//
// Parameters (7):
//   0: v_eq      [km/s]  equatorial rotation velocity
//   1: v_wind    [km/s]  equatorial super-rotation wind
//   2: sigma_inst[km/s]  instrumental broadening width
//   3: u1        [-]     quadratic stellar limb-darkening coeff 1
//   4: u2        [-]     quadratic stellar limb-darkening coeff 2
//   5: Rp_Rs     [-]     planet-to-star radius ratio
//   6: impact_b  [-]     transit impact parameter
class PhaseResolvedBroadening : public Module{
  public:
    PhaseResolvedBroadening(
      const std::vector<std::string>& parameters,
      SpectralGrid* spectral_grid_);
    virtual ~PhaseResolvedBroadening();

    virtual void modifySpectrum(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      std::vector<double>& spectrum);

    virtual void modifySpectrumGPU(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      float* spectrum_gpu);

    void setSpectralGrid(SpectralGrid* grid);

  private:
    SpectralGrid* spectral_grid;

    // Velocity step between neighbouring grid points (km/s), and its minimum.
    // The high-res grid is not uniform in ln(lambda), so the kernel is evaluated
    // at the actual velocity offsets instead of on a fixed pixel grid.
    std::vector<double> dv_kms;
    double dv_min_kms = 0;
    float* dv_gpu = nullptr;

    float* temp_buffer_gpu = nullptr;
    float* table_gpu = nullptr;
    size_t table_gpu_capacity = 0;

    static constexpr int N_ANGLE = 720;
    static constexpr double WIND_LAT_LIMIT_DEG = 25.0;

    void setVelocitySteps();
    void freeDeviceBuffers();

    // Broadening kernel tabulated at velocities v_min + k*step, k = 0..size-1,
    // non-zero within |v| <= half_width.
    void buildBroadeningKernel(
      double v_eq, double v_wind, double sigma_inst,
      double u1, double u2, double Rp_Rs, double impact_b,
      std::vector<double>& table, double& v_min, double& step, double& half_width);
};


}

#endif
