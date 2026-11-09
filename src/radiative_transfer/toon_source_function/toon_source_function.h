// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _toon_source_function_h
#define _toon_source_function_h


#include <vector>
#include <cstddef>

#include "../radiative_transfer.h"
#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../spectral_grid/spectral_grid.h"


namespace bear {


//Thermal emission with multiple scattering following Toon et al. (1989, JGR 94, 16287):
//a hemispheric-mean two-stream solution provides the scattering source function,
//the emergent intensity is then obtained from the formal solution along a set of
//Gauss-Legendre angles (source function technique) and integrated to the flux.
//In the absence of scattering the result is identical to the short characteristics
//solver with the same angular quadrature.
class ToonSourceFunction : public RadiativeTransfer{
  public:
    ToonSourceFunction(
      SpectralGrid* spectral_grid_ptr,
      const size_t nb_angles,
      const size_t nb_grid_points,
      const bool use_gpu);
    virtual ~ToonSourceFunction();

    virtual void calcSpectrum(
      const Atmosphere& atmosphere,
      const std::vector< std::vector<double> >& absorption_coeff,
      const std::vector< std::vector<double> >& scattering_coeff,
      const std::vector< std::vector<double> >& cloud_optical_depth,
      const std::vector< std::vector<double> >& cloud_single_scattering,
      const std::vector< std::vector<double> >& cloud_asym_param,
      const double spectrum_scaling,
      std::vector<double>& spectrum);

    virtual void calcSpectrumGPU(
      const Atmosphere& atmosphere,
      float* absorption_coeff_dev,
      float* scattering_coeff_dev,
      float* cloud_optical_depth,
      float* cloud_single_scattering,
      float* cloud_asym_param,
      const double spectrum_scaling,
      float* model_spectrum_dev);

    static constexpr size_t max_nb_angles = 8;

  private:
    size_t nb_angles = 2;
    size_t nb_grid_points = 0;
    bool gpu_enabled = false;

    //Gauss-Legendre quadrature on [0,1] for the angular integration of the intensity
    std::vector<double> gauss_nodes;
    std::vector<double> gauss_weights;

    //GPU data: quadrature and workspace for the tridiagonal solver
    double* gauss_nodes_dev = nullptr;
    double* gauss_weights_dev = nullptr;
    double* workspace_dev = nullptr;
    size_t workspace_size = 0;

    double calcSpectrum(
      const Atmosphere& atmosphere,
      const std::vector<double>& absorption_coeff,
      const std::vector<double>& scattering_coeff,
      const std::vector<double>& cloud_optical_depth,
      const std::vector<double>& cloud_single_scattering,
      const std::vector<double>& cloud_asym_param,
      const size_t nu_index) const;

    void initDeviceData();
};


}
#endif
