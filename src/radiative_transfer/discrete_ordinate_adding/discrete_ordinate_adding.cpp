// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <vector>
#include <cmath>
#include <omp.h>

#include "discrete_ordinate_adding.h"
#include "dda_layer.h"

#include "../layer_optical_properties.h"
#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../additional/aux_functions.h"
#include "../../spectral_grid/spectral_grid.h"


namespace bear{


DiscreteOrdinateAdding::DiscreteOrdinateAdding(
  SpectralGrid* spectral_grid_ptr,
  const size_t nb_grid_points_,
  const bool use_gpu,
  const bool use_delta_m_)
   : RadiativeTransfer(spectral_grid_ptr)
   , nb_grid_points(nb_grid_points_)
   , gpu_enabled(use_gpu)
   , use_delta_m(use_delta_m_)
{

}



void DiscreteOrdinateAdding::calcSpectrum(
  const Atmosphere& atmosphere,
  const std::vector< std::vector<double> >& absorption_coeff,
  const std::vector< std::vector<double> >& scattering_coeff,
  const std::vector< std::vector<double> >& cloud_optical_depth,
  const std::vector< std::vector<double> >& cloud_single_scattering,
  const std::vector< std::vector<double> >& cloud_asym_param,
  const double spectrum_scaling,
  std::vector<double>& spectrum)
{
  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<spectrum.size(); ++i)
    spectrum[i] = calcSpectrum(
      atmosphere,
      absorption_coeff[i],
      scattering_coeff[i],
      cloud_optical_depth[i],
      cloud_single_scattering[i],
      cloud_asym_param[i],
      i) * spectrum_scaling;
}



//adds the layers from the bottom up for a single spectral point
//returns the upward flux at the top of the atmosphere in W m-2 cm-1
double DiscreteOrdinateAdding::calcSpectrum(
  const Atmosphere& atmosphere,
  const std::vector<double>& absorption_coeff,
  const std::vector<double>& scattering_coeff,
  const std::vector<double>& cloud_optical_depth,
  const std::vector<double>& cloud_single_scattering,
  const std::vector<double>& cloud_asym_param,
  const size_t nu_index) const
{
  const size_t nb_layers = nb_grid_points - 1;
  const double wavenumber = spectral_grid->wavenumber_list[nu_index];

  //lower boundary: black surface emitting at the bottom temperature
  double planck_bot = aux::planckFunctionWavenumber(atmosphere.temperature[0], wavenumber);

  dda::Mat2 reflection = {0.0, 0.0, 0.0, 0.0};
  dda::Vec2 intensity = {planck_bot, planck_bot};

  for (size_t i=0; i<nb_layers; ++i)
  {
    const double planck_top = aux::planckFunctionWavenumber(atmosphere.temperature[i+1], wavenumber);

    const LayerOpticalProperties p = layerOpticalProperties(
      atmosphere.altitude,
      absorption_coeff,
      scattering_coeff,
      cloud_optical_depth,
      cloud_single_scattering,
      cloud_asym_param,
      i);

    double tau = p.optical_depth;
    double omega = p.single_scattering;
    double chi[dda::nb_moments];

    dda::scaledOpticalProperties(tau, omega, p.cloud_scattering_fraction, p.cloud_asym_param, use_delta_m, chi);

    const dda::Layer layer = dda::layerSolution(tau, omega, chi, planck_top, planck_bot);

    dda::addLayer(reflection, intensity, layer);

    planck_bot = planck_top;
  }

  return dda::upwardFlux(intensity) * 1e-3; //in W m-2 cm-1
}


}
