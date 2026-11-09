// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "toon_source_function.h"

#include "../layer_optical_properties.h"
#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../additional/aux_functions.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"
#include "../../spectral_grid/spectral_grid.h"


namespace bear{


namespace {

//numerical safeguards, identical to the ones used in the CUDA kernel
constexpr double omega_max = 1.0 - 1e-6;   //conservative scattering makes the two-stream solution degenerate
constexpr double tau_min_two_stream = 1e-6; //floor of the layer optical depth in the two-stream system
constexpr double mu1 = 0.5;                 //hemispheric mean, Table 1 of Toon et al. (1989)


//Gauss-Legendre nodes and weights on the interval [0,1]
void gaussLegendre01(
  const size_t n,
  std::vector<double>& nodes,
  std::vector<double>& weights)
{
  nodes.assign(n, 0.0);
  weights.assign(n, 0.0);

  for (size_t i = 0; i < n; ++i)
  {
    //initial guess for the i-th root of P_n on [-1,1]
    double x = std::cos(constants::pi * (i + 0.75) / (n + 0.5));
    double dp = 1.0;

    for (int iter = 0; iter < 100; ++iter)
    {
      double p0 = 1.0;
      double p1 = x;

      for (size_t k = 2; k <= n; ++k)
      {
        const double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
        p0 = p1;
        p1 = p2;
      }

      dp = n * (x * p1 - p0) / (x * x - 1.0);
      const double dx = p1 / dp;
      x -= dx;

      if (std::abs(dx) < 1e-15) break;
    }

    //map from [-1,1] to [0,1]; the roots come out in decreasing order
    nodes[n - 1 - i] = 0.5 * (x + 1.0);
    weights[n - 1 - i] = 1.0 / ((1.0 - x * x) * dp * dp);
  }
}


//two-stream and source function quantities of a single layer
struct ToonLayer {
  double tau = 0;        //true optical depth
  double tau_ts = 0;     //optical depth used in the two-stream system (floored)
  double omega = 0;
  double g = 0;
  double b0 = 0;         //Planck function at the layer top
  double b1 = 0;         //Planck gradient dB/dtau (with tau_ts)
  double gamma_sum = 2;  //gamma1 + gamma2
  double lambda = 2;
  double gam = 0;        //Gamma = gamma2 / (gamma1 + lambda)
  double e_lambda = 0;   //exp(-lambda tau_ts)
  double e1 = 0, e2 = 0, e3 = 0, e4 = 0;
  double c_plus_top = 0, c_minus_top = 0, c_plus_bot = 0, c_minus_bot = 0;
};


ToonLayer toonLayer(
  const double tau,
  const double omega_in,
  const double g,
  const double planck_top,
  const double planck_bot)
{
  ToonLayer l;

  l.tau = tau;
  l.tau_ts = std::max(tau, tau_min_two_stream);
  l.omega = std::min(omega_in, omega_max);
  l.g = g;

  l.b0 = planck_top;
  l.b1 = (planck_bot - planck_top) / l.tau_ts;

  //hemispheric mean closure
  const double gamma1 = 2.0 - l.omega * (1.0 + l.g);
  const double gamma2 = l.omega * (1.0 - l.g);

  l.gamma_sum = gamma1 + gamma2;
  l.lambda = std::sqrt((gamma1 - gamma2) * (gamma1 + gamma2));
  l.gam = gamma2 / (gamma1 + l.lambda);
  l.e_lambda = std::exp(-l.lambda * l.tau_ts);

  l.e1 = 1.0 + l.gam * l.e_lambda;
  l.e2 = 1.0 - l.gam * l.e_lambda;
  l.e3 = l.gam + l.e_lambda;
  l.e4 = l.gam - l.e_lambda;

  //particular solution, Eq. 27 of Toon et al. with 2 pi mu1 = pi
  const double c0 = constants::pi * l.b0;
  const double c1 = constants::pi * l.b1 / l.gamma_sum;
  const double c_bot = constants::pi * l.b1 * l.tau_ts;

  l.c_plus_top = c0 + c1;
  l.c_minus_top = c0 - c1;
  l.c_plus_bot = c0 + c_bot + c1;
  l.c_minus_bot = c0 + c_bot - c1;

  return l;
}

}



ToonSourceFunction::ToonSourceFunction(
  SpectralGrid* spectral_grid_ptr,
  const size_t nb_angles_,
  const size_t nb_grid_points_,
  const bool use_gpu)
   : RadiativeTransfer(spectral_grid_ptr)
   , nb_angles(nb_angles_)
   , nb_grid_points(nb_grid_points_)
   , gpu_enabled(use_gpu)
{
  if (nb_angles < 1 || nb_angles > max_nb_angles)
  {
    std::string error_message = "Toon source function RT requires between 1 and "
      + std::to_string(max_nb_angles) + " angles\n";
    throw InvalidInput(std::string ("ToonSourceFunction::ToonSourceFunction"), error_message);
  }

  gaussLegendre01(nb_angles, gauss_nodes, gauss_weights);

  if (gpu_enabled) initDeviceData();
}



void ToonSourceFunction::calcSpectrum(
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



//solves the two-stream problem and the formal solution for a single spectral point
//returns the upward flux at the top of the atmosphere in W m-2 cm-1
double ToonSourceFunction::calcSpectrum(
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

  std::vector<double> planck(nb_grid_points, 0.0);

  for (size_t i=0; i<nb_grid_points; ++i)
    planck[i] = aux::planckFunctionWavenumber(atmosphere.temperature[i], wavenumber);


  //layer quantities, BeAR layer index: 0 = bottom
  std::vector<ToonLayer> layers(nb_layers);

  for (size_t i=0; i<nb_layers; ++i)
  {
    const LayerOpticalProperties p = layerOpticalProperties(
      atmosphere.altitude,
      absorption_coeff,
      scattering_coeff,
      cloud_optical_depth,
      cloud_single_scattering,
      cloud_asym_param,
      i);

    layers[i] = toonLayer(
      p.optical_depth,
      p.single_scattering,
      p.cloud_scattering_fraction * p.cloud_asym_param,
      planck[i+1],
      planck[i]);
  }


  //two-stream solution: tridiagonal system of Toon et al. (1989), Eqs. 39-44
  //rows are ordered from the top of the atmosphere (Toon layer n = 1) to the bottom
  //Toon layer n corresponds to BeAR layer nb_layers - n
  const size_t nb_rows = 2 * nb_layers;
  std::vector<double> c_prime(nb_rows, 0.0);
  std::vector<double> d_prime(nb_rows, 0.0);

  auto thomasStep = [&](const size_t r, const double a, const double b, const double d, const double e)
  {
    if (r == 0)
    {
      c_prime[r] = d / b;
      d_prime[r] = e / b;
    }
    else
    {
      const double denom = b - a * c_prime[r-1];
      c_prime[r] = d / denom;
      d_prime[r] = (e - a * d_prime[r-1]) / denom;
    }
  };

  {
    //top boundary: no downward flux at the top of the atmosphere
    const ToonLayer& top = layers[nb_layers-1];
    thomasStep(0, 0.0, top.e1, -top.e2, -top.c_minus_top);

    //interfaces between Toon layers n and n+1 (rows 2n and 2n+1)
    for (size_t n = 1; n < nb_layers; ++n)
    {
      const ToonLayer& cur = layers[nb_layers - n];
      const ToonLayer& nxt = layers[nb_layers - n - 1];

      const double delta_plus = nxt.c_plus_top - cur.c_plus_bot;
      const double delta_minus = nxt.c_minus_top - cur.c_minus_bot;

      //even row 2n, 0-based index 2n-1
      thomasStep(2*n - 1,
        nxt.e2 * cur.e1 - nxt.e4 * cur.e3,
        nxt.e2 * cur.e2 - nxt.e4 * cur.e4,
        nxt.e4 * nxt.e1 - nxt.e2 * nxt.e3,
        nxt.e2 * delta_plus - nxt.e4 * delta_minus);

      //odd row 2n+1, 0-based index 2n
      thomasStep(2*n,
        cur.e3 * cur.e2 - cur.e1 * cur.e4,
        cur.e1 * nxt.e1 - cur.e3 * nxt.e3,
        cur.e3 * nxt.e4 - cur.e1 * nxt.e2,
        cur.e3 * delta_plus - cur.e1 * delta_minus);
    }

    //bottom boundary: black surface emitting at the bottom temperature
    const ToonLayer& bot = layers[0];
    thomasStep(nb_rows - 1, bot.e1, bot.e2, 0.0, constants::pi * planck[0] - bot.c_plus_bot);
  }


  //back substitution fused with the upward formal solution (bottom to top)
  //Toon layer n uses the rows 2n-1 and 2n (1-based), i.e. 0-based 2n-2 and 2n-1
  std::vector<double> intensity(nb_angles, planck[0]);

  double y_next = 0.0;

  for (size_t i = 0; i < nb_layers; ++i)
  {
    const size_t n = nb_layers - i;   //Toon layer index
    const size_t r2 = 2*n - 1;        //0-based row of Y_2n
    const size_t r1 = 2*n - 2;        //0-based row of Y_2n-1

    const double y2 = (r2 == nb_rows - 1) ? d_prime[r2] : d_prime[r2] - c_prime[r2] * y_next;
    const double y1 = d_prime[r1] - c_prime[r1] * y2;
    y_next = y1;

    const ToonLayer& l = layers[i];

    if (l.tau <= 0.0) continue;

    //coefficients of the exponential terms in the upward and downward two-stream fluxes
    const double a = y1 + y2;
    const double b = y1 - y2;

    const bool scattering = l.tau >= tau_min_two_stream && l.omega > 0.0;

    //scattering source coefficients, Table 3 of Toon et al. (1989), divided by 2 pi
    const double G = (1.0/mu1 - l.lambda) * a / (2.0 * constants::pi);
    const double H = l.gam * (l.lambda + 1.0/mu1) * b / (2.0 * constants::pi);
    const double planck_shift = scattering ? l.b1 * l.omega * l.g / l.gamma_sum : 0.0;

    //Planck gradient of the true layer for the thermal part of the source
    const double b1 = (planck[i] - planck[i+1]) / l.tau;

    for (size_t k = 0; k < nb_angles; ++k)
    {
      const double mu = gauss_nodes[k];
      const double x = l.tau / mu;
      const double e_mu = std::exp(-x);

      //integral of tau exp(-tau/mu) dtau/mu, series for thin layers
      const double lin_term = (x > 1e-3)
        ? -mu * std::expm1(-x) - l.tau * e_mu
        : l.tau * x / 2.0 * (1.0 - 2.0 * x / 3.0);

      double i_up = intensity[k] * e_mu
        + (l.b0 + planck_shift) * (1.0 - e_mu)
        + b1 * lin_term;

      if (scattering)
      {
        //removable singularity at lambda mu = 1
        const double xi = l.lambda * mu - 1.0;
        const double t1 = (std::abs(xi) > 1e-6)
          ? (e_mu - l.e_lambda) / xi
          : e_mu * x * (1.0 - xi * x / 2.0);
        const double t2 = (1.0 - l.e_lambda * e_mu) / (l.lambda * mu + 1.0);

        i_up += G * t1 + H * t2;
      }

      intensity[k] = i_up;
    }
  }


  double flux = 0.0;

  for (size_t k = 0; k < nb_angles; ++k)
    flux += intensity[k] * gauss_nodes[k] * gauss_weights[k];

  return 2.0 * constants::pi * flux * 1e-3; //in W m-2 cm-1
}


}
