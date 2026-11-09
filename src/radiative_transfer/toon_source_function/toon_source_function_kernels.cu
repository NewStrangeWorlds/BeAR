// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <vector>
#include <cmath>

#include "toon_source_function.h"

#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../additional/physical_const.h"
#include "../../CUDA_kernels/error_check.h"
#include "../../CUDA_kernels/planck_function.h"
#include "../../CUDA_kernels/data_management_kernels.h"


namespace bear{


namespace {

//numerical safeguards, identical to the CPU version
constexpr double omega_max = 1.0 - 1e-6;
constexpr double tau_min_two_stream = 1e-6;
constexpr double mu1 = 0.5;


//two-stream and source function quantities of a single layer
struct ToonLayerDev {
  double tau, tau_ts, omega, g, b0, b1;
  double gamma_sum, lambda, gam, e_lambda;
  double e1, e2, e3, e4;
  double c_plus_top, c_minus_top, c_plus_bot, c_minus_bot;
};


//assembles the optical properties of layer i (0 = bottom) and the two-stream quantities
//the computation is done in double precision: the two-stream system becomes
//ill-conditioned for thin or conservatively scattering layers
__device__ __forceinline__ ToonLayerDev toonLayerDev(
  const int i,
  const int tid,
  const int nb_spectral_points,
  const float* __restrict__ absorption_coeff,
  const float* __restrict__ scattering_coeff,
  const float* __restrict__ cloud_optical_depth,
  const float* __restrict__ cloud_single_scattering,
  const float* __restrict__ cloud_asym_param,
  const float* s_altitude,
  const double planck_top,
  const double planck_bot)
{
  ToonLayerDev l;

  const int idx0 = i * nb_spectral_points + tid;
  const int idx1 = (i + 1) * nb_spectral_points + tid;

  const double dz = s_altitude[i+1] - s_altitude[i];

  const double absorption_depth = dz * (absorption_coeff[idx1] + absorption_coeff[idx0]) * 0.5;
  const double gas_scattering_depth = dz * (scattering_coeff[idx1] + scattering_coeff[idx0]) * 0.5;

  double tau = absorption_depth + gas_scattering_depth;
  double cloud_scattering_depth = 0.0;
  double g_cloud = 0.0;

  if (cloud_optical_depth != nullptr)
  {
    const double cloud_tau = cloud_optical_depth[idx0];
    tau += cloud_tau;

    if (cloud_single_scattering != nullptr)
      cloud_scattering_depth = cloud_tau * cloud_single_scattering[idx0];

    if (cloud_asym_param != nullptr)
      g_cloud = cloud_asym_param[idx0];
  }

  if (tau < 0.0) tau = 0.0;

  const double scattering_depth = gas_scattering_depth + cloud_scattering_depth;

  double omega = (tau > 0.0) ? scattering_depth / tau : 0.0;
  if (omega > omega_max) omega = omega_max;

  //Rayleigh scattering has a zero asymmetry parameter
  const double g = (scattering_depth > 0.0) ? cloud_scattering_depth / scattering_depth * g_cloud : 0.0;

  l.tau = tau;
  l.tau_ts = fmax(tau, tau_min_two_stream);
  l.omega = omega;
  l.g = g;

  l.b0 = planck_top;
  l.b1 = (planck_bot - planck_top) / l.tau_ts;

  //hemispheric mean closure
  const double gamma1 = 2.0 - omega * (1.0 + g);
  const double gamma2 = omega * (1.0 - g);

  l.gamma_sum = gamma1 + gamma2;
  l.lambda = sqrt((gamma1 - gamma2) * (gamma1 + gamma2));
  l.gam = gamma2 / (gamma1 + l.lambda);
  l.e_lambda = exp(-l.lambda * l.tau_ts);

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



//forward elimination step of the Thomas algorithm for row r
//the modified coefficients are stored in the workspace, indexed [row][spectral point]
__device__ __forceinline__ void thomasStep(
  double* __restrict__ c_prime,
  double* __restrict__ d_prime,
  const int r,
  const int tid,
  const int nb_spectral_points,
  const double a,
  const double b,
  const double d,
  const double e)
{
  const int idx = r * nb_spectral_points + tid;

  if (r == 0)
  {
    c_prime[idx] = d / b;
    d_prime[idx] = e / b;
  }
  else
  {
    const int idx_prev = idx - nb_spectral_points;
    const double denom = b - a * c_prime[idx_prev];
    c_prime[idx] = d / denom;
    d_prime[idx] = (e - a * d_prime[idx_prev]) / denom;
  }
}

}



//one thread per spectral point; the temperature and altitude profiles are staged in shared memory
//pass 1 (top to bottom): assembles and eliminates the tridiagonal two-stream system (Toon et al. 1989, Eqs. 39-44)
//pass 2 (bottom to top): back substitution fused with the formal solution along the Gauss angles
__global__
void toonSourceFunctionDev(
  float* __restrict__ model_spectrum_gpu,
  const float* __restrict__ absorption_coeff_dev,
  const float* __restrict__ scattering_coeff_dev,
  const float* __restrict__ cloud_optical_depth_dev,
  const float* __restrict__ cloud_single_scattering_dev,
  const float* __restrict__ cloud_asym_param_dev,
  const double* __restrict__ wavenumber_list_dev,
  const float* __restrict__ temperature_dev,
  const float* __restrict__ altitude_dev,
  const double* __restrict__ gauss_nodes_dev,
  const double* __restrict__ gauss_weights_dev,
  const int nb_angles,
  double* __restrict__ workspace,
  const double spectrum_scaling,
  const int nb_spectral_points,
  const int nb_grid_points)
{
  constexpr int max_nb_angles = ToonSourceFunction::max_nb_angles;

  extern __shared__ float shared_data[];
  float* s_altitude = shared_data;
  float* s_temperature = &shared_data[nb_grid_points];

  for (int i = threadIdx.x; i < nb_grid_points; i += blockDim.x)
  {
    s_altitude[i] = altitude_dev[i];
    s_temperature[i] = temperature_dev[i];
  }

  __syncthreads();

  const int nb_layers = nb_grid_points - 1;
  const int nb_rows = 2 * nb_layers;

  double* c_prime = workspace;
  double* d_prime = workspace + static_cast<size_t>(nb_rows) * nb_spectral_points;

  double mu[max_nb_angles];
  double weight[max_nb_angles];

  #pragma unroll
  for (int k = 0; k < max_nb_angles; ++k)
  {
    mu[k] = (k < nb_angles) ? gauss_nodes_dev[k] : 1.0;
    weight[k] = (k < nb_angles) ? gauss_weights_dev[k] : 0.0;
  }


  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    const double wavenumber = wavenumber_list_dev[tid];

    //---- pass 1: top to bottom, Toon layer n = 1 is the top layer (BeAR layer nb_layers - 1)
    double planck_top = planckFunction(static_cast<double>(s_temperature[nb_layers]), wavenumber);
    double planck_bot = planckFunction(static_cast<double>(s_temperature[nb_layers-1]), wavenumber);

    ToonLayerDev cur = toonLayerDev(
      nb_layers - 1, tid, nb_spectral_points,
      absorption_coeff_dev, scattering_coeff_dev,
      cloud_optical_depth_dev, cloud_single_scattering_dev, cloud_asym_param_dev,
      s_altitude, planck_top, planck_bot);

    //top boundary: no downward flux at the top of the atmosphere
    thomasStep(c_prime, d_prime, 0, tid, nb_spectral_points,
      0.0, cur.e1, -cur.e2, -cur.c_minus_top);

    for (int n = 2; n <= nb_layers; ++n)
    {
      const int i = nb_layers - n;

      planck_top = planck_bot;
      planck_bot = planckFunction(static_cast<double>(s_temperature[i]), wavenumber);

      const ToonLayerDev nxt = toonLayerDev(
        i, tid, nb_spectral_points,
        absorption_coeff_dev, scattering_coeff_dev,
        cloud_optical_depth_dev, cloud_single_scattering_dev, cloud_asym_param_dev,
        s_altitude, planck_top, planck_bot);

      const double delta_plus = nxt.c_plus_top - cur.c_plus_bot;
      const double delta_minus = nxt.c_minus_top - cur.c_minus_bot;

      //interface between Toon layers n-1 and n: even row 2(n-1) and odd row 2(n-1)+1 (1-based)
      thomasStep(c_prime, d_prime, 2*n - 3, tid, nb_spectral_points,
        nxt.e2 * cur.e1 - nxt.e4 * cur.e3,
        nxt.e2 * cur.e2 - nxt.e4 * cur.e4,
        nxt.e4 * nxt.e1 - nxt.e2 * nxt.e3,
        nxt.e2 * delta_plus - nxt.e4 * delta_minus);

      thomasStep(c_prime, d_prime, 2*n - 2, tid, nb_spectral_points,
        cur.e3 * cur.e2 - cur.e1 * cur.e4,
        cur.e1 * nxt.e1 - cur.e3 * nxt.e3,
        cur.e3 * nxt.e4 - cur.e1 * nxt.e2,
        cur.e3 * delta_plus - cur.e1 * delta_minus);

      cur = nxt;
    }

    //bottom boundary: black surface emitting at the bottom temperature
    //planck_bot now holds the Planck function at the bottom level
    const double planck_surface = planck_bot;

    thomasStep(c_prime, d_prime, nb_rows - 1, tid, nb_spectral_points,
      cur.e1, cur.e2, 0.0, constants::pi * planck_surface - cur.c_plus_bot);


    //---- pass 2: bottom to top, back substitution fused with the formal solution
    double intensity[max_nb_angles];

    #pragma unroll
    for (int k = 0; k < max_nb_angles; ++k)
      intensity[k] = planck_surface;

    double y_next = 0.0;

    for (int i = 0; i < nb_layers; ++i)
    {
      const int n = nb_layers - i;   //Toon layer index
      const int r2 = 2*n - 1;        //0-based row of Y_2n
      const int r1 = 2*n - 2;        //0-based row of Y_2n-1

      const int idx2 = r2 * nb_spectral_points + tid;
      const int idx1 = r1 * nb_spectral_points + tid;

      const double y2 = (r2 == nb_rows - 1) ? d_prime[idx2] : d_prime[idx2] - c_prime[idx2] * y_next;
      const double y1 = d_prime[idx1] - c_prime[idx1] * y2;
      y_next = y1;

      planck_top = planckFunction(static_cast<double>(s_temperature[i+1]), wavenumber);

      const ToonLayerDev l = toonLayerDev(
        i, tid, nb_spectral_points,
        absorption_coeff_dev, scattering_coeff_dev,
        cloud_optical_depth_dev, cloud_single_scattering_dev, cloud_asym_param_dev,
        s_altitude, planck_top, planck_bot);

      const double planck_layer_bot = planck_bot;
      planck_bot = planck_top;

      if (l.tau <= 0.0) continue;

      const double a = y1 + y2;
      const double b = y1 - y2;

      const bool scattering = l.tau >= tau_min_two_stream && l.omega > 0.0;

      //scattering source coefficients, Table 3 of Toon et al. (1989), divided by 2 pi
      const double G = (1.0/mu1 - l.lambda) * a / (2.0 * constants::pi);
      const double H = l.gam * (l.lambda + 1.0/mu1) * b / (2.0 * constants::pi);
      const double planck_shift = scattering ? l.b1 * l.omega * l.g / l.gamma_sum : 0.0;

      //Planck gradient of the true layer for the thermal part of the source
      const double b1 = (planck_layer_bot - planck_top) / l.tau;

      #pragma unroll
      for (int k = 0; k < max_nb_angles; ++k)
      {
        if (k >= nb_angles) break;

        const double x = l.tau / mu[k];
        const double e_mu = exp(-x);

        //integral of tau exp(-tau/mu) dtau/mu, series for thin layers
        const double lin_term = (x > 1e-3)
          ? -mu[k] * expm1(-x) - l.tau * e_mu
          : l.tau * x / 2.0 * (1.0 - 2.0 * x / 3.0);

        double i_up = intensity[k] * e_mu
          + (l.b0 + planck_shift) * (1.0 - e_mu)
          + b1 * lin_term;

        if (scattering)
        {
          //removable singularity at lambda mu = 1
          const double xi = l.lambda * mu[k] - 1.0;
          const double t1 = (fabs(xi) > 1e-6)
            ? (e_mu - l.e_lambda) / xi
            : e_mu * x * (1.0 - xi * x / 2.0);
          const double t2 = (1.0 - l.e_lambda * e_mu) / (l.lambda * mu[k] + 1.0);

          i_up += G * t1 + H * t2;
        }

        intensity[k] = i_up;
      }
    }


    double flux = 0.0;

    #pragma unroll
    for (int k = 0; k < max_nb_angles; ++k)
      flux += intensity[k] * mu[k] * weight[k];

    model_spectrum_gpu[tid] = static_cast<float>(2.0 * constants::pi * flux * 1e-3 * spectrum_scaling);
  }
}



ToonSourceFunction::~ToonSourceFunction()
{
  if (gauss_nodes_dev != nullptr) deleteFromDevice(gauss_nodes_dev);
  if (gauss_weights_dev != nullptr) deleteFromDevice(gauss_weights_dev);
  if (workspace_dev != nullptr) deleteFromDevice(workspace_dev);
}



void ToonSourceFunction::initDeviceData()
{
  allocateOnDevice(gauss_nodes_dev, nb_angles);
  allocateOnDevice(gauss_weights_dev, nb_angles);

  moveToDevice(gauss_nodes_dev, gauss_nodes);
  moveToDevice(gauss_weights_dev, gauss_weights);
}



void ToonSourceFunction::calcSpectrumGPU(
  const Atmosphere& atmosphere,
  float* absorption_coeff_dev,
  float* scattering_coeff_dev,
  float* cloud_optical_depth_dev,
  float* cloud_single_scattering_dev,
  float* cloud_asym_param_dev,
  const double spectrum_scaling,
  float* model_spectrum_dev)
{
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  const size_t nb_layers = nb_grid_points - 1;

  //workspace for the Thomas algorithm: two arrays of 2 nb_layers rows per spectral point
  const size_t required_size = 2 * 2 * nb_layers * nb_spectral_points;

  if (workspace_dev == nullptr || workspace_size != required_size)
  {
    if (workspace_dev != nullptr) deleteFromDevice(workspace_dev);

    allocateOnDevice(workspace_dev, required_size);
    workspace_size = required_size;
  }

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  toonSourceFunctionDev<<<blocks, threads, 2 * nb_grid_points * sizeof(float)>>>(
    model_spectrum_dev,
    absorption_coeff_dev,
    scattering_coeff_dev,
    cloud_optical_depth_dev,
    cloud_single_scattering_dev,
    cloud_asym_param_dev,
    spectral_grid->wavenumber_list_gpu,
    atmosphere.temperature_dev,
    atmosphere.altitude_dev,
    gauss_nodes_dev,
    gauss_weights_dev,
    static_cast<int>(nb_angles),
    workspace_dev,
    spectrum_scaling,
    static_cast<int>(nb_spectral_points),
    static_cast<int>(nb_grid_points));

  CUDA_CHECK_AFTER_KERNEL();
}


}
