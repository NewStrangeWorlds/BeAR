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

//formal solution (see the comment on the second pass of the kernel): with lambda mu
//below this value the integrals of cosh(lambda t) and sinh(lambda t) along the ray are
//expanded in powers of lambda, above it their closed forms are used
constexpr float lambda_mu_series = 0.05f;

//integrals of t^n exp(-t/mu): series below this value of tau/mu, closed form above
constexpr float x_series = 1.0f;

//number of terms of the series of the lower incomplete gamma function
constexpr int nb_series_terms = 10;


//1/x in double precision from a single-precision seed and two Newton steps
//(relative error of a few 1e-16); the FP64 division sequence is several times slower
__device__ __forceinline__ double reciprocal(const double x)
{
  double r = static_cast<double>(__frcp_rn(static_cast<float>(x)));
  r = fma(r, fma(-x, r, 1.0), r);
  r = fma(r, fma(-x, r, 1.0), r);
  return r;
}



//optical properties of layer i (0 = bottom), assembled in single precision from the
//level-based gas coefficients and the layer-based cloud properties
struct LayerPropsDev {
  float tau, omega, one_minus_omega, g;
};


__device__ __forceinline__ LayerPropsDev layerPropsDev(
  const int i,
  const int tid,
  const int nb_spectral_points,
  const float* __restrict__ absorption_coeff,
  const float* __restrict__ scattering_coeff,
  const float* __restrict__ cloud_optical_depth,
  const float* __restrict__ cloud_single_scattering,
  const float* __restrict__ cloud_asym_param,
  const float* s_altitude)
{
  const int idx0 = i * nb_spectral_points + tid;
  const int idx1 = (i + 1) * nb_spectral_points + tid;

  const float dz = s_altitude[i+1] - s_altitude[i];

  const float absorption_depth = dz * (absorption_coeff[idx1] + absorption_coeff[idx0]) * 0.5f;
  const float gas_scattering_depth = dz * (scattering_coeff[idx1] + scattering_coeff[idx0]) * 0.5f;

  float tau = absorption_depth + gas_scattering_depth;
  float cloud_scattering_depth = 0.0f;
  float cloud_absorption_depth = 0.0f;
  float g_cloud = 0.0f;

  if (cloud_optical_depth != nullptr)
  {
    const float cloud_tau = cloud_optical_depth[idx0];
    tau += cloud_tau;
    cloud_absorption_depth = cloud_tau;

    if (cloud_single_scattering != nullptr)
    {
      const float cloud_ssa = cloud_single_scattering[idx0];
      cloud_scattering_depth = cloud_tau * cloud_ssa;
      cloud_absorption_depth = cloud_tau * (1.0f - cloud_ssa);
    }

    if (cloud_asym_param != nullptr)
      g_cloud = cloud_asym_param[idx0];
  }

  if (tau < 0.0f) tau = 0.0f;

  const float scattering_depth = gas_scattering_depth + cloud_scattering_depth;

  LayerPropsDev p;

  p.tau = tau;
  p.omega = (tau > 0.0f) ? fminf(scattering_depth / tau, 1.0f) : 0.0f;

  //the absorbed fraction assembled from the absorption directly: 1 - omega loses its
  //relative accuracy in single precision for nearly conservative scattering
  p.one_minus_omega = (tau > 0.0f) ? fminf(fmaxf((absorption_depth + cloud_absorption_depth) / tau, 0.0f), 1.0f) : 1.0f;

  //Rayleigh scattering has a zero asymmetry parameter
  p.g = (scattering_depth > 0.0f) ? cloud_scattering_depth / scattering_depth * g_cloud : 0.0f;

  return p;
}



//two-stream quantities of a layer (Toon et al. 1989, hemispheric mean closure)
//the two-stream system is ill-conditioned for thin and for conservatively scattering
//layers; its coefficients and its solution are therefore kept in double precision,
//with the same expressions as in the CPU version
struct TwoStreamDev {
  double tau_ts, lambda, gam, e_lambda;
  double e1, e2, e3, e4;
  double c_plus_top, c_minus_top, c_plus_bot, c_minus_bot;
};


__device__ __forceinline__ TwoStreamDev twoStreamDev(
  const LayerPropsDev& p,
  const float planck_top,
  const float planck_diff)
{
  TwoStreamDev l;

  //omega and 1 - omega as a consistent pair: the accurately known one of the two single-precision
  //values defines the other (gamma1 - gamma2 = 2 (1 - omega), and Gamma depends on both)
  double omega, one_minus_omega;

  if (p.one_minus_omega < 0.5f)
  {
    one_minus_omega = fmax(static_cast<double>(p.one_minus_omega), 1.0 - omega_max);
    omega = 1.0 - one_minus_omega;
  }
  else
  {
    omega = static_cast<double>(p.omega);
    one_minus_omega = 1.0 - omega;
  }

  const double g = p.g;

  l.tau_ts = fmax(static_cast<double>(p.tau), tau_min_two_stream);

  //the Planck gradient from the difference of the level values evaluated with full
  //relative accuracy (the two-stream solution of a diffusive layer depends on the gradient)
  const double b0 = planck_top;
  const double b1 = static_cast<double>(planck_diff) * reciprocal(l.tau_ts);

  const double gamma1 = 2.0 - omega * (1.0 + g);
  const double gamma2 = omega * (1.0 - g);
  const double gamma_sum = gamma1 + gamma2;

  //gamma1 - gamma2 = 2 (1 - omega)
  l.lambda = sqrt(2.0 * one_minus_omega * gamma_sum);
  l.gam = gamma2 * reciprocal(gamma1 + l.lambda);
  l.e_lambda = exp(-l.lambda * l.tau_ts);

  l.e1 = 1.0 + l.gam * l.e_lambda;
  l.e2 = 1.0 - l.gam * l.e_lambda;
  l.e3 = l.gam + l.e_lambda;
  l.e4 = l.gam - l.e_lambda;

  //particular solution, Eq. 27 of Toon et al. with 2 pi mu1 = pi
  const double c0 = constants::pi * b0;
  const double c1 = constants::pi * b1 * reciprocal(gamma_sum);
  const double c_bot = constants::pi * b1 * l.tau_ts;

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
    const double rb = reciprocal(b);
    c_prime[idx] = d * rb;
    d_prime[idx] = e * rb;
  }
  else
  {
    const int idx_prev = idx - nb_spectral_points;
    const double rdenom = reciprocal(b - a * c_prime[idx_prev]);
    c_prime[idx] = d * rdenom;
    d_prime[idx] = (e - a * d_prime[idx_prev]) * rdenom;
  }
}



//integrals I_n = int_0^tau t^n exp(-t/mu) dt/mu = mu^n gamma(n+1, x) for n = 1 .. n_max,
//with x = tau/mu and the lower incomplete gamma function gamma; e_mu = exp(-x)
//small x: gamma(n+1, x) = x^(n+1) e^-x sum_j x^j / ((n+1)(n+2)...(n+1+j))
//large x: gamma(n+1, x) = n! (1 - e^-x sum_(k<=n) x^k/k!)
template <int n_max>
__device__ __forceinline__ void thinLayerIntegrals(
  const float x,
  const float mu,
  const float e_mu,
  float integrals[6])
{
  float mu_n = 1.0f;

  if (x < x_series)
  {
    float x_np1 = x * x;

    #pragma unroll
    for (int n = 1; n <= n_max; ++n)
    {
      float term = 1.0f / static_cast<float>(n + 1);
      float sum = term;

      #pragma unroll
      for (int j = 1; j <= nb_series_terms; ++j)
      {
        term *= x / static_cast<float>(n + 1 + j);
        sum += term;
      }

      mu_n *= mu;
      integrals[n] = mu_n * x_np1 * e_mu * sum;
      x_np1 *= x;
    }
  }
  else
  {
    float poly = 1.0f;
    float x_k = 1.0f;
    float factorial = 1.0f;

    #pragma unroll
    for (int n = 1; n <= n_max; ++n)
    {
      x_k *= x;
      factorial *= static_cast<float>(n);
      poly += x_k / factorial;

      mu_n *= mu;
      integrals[n] = mu_n * factorial * (1.0f - e_mu * poly);
    }
  }
}

}



//one thread per spectral point; the temperature and altitude profiles are staged in shared memory
//pass 1 (top to bottom): assembles and eliminates the tridiagonal two-stream system (Toon et al. 1989, Eqs. 39-44)
//pass 2 (bottom to top): back substitution fused with the formal solution along the Gauss angles
//
//Precision: the Planck function, the layer optical properties and the formal solution are
//evaluated in single precision, the two-stream system in double precision (see twoStreamDev).
//The scattering source of a layer, G exp(-lambda (tau - t)) + H exp(-lambda t), has coefficients
//G and H that grow like 1/tau for thin layers and like 1/sqrt(1 - omega) for conservative
//scattering, while their combinations U = G e_lambda + H and V = lambda (G e_lambda - H) stay
//finite. U and V are formed in double precision and the source is integrated along the ray as
//  S(t) = B0 + B1 t + U cosh(lambda t) + V sinh(lambda t)/lambda,
//with the integrals of the last two terms expanded in powers of lambda where lambda mu is small
//and in closed form otherwise. Only for lambda mu >= 1, i.e. weakly scattering layers whose
//G and H are small, is the original form of the integrals used directly.
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

  //workspace: the eliminated tridiagonal system and the two-stream quantities needed again in pass 2
  double* c_prime = workspace;
  double* d_prime = c_prime + static_cast<size_t>(nb_rows) * nb_spectral_points;
  double* lambda_ws = d_prime + static_cast<size_t>(nb_rows) * nb_spectral_points;
  double* gam_ws = lambda_ws + static_cast<size_t>(nb_layers) * nb_spectral_points;
  double* e_lambda_ws = gam_ws + static_cast<size_t>(nb_layers) * nb_spectral_points;

  float mu[max_nb_angles];
  float weight[max_nb_angles];

  #pragma unroll
  for (int k = 0; k < max_nb_angles; ++k)
  {
    mu[k] = (k < nb_angles) ? static_cast<float>(gauss_nodes_dev[k]) : 1.0f;
    weight[k] = (k < nb_angles) ? static_cast<float>(gauss_weights_dev[k]) : 0.0f;
  }

  constexpr float two_pi = 2.0f * constants::pi;


  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    const float wavenumber = static_cast<float>(wavenumber_list_dev[tid]);
    const float wavenumber_cube = wavenumber * wavenumber * wavenumber;

    //---- pass 1: top to bottom, Toon layer n = 1 is the top layer (BeAR layer nb_layers - 1)
    float planck_top = planckFunctionAccurate(s_temperature[nb_layers], wavenumber_cube, wavenumber);
    float planck_bot = planckFunctionAccurate(s_temperature[nb_layers-1], wavenumber_cube, wavenumber);
    float planck_diff = planckDifferenceAccurate(s_temperature[nb_layers-1], s_temperature[nb_layers], wavenumber_cube, wavenumber);

    TwoStreamDev cur = twoStreamDev(
      layerPropsDev(
        nb_layers - 1, tid, nb_spectral_points,
        absorption_coeff_dev, scattering_coeff_dev,
        cloud_optical_depth_dev, cloud_single_scattering_dev, cloud_asym_param_dev,
        s_altitude),
      planck_top, planck_diff);

    //top boundary: no downward flux at the top of the atmosphere
    thomasStep(c_prime, d_prime, 0, tid, nb_spectral_points,
      0.0, cur.e1, -cur.e2, -cur.c_minus_top);

    for (int n = 2; n <= nb_layers; ++n)
    {
      const int i = nb_layers - n;

      const int idx_layer = (i + 1) * nb_spectral_points + tid;
      lambda_ws[idx_layer] = cur.lambda;
      gam_ws[idx_layer] = cur.gam;
      e_lambda_ws[idx_layer] = cur.e_lambda;

      planck_top = planck_bot;
      planck_bot = planckFunctionAccurate(s_temperature[i], wavenumber_cube, wavenumber);
      planck_diff = planckDifferenceAccurate(s_temperature[i], s_temperature[i+1], wavenumber_cube, wavenumber);

      const TwoStreamDev nxt = twoStreamDev(
        layerPropsDev(
          i, tid, nb_spectral_points,
          absorption_coeff_dev, scattering_coeff_dev,
          cloud_optical_depth_dev, cloud_single_scattering_dev, cloud_asym_param_dev,
          s_altitude),
        planck_top, planck_diff);

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

    lambda_ws[tid] = cur.lambda;
    gam_ws[tid] = cur.gam;
    e_lambda_ws[tid] = cur.e_lambda;

    //bottom boundary: black surface emitting at the bottom temperature
    //planck_bot now holds the Planck function at the bottom level
    const float planck_surface = planck_bot;

    thomasStep(c_prime, d_prime, nb_rows - 1, tid, nb_spectral_points,
      cur.e1, cur.e2, 0.0, constants::pi * static_cast<double>(planck_surface) - cur.c_plus_bot);


    //---- pass 2: bottom to top, back substitution fused with the formal solution
    float intensity[max_nb_angles];

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

      planck_top = planckFunctionAccurate(s_temperature[i+1], wavenumber_cube, wavenumber);
      planck_diff = planckDifferenceAccurate(s_temperature[i], s_temperature[i+1], wavenumber_cube, wavenumber);

      const LayerPropsDev p = layerPropsDev(
        i, tid, nb_spectral_points,
        absorption_coeff_dev, scattering_coeff_dev,
        cloud_optical_depth_dev, cloud_single_scattering_dev, cloud_asym_param_dev,
        s_altitude);

      planck_bot = planck_top;

      if (p.tau <= 0.0f) continue;

      const int idx_layer = i * nb_spectral_points + tid;
      const double lambda = lambda_ws[idx_layer];
      const double gam = gam_ws[idx_layer];
      const double e_lambda = e_lambda_ws[idx_layer];

      const bool scattering = p.tau >= static_cast<float>(tau_min_two_stream) && p.omega > 0.0f;

      const float tau_ts = fmaxf(p.tau, static_cast<float>(tau_min_two_stream));
      const float omega = fminf(p.omega, static_cast<float>(omega_max));
      const float gamma_sum = 2.0f - 2.0f * omega * p.g;

      //Planck gradients: with the floored optical depth for the two-stream source,
      //with the true one for the thermal emission
      const float b1_ts = planck_diff / tau_ts;
      const float b1 = planck_diff / p.tau;
      const float planck_shift = scattering ? b1_ts * omega * p.g / gamma_sum : 0.0f;

      //scattering source coefficients, Table 3 of Toon et al. (1989), divided by 2 pi
      const double a = y1 + y2;
      const double b = y1 - y2;
      const double G = (1.0/mu1 - lambda) * a / (2.0 * constants::pi);
      const double H = gam * (lambda + 1.0/mu1) * b / (2.0 * constants::pi);

      //source coefficients at the layer top, formed in double precision
      const double G_top = G * e_lambda;
      const float U = scattering ? static_cast<float>(G_top + H) : 0.0f;
      const float V = scattering ? static_cast<float>(lambda * (G_top - H)) : 0.0f;
      const float G_f = static_cast<float>(G);
      const float H_f = static_cast<float>(H);
      const float lambda_f = static_cast<float>(lambda);
      const float lambda2 = lambda_f * lambda_f;
      const float e_lambda_f = static_cast<float>(e_lambda);
      const float source_0 = planck_top + planck_shift;

      #pragma unroll
      for (int k = 0; k < max_nb_angles; ++k)
      {
        if (k >= nb_angles) break;

        const float x = p.tau / mu[k];
        const float e_mu = expf(-x);
        const float d_mu = -expm1f(-x);
        const float lambda_mu = lambda_f * mu[k];

        float integrals[6];

        //thermal part: B0 (1 - e_mu) + B1 int t exp(-t/mu) dt/mu
        if (scattering && lambda_mu < lambda_mu_series)
          thinLayerIntegrals<5>(x, mu[k], e_mu, integrals);
        else
          thinLayerIntegrals<1>(x, mu[k], e_mu, integrals);

        float i_up = intensity[k] * e_mu + source_0 * d_mu + b1 * integrals[1];

        if (scattering)
        {
          if (lambda_mu < lambda_mu_series)
          {
            //int cosh(lambda t) e^(-t/mu) dt/mu and int sinh(lambda t)/lambda e^(-t/mu) dt/mu as series
            const float cosh_int = d_mu + lambda2 * (integrals[2] * (1.0f/2.0f) + lambda2 * integrals[4] * (1.0f/24.0f));
            const float sinh_int = integrals[1] + lambda2 * (integrals[3] * (1.0f/6.0f) + lambda2 * integrals[5] * (1.0f/120.0f));

            i_up += U * cosh_int + V * sinh_int;
          }
          else if (lambda_mu < 1.0f)
          {
            //closed forms from the integrals of exp(+-lambda t - t/mu): x expm1(y)/y with y = (lambda mu - 1) x
            //and (1 - exp(-(lambda mu + 1) x))/(lambda mu + 1)
            const float y = (lambda_mu - 1.0f) * x;
            const float int_plus = x * ((y != 0.0f) ? expm1f(y) / y : 1.0f);
            const float int_minus = -expm1f(-(lambda_mu + 1.0f) * x) / (lambda_mu + 1.0f);

            const float cosh_int = 0.5f * (int_plus + int_minus);
            const float sinh_int = 0.5f * (int_plus - int_minus) / lambda_f;

            i_up += U * cosh_int + V * sinh_int;
          }
          else
          {
            //weakly scattering layer: G and H are small, the original integrals can be used
            //removable singularity at lambda mu = 1: with y = (lambda mu - 1) x the first integral
            //is e_mu x (1 - e^-y)/y, used where the direct difference of the exponentials cancels
            const float xi = lambda_mu - 1.0f;
            const float y = xi * x;
            const float t1 = (fabsf(y) > 0.1f)
              ? (e_mu - e_lambda_f) / xi
              : e_mu * x * ((y != 0.0f) ? -expm1f(-y) / y : 1.0f);
            const float t2 = (1.0f - e_lambda_f * e_mu) / (lambda_mu + 1.0f);

            i_up += G_f * t1 + H_f * t2;
          }
        }

        intensity[k] = i_up;
      }
    }


    float flux = 0.0f;

    #pragma unroll
    for (int k = 0; k < max_nb_angles; ++k)
      flux += intensity[k] * mu[k] * weight[k];

    model_spectrum_gpu[tid] = static_cast<float>(two_pi * flux * 1e-3 * spectrum_scaling);
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

  //workspace: two arrays of 2 nb_layers rows for the Thomas algorithm
  //and three per-layer arrays of two-stream quantities, all per spectral point
  const size_t required_size = (2 * 2 + 3) * nb_layers * nb_spectral_points;

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
