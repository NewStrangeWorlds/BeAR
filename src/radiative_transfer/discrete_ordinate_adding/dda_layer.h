// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _dda_layer_h
#define _dda_layer_h

#include <cmath>


// Included by both the CPU solver (host compiler) and the CUDA kernel (nvcc),
// so that both paths use exactly the same four-stream layer algebra.
#ifdef __CUDACC__
  #define BEAR_HOST_DEVICE __host__ __device__
#else
  #define BEAR_HOST_DEVICE
#endif


namespace bear {
namespace dda {


//Four-stream discrete ordinate adding method for thermal emission
//(Zhang et al. 2013, JQSRT 122, 64; Exo-FMS lw_DDA).
//Each homogeneous layer is solved analytically for its 2x2 reflection and
//transmission matrices and its thermal sources; the layers are then combined
//with the adding method from the bottom up.
//Double-Gauss quadrature: Gauss-Legendre nodes on [0,1] in each hemisphere,
//the same angles as used by the short characteristics solver.

constexpr double mu0 = 0.2113248654051871;
constexpr double mu1 = 0.7886751345948129;
constexpr double weight0 = 0.5;
constexpr double weight1 = 0.5;

constexpr int nb_moments = 4;                //Legendre moments used by the solver (l = 0..3)
constexpr double omega_max = 1.0 - 1e-6;     //conservative scattering makes the eigen problem degenerate
constexpr double tau_min = 1e-10;            //thinner layers are treated as transparent
constexpr double tiny = 1e-300;


struct Mat2 { double a00, a01, a10, a11; };
struct Vec2 { double v0, v1; };


BEAR_HOST_DEVICE inline Mat2 mul(const Mat2& A, const Mat2& B)
{
  return {A.a00*B.a00 + A.a01*B.a10, A.a00*B.a01 + A.a01*B.a11,
          A.a10*B.a00 + A.a11*B.a10, A.a10*B.a01 + A.a11*B.a11};
}

BEAR_HOST_DEVICE inline Vec2 mul(const Mat2& A, const Vec2& v)
{
  return {A.a00*v.v0 + A.a01*v.v1, A.a10*v.v0 + A.a11*v.v1};
}

BEAR_HOST_DEVICE inline Mat2 add(const Mat2& A, const Mat2& B)
{
  return {A.a00 + B.a00, A.a01 + B.a01, A.a10 + B.a10, A.a11 + B.a11};
}

BEAR_HOST_DEVICE inline Mat2 sub(const Mat2& A, const Mat2& B)
{
  return {A.a00 - B.a00, A.a01 - B.a01, A.a10 - B.a10, A.a11 - B.a11};
}

BEAR_HOST_DEVICE inline Vec2 add(const Vec2& a, const Vec2& b)
{
  return {a.v0 + b.v0, a.v1 + b.v1};
}

BEAR_HOST_DEVICE inline Mat2 inverse(const Mat2& A)
{
  const double rdet = 1.0 / (A.a00*A.a11 - A.a01*A.a10);
  return {A.a11*rdet, -A.a01*rdet, -A.a10*rdet, A.a00*rdet};
}



//Legendre polynomials P_0 ... P_3 at x
BEAR_HOST_DEVICE inline void legendre(const double x, double p[nb_moments])
{
  p[0] = 1.0;
  p[1] = x;
  p[2] = 0.5 * (3.0*x*x - 1.0);
  p[3] = 0.5 * (5.0*x*x*x - 3.0*x);
}



//reflection and transmission matrices of a homogeneous layer, equal from both sides,
//and its thermal sources: Eu is the intensity emitted upward at the layer top,
//Ed the one emitted downward at the layer bottom (per stream)
struct Layer {
  Mat2 R, T;
  Vec2 Eu, Ed;
};



//Phase function moments of a layer with Rayleigh scattering by the gas and
//Henyey-Greenstein scattering by the particles, mixed by the particle fraction
//of the scattering optical depth. With delta-M scaling the truncated moment
//chi_4 is removed from the phase function, and tau and omega are scaled accordingly.
BEAR_HOST_DEVICE inline void scaledOpticalProperties(
  double& tau,
  double& omega,
  const double particle_fraction,
  const double g,
  const bool use_delta_m,
  double chi[nb_moments])
{
  if (omega > omega_max) omega = omega_max;

  const double gas_fraction = 1.0 - particle_fraction;
  const double g2 = g*g;

  //Rayleigh: chi_0 = 1, chi_2 = 0.1; Henyey-Greenstein: chi_l = g^l
  chi[0] = 1.0;
  chi[1] = particle_fraction * g;
  chi[2] = gas_fraction * 0.1 + particle_fraction * g2;
  chi[3] = particle_fraction * g2 * g;

  const double chi4 = particle_fraction * g2 * g2;

  const bool scale = use_delta_m && omega > 0.0 && chi4 > 1e-12 && chi4 < 1.0 - 1e-12;

  if (!scale) return;

  const double f = chi4;
  const double rf = 1.0 / (1.0 - f);

  for (int l = 1; l < nb_moments; ++l)
    chi[l] = (chi[l] - f) * rf;

  tau *= 1.0 - omega * f;
  omega = omega * (1.0 - f) / (1.0 - omega * f);
}



//decaying mode I+ = gp exp(-k tau), I- = gm exp(-k tau) of the four-stream
//system for the eigenvalue k^2 of Q = (alpha + beta)(alpha - beta)
BEAR_HOST_DEVICE inline void decayingMode(
  const Mat2& Q,
  const Mat2& amb,
  const double k,
  Vec2& gp,
  Vec2& gm)
{
  const double lam = k*k;

  //eigenvector of Q from the row with the larger residual norm
  const Vec2 v1 = {Q.a01, lam - Q.a00};
  const Vec2 v2 = {lam - Q.a11, Q.a10};

  const double n1 = v1.v0*v1.v0 + v1.v1*v1.v1;
  const double n2 = v2.v0*v2.v0 + v2.v1*v2.v1;

  Vec2 X;

  if (n1 >= n2)
  {
    const double s = 1.0 / sqrt(n1 > tiny ? n1 : tiny);
    X = {v1.v0*s, v1.v1*s};
  }
  else
  {
    const double s = 1.0 / sqrt(n2 > tiny ? n2 : tiny);
    X = {v2.v0*s, v2.v1*s};
  }

  const Vec2 Y = mul(amb, X);

  gp = {0.5*(X.v0 - Y.v0/k), 0.5*(X.v1 - Y.v1/k)};
  gm = {0.5*(X.v0 + Y.v0/k), 0.5*(X.v1 + Y.v1/k)};
}



//analytic four-stream solution of a single homogeneous layer with
//(scaled) optical depth tau, single scattering albedo omega, phase function moments chi
//and a Planck function varying linearly in tau between planck_top and planck_bot
BEAR_HOST_DEVICE inline Layer layerSolution(
  const double tau,
  const double omega,
  const double chi[nb_moments],
  const double planck_top,
  const double planck_bot)
{
  Layer l;

  if (tau <= tau_min)
  {
    l.R = {0.0, 0.0, 0.0, 0.0};
    l.T = {1.0, 0.0, 0.0, 1.0};
    l.Eu = {0.0, 0.0};
    l.Ed = {0.0, 0.0};

    return l;
  }

  //phase matrices Ws_ij = w_j/2 sum_l chi_l (2l+1) P_l(mu_i) P_l(mu_j) for scattering
  //within a hemisphere and Wd with P_l(-mu_j) = (-1)^l P_l(mu_j) across the hemispheres
  double p0[nb_moments], p1[nb_moments];
  legendre(mu0, p0);
  legendre(mu1, p1);

  Mat2 Ws = {0.0, 0.0, 0.0, 0.0};
  Mat2 Wd = {0.0, 0.0, 0.0, 0.0};

  for (int m = 0; m < nb_moments; ++m)
  {
    const double c = chi[m] * (2.0*m + 1.0) * 0.5 * weight0; //both weights are equal
    const double s = (m % 2 == 0) ? 1.0 : -1.0;

    const double c00 = c * p0[m] * p0[m];
    const double c01 = c * p0[m] * p1[m];
    const double c11 = c * p1[m] * p1[m];

    Ws.a00 += c00; Ws.a01 += c01; Ws.a10 += c01; Ws.a11 += c11;
    Wd.a00 += s*c00; Wd.a01 += s*c01; Wd.a10 += s*c01; Wd.a11 += s*c11;
  }

  const double r0 = 1.0 / mu0;
  const double r1 = 1.0 / mu1;

  //dI+/dtau = alpha I+ - beta I-, dI-/dtau = beta I+ - alpha I- (homogeneous part)
  const Mat2 al = {(1.0 - omega*Ws.a00)*r0, -omega*Ws.a01*r0, -omega*Ws.a10*r1, (1.0 - omega*Ws.a11)*r1};
  const Mat2 bt = {omega*Wd.a00*r0, omega*Wd.a01*r0, omega*Wd.a10*r1, omega*Wd.a11*r1};

  const Mat2 amb = sub(al, bt);
  const Mat2 apb = add(al, bt);
  const Mat2 Q = mul(apb, amb);

  //eigenvalues k^2 of Q; the smaller root via det/lambda1 to avoid cancellation
  const double tr = Q.a00 + Q.a11;
  const double det = Q.a00*Q.a11 - Q.a01*Q.a10;
  const double disc = tr*tr - 4.0*det;
  const double lam1 = 0.5 * (tr + sqrt(disc > tiny ? disc : tiny));
  const double k1 = sqrt(lam1);
  const double k2 = sqrt(det/lam1 > tiny ? det/lam1 : tiny);

  Vec2 gp1, gm1, gp2, gm2;
  decayingMode(Q, amb, k1, gp1, gm1);
  decayingMode(Q, amb, k2, gp2, gm2);

  //columns are the modes
  const Mat2 Gp = {gp1.v0, gp2.v0, gp1.v1, gp2.v1};
  const Mat2 Gm = {gm1.v0, gm2.v0, gm1.v1, gm2.v1};

  const double e1 = exp(-k1 * tau);
  const double e2 = exp(-k2 * tau);

  const Mat2 GpE = {Gp.a00*e1, Gp.a01*e2, Gp.a10*e1, Gp.a11*e2};
  const Mat2 GmE = {Gm.a00*e1, Gm.a01*e2, Gm.a10*e1, Gm.a11*e2};

  //R + T = (Gp + Gm E)(Gm + Gp E)^-1, R - T = (Gp - Gm E)(Gm - Gp E)^-1
  const Mat2 P1 = mul(add(Gp, GmE), inverse(add(Gm, GpE)));
  const Mat2 P2 = mul(sub(Gp, GmE), inverse(sub(Gm, GpE)));

  l.R = {0.5*(P1.a00 + P2.a00), 0.5*(P1.a01 + P2.a01), 0.5*(P1.a10 + P2.a10), 0.5*(P1.a11 + P2.a11)};
  l.T = {0.5*(P1.a00 - P2.a00), 0.5*(P1.a01 - P2.a01), 0.5*(P1.a10 - P2.a10), 0.5*(P1.a11 - P2.a11)};

  //thermal particular solution for B linear in tau: I+/- = B(tau) +/- d B1
  //with (1 - omega (Ws - Wd)) d = mu
  const double b1 = (planck_bot - planck_top) / tau;

  const Mat2 Wm = sub(Ws, Wd);
  const Mat2 M = {1.0 - omega*Wm.a00, -omega*Wm.a01, -omega*Wm.a10, 1.0 - omega*Wm.a11};
  const Vec2 mu_vec = {mu0, mu1};
  const Vec2 d = mul(inverse(M), mu_vec);
  const Vec2 P2d = mul(P2, d);

  const double sum_b = planck_top + planck_bot;
  const double dif_b = planck_top - planck_bot;

  const double sum_o0 = sum_b * (1.0 - P1.a00 - P1.a01);
  const double sum_o1 = sum_b * (1.0 - P1.a10 - P1.a11);
  const double dif_o0 = dif_b * (1.0 - P2.a00 - P2.a01) + 2.0 * b1 * (d.v0 + P2d.v0);
  const double dif_o1 = dif_b * (1.0 - P2.a10 - P2.a11) + 2.0 * b1 * (d.v1 + P2d.v1);

  l.Eu = {0.5*(sum_o0 + dif_o0), 0.5*(sum_o1 + dif_o1)};
  l.Ed = {0.5*(sum_o0 - dif_o0), 0.5*(sum_o1 - dif_o1)};

  return l;
}



//adds a layer on top of the slab below it, described by its reflection matrix Rb
//and the upward intensity Su it emits at its top; both are updated in place
BEAR_HOST_DEVICE inline void addLayer(Mat2& Rb, Vec2& Su, const Layer& l)
{
  const Mat2 RR = mul(Rb, l.R);
  const Mat2 one_minus_RR = {1.0 - RR.a00, -RR.a01, -RR.a10, 1.0 - RR.a11};
  const Mat2 A = mul(l.T, inverse(one_minus_RR));

  const Vec2 src = add(l.Eu, mul(A, mul(Rb, l.Ed)));
  const Vec2 Su_new = add(src, mul(A, Su));

  Rb = add(l.R, mul(A, mul(Rb, l.T)));
  Su = Su_new;
}



//upward flux at the top from the upward intensities per stream, F = 2 pi sum_i w_i mu_i I_i
BEAR_HOST_DEVICE inline double upwardFlux(const Vec2& I)
{
  return 2.0 * 3.14159265358979323846 * (weight0*mu0*I.v0 + weight1*mu1*I.v1);
}


}
}

#endif
