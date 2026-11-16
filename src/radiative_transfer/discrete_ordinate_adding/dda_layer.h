// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _dda_layer_h
#define _dda_layer_h

#include <cmath>


// Included by both the CPU solver (host compiler) and the CUDA kernel (nvcc),
// so that both paths use exactly the same four-stream layer algebra.
// The algebra is templated on the floating-point type: the CPU solver uses
// double, the GPU kernel float. All expressions are arranged so that no
// quantity of order tau or (1 - omega) is obtained as a difference of
// quantities of order unity, which keeps the float version accurate to a few
// 1e-6 for thin layers and for (nearly) conservative scattering.
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
constexpr double one_minus_omega_min = 1e-6; //the corresponding floor of the absorbed fraction
constexpr double tau_min = 1e-10;            //thinner layers are treated as transparent
constexpr double tau_thin = 0.1;             //k1 tau below this: product form of the thermal source


//precision-dependent constants and math functions
template <typename T> struct Traits;

template <> struct Traits<double> {
  BEAR_HOST_DEVICE static double tiny() { return 1e-300; }
  BEAR_HOST_DEVICE static double sqrt(const double x) { return ::sqrt(x); }
  BEAR_HOST_DEVICE static double expm1(const double x) { return ::expm1(x); }
};

template <> struct Traits<float> {
  BEAR_HOST_DEVICE static float tiny() { return 1e-30f; }
  BEAR_HOST_DEVICE static float sqrt(const float x) { return ::sqrtf(x); }
  BEAR_HOST_DEVICE static float expm1(const float x) { return ::expm1f(x); }
};



template <typename T> struct Mat2 { T a00, a01, a10, a11; };
template <typename T> struct Vec2 { T v0, v1; };


template <typename T>
BEAR_HOST_DEVICE inline Mat2<T> mul(const Mat2<T>& A, const Mat2<T>& B)
{
  return {A.a00*B.a00 + A.a01*B.a10, A.a00*B.a01 + A.a01*B.a11,
          A.a10*B.a00 + A.a11*B.a10, A.a10*B.a01 + A.a11*B.a11};
}

template <typename T>
BEAR_HOST_DEVICE inline Vec2<T> mul(const Mat2<T>& A, const Vec2<T>& v)
{
  return {A.a00*v.v0 + A.a01*v.v1, A.a10*v.v0 + A.a11*v.v1};
}

template <typename T>
BEAR_HOST_DEVICE inline Mat2<T> add(const Mat2<T>& A, const Mat2<T>& B)
{
  return {A.a00 + B.a00, A.a01 + B.a01, A.a10 + B.a10, A.a11 + B.a11};
}

template <typename T>
BEAR_HOST_DEVICE inline Mat2<T> sub(const Mat2<T>& A, const Mat2<T>& B)
{
  return {A.a00 - B.a00, A.a01 - B.a01, A.a10 - B.a10, A.a11 - B.a11};
}

template <typename T>
BEAR_HOST_DEVICE inline Vec2<T> add(const Vec2<T>& a, const Vec2<T>& b)
{
  return {a.v0 + b.v0, a.v1 + b.v1};
}

template <typename T>
BEAR_HOST_DEVICE inline Mat2<T> scale(const Mat2<T>& A, const T s)
{
  return {A.a00*s, A.a01*s, A.a10*s, A.a11*s};
}

//scales the columns of A: (A diag(s))_ij = A_ij s_j
template <typename T>
BEAR_HOST_DEVICE inline Mat2<T> scaleColumns(const Mat2<T>& A, const Vec2<T>& s)
{
  return {A.a00*s.v0, A.a01*s.v1, A.a10*s.v0, A.a11*s.v1};
}

//scales the rows of A: (diag(s) A)_ij = s_i A_ij
template <typename T>
BEAR_HOST_DEVICE inline Mat2<T> scaleRows(const Vec2<T>& s, const Mat2<T>& A)
{
  return {A.a00*s.v0, A.a01*s.v0, A.a10*s.v1, A.a11*s.v1};
}

template <typename T>
BEAR_HOST_DEVICE inline Vec2<T> rowSums(const Mat2<T>& A)
{
  return {A.a00 + A.a01, A.a10 + A.a11};
}

template <typename T>
BEAR_HOST_DEVICE inline Mat2<T> inverse(const Mat2<T>& A)
{
  const T rdet = T(1) / (A.a00*A.a11 - A.a01*A.a10);
  return {A.a11*rdet, -A.a01*rdet, -A.a10*rdet, A.a00*rdet};
}

template <typename T>
BEAR_HOST_DEVICE inline Mat2<T> identity()
{
  return {T(1), T(0), T(0), T(1)};
}



//Legendre polynomials P_0 ... P_3 at x
template <typename T>
BEAR_HOST_DEVICE inline void legendre(const T x, T p[nb_moments])
{
  p[0] = T(1);
  p[1] = x;
  p[2] = T(0.5) * (T(3)*x*x - T(1));
  p[3] = T(0.5) * (T(5)*x*x*x - T(3)*x);
}



//reflection and transmission matrices of a homogeneous layer, equal from both sides,
//and its thermal sources: Eu is the intensity emitted upward at the layer top,
//Ed the one emitted downward at the layer bottom (per stream)
template <typename T>
struct Layer {
  Mat2<T> R, T_;
  Vec2<T> Eu, Ed;
};



//Phase function moments of a layer with Rayleigh scattering by the gas and
//Henyey-Greenstein scattering by the particles, mixed by the particle fraction
//of the scattering optical depth. With delta-M scaling the truncated moment
//chi_4 is removed from the phase function, and tau, omega and 1 - omega are
//scaled accordingly.
//one_minus_omega is the absorbed fraction of the extinction, assembled from the
//absorption directly rather than as 1 - omega, so that it stays accurate for
//nearly conservative scattering.
template <typename T>
BEAR_HOST_DEVICE inline void scaledOpticalProperties(
  T& tau,
  T& omega,
  T& one_minus_omega,
  const T particle_fraction,
  const T g,
  const bool use_delta_m,
  T chi[nb_moments])
{
  if (omega > T(omega_max)) omega = T(omega_max);
  if (one_minus_omega < T(one_minus_omega_min)) one_minus_omega = T(one_minus_omega_min);

  const T gas_fraction = T(1) - particle_fraction;
  const T g2 = g*g;

  //Rayleigh: chi_0 = 1, chi_2 = 0.1; Henyey-Greenstein: chi_l = g^l
  chi[0] = T(1);
  chi[1] = particle_fraction * g;
  chi[2] = gas_fraction * T(0.1) + particle_fraction * g2;
  chi[3] = particle_fraction * g2 * g;

  const T chi4 = particle_fraction * g2 * g2;

  const bool scale = use_delta_m && omega > T(0) && chi4 > T(1e-12) && chi4 < T(1) - T(1e-12);

  if (!scale) return;

  const T f = chi4;
  const T rf = T(1) / (T(1) - f);
  const T rof = T(1) / (T(1) - omega * f);

  for (int l = 1; l < nb_moments; ++l)
    chi[l] = (chi[l] - f) * rf;

  tau *= T(1) - omega * f;
  omega = omega * (T(1) - f) * rof;
  one_minus_omega = one_minus_omega * rof;
}



//eigenvector X of Q for the eigenvalue lam, taken from the row with the larger
//residual norm, and the corresponding Y = (alpha - beta) X = lam (alpha + beta)^-1 X.
//The second form avoids the cancellation in (alpha - beta) X for the slowly
//decaying mode of a nearly conservative layer.
template <typename T>
BEAR_HOST_DEVICE inline void eigenMode(
  const Mat2<T>& Q,
  const Mat2<T>& apb_inv,
  const T lam,
  Vec2<T>& X,
  Vec2<T>& Y)
{
  const T tiny = Traits<T>::tiny();

  const Vec2<T> v1 = {Q.a01, lam - Q.a00};
  const Vec2<T> v2 = {lam - Q.a11, Q.a10};

  const T n1 = v1.v0*v1.v0 + v1.v1*v1.v1;
  const T n2 = v2.v0*v2.v0 + v2.v1*v2.v1;

  if (n1 >= n2)
  {
    const T s = T(1) / Traits<T>::sqrt(n1 > tiny ? n1 : tiny);
    X = {v1.v0*s, v1.v1*s};
  }
  else
  {
    const T s = T(1) / Traits<T>::sqrt(n2 > tiny ? n2 : tiny);
    X = {v2.v0*s, v2.v1*s};
  }

  const Vec2<T> w = mul(apb_inv, X);
  Y = {lam * w.v0, lam * w.v1};
}



//x - (1 - exp(-x)) for x = k tau and d = 1 - exp(-x), accurate for small arguments
template <typename T>
BEAR_HOST_DEVICE inline T expm1Remainder(const T x, const T d)
{
  if (x < T(0.5))
    return x*x * (T(1./2.) + x*(T(-1./6.) + x*(T(1./24.) + x*(T(-1./120.) + x*(T(1./720.)
         + x*(T(-1./5040.) + x*(T(1./40320.) - x*T(1./362880.))))))));

  return x - d;
}



//analytic four-stream solution of a single homogeneous layer with
//(scaled) optical depth tau, single scattering albedo omega (and absorbed fraction
//one_minus_omega), phase function moments chi and a Planck function varying
//linearly in tau between planck_top and planck_bot.
//
//With the modes X e^{-k tau} of the system and Y = (alpha - beta) X, the
//matrices Gp = (X - Y/k)/2 and Gm = (X + Y/k)/2 of Zhang et al. give
//  R + T = (Gp + Gm E)(Gm + Gp E)^-1 = I - A1,   A1 = Delta D (X - Gp D)^-1,
//  R - T = (Gp - Gm E)(Gm - Gp E)^-1 = A2 - I,   A2 = X D (Delta + Gp D)^-1,
//with Delta = Gm - Gp = Y/k, E = exp(-k tau) and D = I - E; A1 and A2 are of
//order tau and are evaluated from products of small quantities only.
//For thin layers A2 is written as tau (alpha + beta) + A2', where A2' is of
//order tau^2, since the thermal source needs A2 d - tau 1 to that order.
template <typename T>
BEAR_HOST_DEVICE inline Layer<T> layerSolution(
  const T tau,
  const T omega,
  const T one_minus_omega,
  const T chi[nb_moments],
  const T planck_top,
  const T planck_bot)
{
  Layer<T> l;

  if (tau <= T(tau_min))
  {
    l.R = {T(0), T(0), T(0), T(0)};
    l.T_ = identity<T>();
    l.Eu = {T(0), T(0)};
    l.Ed = {T(0), T(0)};

    return l;
  }

  const T tiny = Traits<T>::tiny();

  //phase matrices Ws_ij = w_j/2 sum_l chi_l (2l+1) P_l(mu_i) P_l(mu_j) for scattering
  //within a hemisphere and Wd with P_l(-mu_j) = (-1)^l P_l(mu_j) across the hemispheres
  T p0[nb_moments], p1[nb_moments];
  legendre(T(mu0), p0);
  legendre(T(mu1), p1);

  Mat2<T> Ws = {T(0), T(0), T(0), T(0)};
  Mat2<T> Wd = {T(0), T(0), T(0), T(0)};

  for (int m = 0; m < nb_moments; ++m)
  {
    const T c = chi[m] * T(2*m + 1) * T(0.5) * T(weight0); //both weights are equal
    const T s = (m % 2 == 0) ? T(1) : T(-1);

    const T c00 = c * p0[m] * p0[m];
    const T c01 = c * p0[m] * p1[m];
    const T c11 = c * p1[m] * p1[m];

    Ws.a00 += c00; Ws.a01 += c01; Ws.a10 += c01; Ws.a11 += c11;
    Wd.a00 += s*c00; Wd.a01 += s*c01; Wd.a10 += s*c01; Wd.a11 += s*c11;
  }

  const T r0 = T(1) / T(mu0);
  const T r1 = T(1) / T(mu1);

  //dI+/dtau = alpha I+ - beta I-, dI-/dtau = beta I+ - alpha I- (homogeneous part)
  const Mat2<T> al = {(T(1) - omega*Ws.a00)*r0, -omega*Ws.a01*r0, -omega*Ws.a10*r1, (T(1) - omega*Ws.a11)*r1};
  const Mat2<T> bt = {omega*Wd.a00*r0, omega*Wd.a01*r0, omega*Wd.a10*r1, omega*Wd.a11*r1};

  const Mat2<T> amb = sub(al, bt);
  const Mat2<T> apb = add(al, bt);
  const Mat2<T> Q = mul(apb, amb);
  const Mat2<T> apb_inv = inverse(apb);

  //eigenvalues k^2 of Q: the larger root from the trace, the smaller via det/lambda1.
  //det(alpha - beta) vanishes for conservative scattering; with N = Ws + Wd, whose
  //rows sum to one, det(1 - omega N) = (1 - omega)(1 - omega (N_00 - N_10)).
  const T tr = Q.a00 + Q.a11;
  const T det_apb = apb.a00*apb.a11 - apb.a01*apb.a10;
  const T det_amb = one_minus_omega * (T(1) - omega*(Ws.a00 + Wd.a00 - Ws.a10 - Wd.a10)) * r0 * r1;
  const T det = det_apb * det_amb;
  const T disc = tr*tr - T(4)*det;
  const T lam1 = T(0.5) * (tr + Traits<T>::sqrt(disc > tiny ? disc : tiny));
  const T lam2 = det/lam1 > tiny ? det/lam1 : tiny;
  const T k1 = Traits<T>::sqrt(lam1);
  const T k2 = Traits<T>::sqrt(lam2);

  Vec2<T> X1, Y1, X2, Y2;
  eigenMode(Q, apb_inv, lam1, X1, Y1);
  eigenMode(Q, apb_inv, lam2, X2, Y2);

  //columns are the modes
  const Mat2<T> X = {X1.v0, X2.v0, X1.v1, X2.v1};
  const Mat2<T> Y = {Y1.v0, Y2.v0, Y1.v1, Y2.v1};
  const Vec2<T> k = {k1, k2};
  const Vec2<T> rk = {T(1)/k1, T(1)/k2};

  const Mat2<T> Delta = scaleColumns(Y, rk);                     //Gm - Gp = Y/k
  const Mat2<T> Gp = scale(sub(X, Delta), T(0.5));

  const Vec2<T> D = {-Traits<T>::expm1(-k1 * tau), -Traits<T>::expm1(-k2 * tau)};  //1 - exp(-k tau)

  const Mat2<T> GpD = scaleColumns(Gp, D);
  const Mat2<T> XD = scaleColumns(X, D);

  const Mat2<T> A1 = mul(scaleColumns(Delta, D), inverse(sub(X, GpD)));

  const bool thin = k1 * tau < T(tau_thin);

  Mat2<T> A2, A2p;

  if (thin)
  {
    //A2 = tau (alpha + beta) + A2' with
    //A2' = -X D2 K Y^-1 - X D K Y^-1 V (1 + V)^-1,  V = Gp D K Y^-1,  D2 = K tau - D
    const Vec2<T> D2 = {expm1Remainder(k1 * tau, D.v0), expm1Remainder(k2 * tau, D.v1)};
    const Mat2<T> KYi = scaleRows(k, inverse(Y));
    const Mat2<T> V = mul(GpD, KYi);
    const Mat2<T> VV = mul(V, inverse(add(identity<T>(), V)));

    A2p = sub(scale(mul(scaleColumns(X, D2), KYi), T(-1)), mul(mul(XD, KYi), VV));
    A2 = add(scale(apb, tau), A2p);
  }
  else
  {
    A2 = mul(XD, inverse(add(Delta, GpD)));
    A2p = A2; //unused
  }

  l.R = scale(sub(A2, A1), T(0.5));
  l.T_ = sub(identity<T>(), scale(add(A1, A2), T(0.5)));

  //thermal particular solution for B linear in tau: I+/- = B(tau) +/- d B1
  //with (1 - omega (Ws - Wd)) d = mu, i.e. d = (alpha + beta)^-1 (1, 1)
  const T b1 = (planck_bot - planck_top) / tau;
  const Vec2<T> d = mul(apb_inv, Vec2<T>{T(1), T(1)});

  //symmetric part: emission of the mean Planck function, 1 - (R + T) = A1
  const T sum_b = planck_top + planck_bot;
  const Vec2<T> a1_sum = rowSums(A1);
  const Vec2<T> sum_o = {sum_b * a1_sum.v0, sum_b * a1_sum.v1};

  //antisymmetric part: (B_top - B_bot)(1 - (R - T)) + 2 B1 (1 + R - T) d
  //= B1 [2 (A2 d - tau 1) + tau A2 1], with A2 d - tau 1 = A2' d for thin layers
  Vec2<T> dif_o;

  if (thin)
  {
    const Vec2<T> a2p_d = mul(A2p, d);
    const Vec2<T> a2_sum = rowSums(A2);
    dif_o = {b1 * (T(2) * a2p_d.v0 + tau * a2_sum.v0), b1 * (T(2) * a2p_d.v1 + tau * a2_sum.v1)};
  }
  else
  {
    const T dif_b = planck_top - planck_bot;
    const Vec2<T> a2_d = mul(A2, d);
    const Vec2<T> a2_sum = rowSums(A2);
    dif_o = {dif_b * (T(2) - a2_sum.v0) + T(2) * b1 * a2_d.v0, dif_b * (T(2) - a2_sum.v1) + T(2) * b1 * a2_d.v1};
  }

  l.Eu = {T(0.5)*(sum_o.v0 + dif_o.v0), T(0.5)*(sum_o.v1 + dif_o.v1)};
  l.Ed = {T(0.5)*(sum_o.v0 - dif_o.v0), T(0.5)*(sum_o.v1 - dif_o.v1)};

  return l;
}



//adds a layer on top of the slab below it, described by its reflection matrix Rb
//and the upward intensity Su it emits at its top; both are updated in place
template <typename T>
BEAR_HOST_DEVICE inline void addLayer(Mat2<T>& Rb, Vec2<T>& Su, const Layer<T>& l)
{
  const Mat2<T> RR = mul(Rb, l.R);
  const Mat2<T> one_minus_RR = {T(1) - RR.a00, -RR.a01, -RR.a10, T(1) - RR.a11};
  const Mat2<T> A = mul(l.T_, inverse(one_minus_RR));

  const Vec2<T> src = add(l.Eu, mul(A, mul(Rb, l.Ed)));
  const Vec2<T> Su_new = add(src, mul(A, Su));

  Rb = add(l.R, mul(A, mul(Rb, l.T_)));
  Su = Su_new;
}



//upward flux at the top from the upward intensities per stream, F = 2 pi sum_i w_i mu_i I_i
template <typename T>
BEAR_HOST_DEVICE inline T upwardFlux(const Vec2<T>& I)
{
  return T(2) * T(3.14159265358979323846) * (T(weight0)*T(mu0)*I.v0 + T(weight1)*T(mu1)*I.v1);
}


}
}

#endif
