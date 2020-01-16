/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2019 Daniel Kitzmann
*
* Helios-r2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Helios-r2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* Helios-r2 directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#include "aux_functions.h"
#include "physical_const.h"

#include <cmath>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


namespace helios{ namespace aux{


double voigt(double xx, double sigma, double lg, int r)
{
   // Computation of Voigt function (normalised).
   // Voigt is a convolution of 
   // gauss(xx) = 1/(sqrt(2*pi)*sigma) * exp(xx*xx/(2*sigma*sigma)
   // and
   // lorentz(xx) = (1/pi) * (lg/2) / (xx*xx + lg*lg/4)
   // functions.
   //
   // The Voigt function is known to be the real part of Faddeeva function also
   // called complex error function [2].
   //
   // The algoritm was developed by J. Humlicek [1].
   // This code is based on fortran code presented by R. J. Wells [2].
   // Translated and adapted by Miha D. Puc
   //
   // To calculate the Faddeeva function with relative error less than 10^(-r).
   // r can be set by the the user subject to the constraints 2 <= r <= 5.
   //
   // [1] J. Humlicek, JQSRT, 21, 437 (1982).
   // [2] R.J. Wells "Rapid Approximation to the Voigt/Faddeeva Function and its
   // Derivatives" JQSRT 62 (1999), pp 29-48.
   // http://www-atm.physics.ox.ac.uk/user/wells/voigt.html

   if ((sigma < 0 || lg < 0) || (sigma==0 && lg==0)) {
      return 0;  // Not meant to be for those who want to be thinner than 0
   }

   if (sigma == 0) {
      return lg * 0.159154943  / (xx*xx + lg*lg /4); //pure Lorentz
   }

   if (lg == 0) {   //pure gauss
      return 0.39894228 / sigma * std::exp(-xx*xx / (2*sigma*sigma));
   }

   double x, y, k;
   x = xx / sigma / 1.41421356;
   y = lg / 2. / sigma / 1.41421356;

   //x = xx / sigma / sqrt(log(2.0));
   //y = lg / sigma / sqrt(log(2.0));

   double r0, r1;

   if (r < 2) r = 2;
   if (r > 5) r = 5;

   r0=1.51 * exp(1.144 * (double)r);
   r1=1.60 * exp(0.554 * (double)r);

   // Constants

   const double rrtpi = 0.56418958;  // 1/SQRT(pi)

   double y0, y0py0, y0q;                      // for CPF12 algorithm
   y0 = 1.5;
   y0py0 = y0 + y0;
   y0q = y0 * y0;

   double c[6] = { 1.0117281, -0.75197147, 0.012557727, 0.010022008, -0.00024206814, 0.00000050084806};
   double s[6] = { 1.393237, 0.23115241, -0.15535147, 0.0062183662, 0.000091908299, -0.00000062752596};
   double t[6] = { 0.31424038, 0.94778839, 1.5976826, 2.2795071, 3.0206370, 3.8897249};

   // Local variables

   int j;                                        // Loop variables
   int rg1, rg2, rg3;                            // y polynomial flags
   double abx, xq, yq, yrrtpi;                 // --x--, x^2, y^2, y/SQRT(pi)
   double xlim0, xlim1, xlim2, xlim3, xlim4;   // --x-- on region boundaries
   double a0=0, d0=0, d2=0, e0=0, e2=0, e4=0, h0=0, h2=0, h4=0, h6=0;// W4 temporary variables
   double p0=0, p2=0, p4=0, p6=0, p8=0, z0=0, z2=0, z4=0, z6=0, z8=0;
   double xp[6], xm[6], yp[6], ym[6];          // CPF12 temporary values
   double mq[6], pq[6], mf[6], pf[6];
   double d, yf, ypy0, ypy0q;

   //***** Start of executable code *****************************************

   rg1 = 1;  // Set flags
   rg2 = 1;
   rg3 = 1;
   yq = y * y;  // y^2
   yrrtpi = y * rrtpi;  // y/SQRT(pi)

   // Region boundaries when both k and L are required or when R<>4

   xlim0 = r0 - y;
   xlim1 = r1 - y;
   xlim3 = 3.097 * y - 0.45;
   xlim2 = 6.8 - y;
   xlim4 = 18.1 * y + 1.65;
   if ( y <= 1e-6 ) {                      // When y<10^-6 avoid W4 algorithm
      xlim1 = xlim0;
      xlim2 = xlim0;
   }

   abx = fabs(x);                                // |x|
   xq = abx * abx;                               // x^2
   if ( abx > xlim0 ) {                          // Region 0 algorithm
      k = yrrtpi / (xq + yq);
   } else if ( abx > xlim1 ) {                   // Humlicek W4 Region 1
      if ( rg1 != 0 ) {                          // First point in Region 1
         rg1 = 0;
         a0 = yq + 0.5;                          // Region 1 y-dependents
         d0 = a0*a0;
         d2 = yq + yq - 1.0;
      }
      d = rrtpi / (d0 + xq*(d2 + xq));
      k = d * y * (a0 + xq);
   } else if ( abx > xlim2 ) {                   // Humlicek W4 Region 2
      if ( rg2 != 0 ) {                          // First point in Region 2
         rg2 = 0;
         h0 = 0.5625 + yq * (4.5 + yq * (10.5 + yq * (6.0 + yq)));
                                                 // Region 2 y-dependents
         h2 = -4.5 + yq * (9.0 + yq * ( 6.0 + yq * 4.0));
         h4 = 10.5 - yq * (6.0 - yq * 6.0);
         h6 = -6.0 + yq * 4.0;
         e0 = 1.875 + yq * (8.25 + yq * (5.5 + yq));
         e2 = 5.25 + yq * (1.0 + yq * 3.0);
         e4 = 0.75 * h6;
      }
      d = rrtpi / (h0 + xq * (h2 + xq * (h4 + xq * (h6 + xq))));
      k = d * y * (e0 + xq * (e2 + xq * (e4 + xq)));
   } else if ( abx < xlim3 ) {                   // Humlicek W4 Region 3
      if ( rg3 != 0 ) {                          // First point in Region 3
         rg3 = 0;
         z0 = 272.1014 + y * (1280.829 + y *
                              (2802.870 + y *
                               (3764.966 + y *
                                (3447.629 + y *
                                 (2256.981 + y *
                                  (1074.409 + y *
                                   (369.1989  + y *
                                    (88.26741 + y *
                                     (13.39880 + y)
                                     ))))))));   // Region 3 y-dependents
         z2 = 211.678 + y * (902.3066 + y *
                             (1758.336 + y *
                              (2037.310 + y *
                               (1549.675 + y *
                                (793.4273 + y *
                                 (266.2987 + y *
                                  (53.59518 + y * 5.0)
                                  ))))));
         z4 = 78.86585 + y * (308.1852 + y *
                              (497.3014 + y *
                               (479.2576 + y *
                                (269.2916 + y *
                                 (80.39278 + y * 10.0)
                                 ))));
         z6 = 22.03523 + y * (55.02933 + y *
                              (92.75679 + y *
                               (53.59518 + y * 10.0)
                               ));
         z8 = 1.496460 + y * (13.39880 + y * 5.0);
         p0 = 153.5168 + y * (549.3954 + y *
                              (919.4955 + y *
                               (946.8970 + y *
                                (662.8097 + y *
                                 (328.2151 + y *
                                  (115.3772 + y *
                                   (27.93941 + y *
                                    (4.264678 + y * 0.3183291)
                                    )))))));
         p2 = -34.16955 + y * (-1.322256+ y *
                               (124.5975 + y *
                                (189.7730 + y *
                                 (139.4665 + y *
                                  (56.81652 + y *
                                   (12.79458 + y * 1.2733163)
                                   )))));
         p4 = 2.584042 + y * (10.46332 + y *
                              (24.01655 + y *
                               (29.81482 + y *
                                (12.79568 + y * 1.9099744)
                                )));
         p6 = -0.07272979 + y * (0.9377051 + y *
                                 (4.266322 + y * 1.273316));
         p8 = 0.0005480304 + y * 0.3183291;
      }
      d = 1.7724538 / (z0 + xq * (z2 + xq * (z4 + xq * (z6 + xq * (z8 + xq)))));
      k = d * (p0 + xq * (p2 + xq * (p4 + xq * (p6 + xq * p8))));
   } else {                             // Humlicek CPF12 algorithm
      ypy0 = y + y0;
      ypy0q = ypy0 * ypy0;
      k = 0.0;
      for (j = 0; j <= 5; j++) {
         d = x - t[j];
         mq[j] = d * d;
         mf[j] = 1.0 / (mq[j] + ypy0q);
         xm[j] = mf[j] * d;
         ym[j] = mf[j] * ypy0;
         d = x + t[j];
         pq[j] = d * d;
         pf[j] = 1.0 / (pq[j] + ypy0q);
         xp[j] = pf[j] * d;
         yp[j] = pf[j] * ypy0;
      }
      if ( abx <= xlim4 ) {                      // Humlicek CPF12 Region I
         for (j = 0; j <= 5; j++) {
            k = k + c[j]*(ym[j]+yp[j]) - s[j]*(xm[j]-xp[j]) ;
         }
      } else {                                   // Humlicek CPF12 Region II
         yf = y + y0py0;
         for ( j = 0; j <= 5; j++) {
            k = k + (c[j] *
                 (mq[j] * mf[j] - y0 * ym[j])
                    + s[j] * yf * xm[j]) / (mq[j]+y0q)
                 + (c[j] * (pq[j] * pf[j] - y0 * yp[j])
                   - s[j] * yf * xp[j]) / (pq[j]+y0q);
         }
         k = y * k + exp( -xq );
      }
   }

   return k / 2.506628 / sigma; // Normalize by dividing by sqrt(2*pi)*sigma
}





double planckFunctionWavenumber(const double temperature, const double wavenumber)
{

  return 2. * CONST_H * CONST_C * CONST_C * wavenumber*wavenumber*wavenumber / ( exp(CONST_H * wavenumber * CONST_C / CONST_K / temperature) - 1.0);


}



double linearInterpolation(const double x1, const double x2, const double y1, const double y2, const double x)
{

  return y1 + (y2 - y1) * (x - x1)/(x2 - x1);

}



double voigtProfile(const double x, const double gaussian_width, const double lorentz_width)
{


  double y = voigt(x, gaussian_width, lorentz_width, 5); 

  
  return y;
}



std::vector<double> interpolateToWavenumberGrid(const std::vector<double>& wavenumber_data, const std::vector<double>& data,
                                                const std::vector<double>& wavenumber_grid,
                                                const double outside_range_value,
                                                const bool interpolate_log)
{
  unsigned int nb_data_points = wavenumber_data.size();


  gsl_interp *spline = gsl_interp_alloc (gsl_interp_linear, nb_data_points);

  double *xa = new double[nb_data_points];
  double *ya = new double[nb_data_points];


  for (unsigned int i=0; i<nb_data_points; i++)
  {
    xa[i] = wavenumber_data[i];

    if (interpolate_log)
      ya[i] = log10(data[i]);
    else
      ya[i] = data[i];
  }



  gsl_interp_init (spline, xa, ya, nb_data_points);

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();



  std::vector<double> result(wavenumber_grid.size(), 0.0);


  for (unsigned int i=0; i<wavenumber_grid.size(); ++i)
  {
    if (wavenumber_grid[i] > wavenumber_data.back())
    {
      result[i] = outside_range_value;

      continue;
    }

    if (wavenumber_grid[i] < wavenumber_data.front())
    {
      result[i] = outside_range_value;

      continue;
    }

    if (interpolate_log)
      result[i] = pow(10., gsl_interp_eval (spline, xa, ya, wavenumber_grid[i],acc));
    else
      result[i] = gsl_interp_eval (spline, xa, ya, wavenumber_grid[i],acc);
  }


  gsl_interp_accel_free (acc);
  gsl_interp_free (spline);


  delete[] xa;
  delete[] ya;


  return result;
}




double normalDistribution(const double mu, const double sigma, const double x)
{

   const double c = std::sqrt(1.0/(2.0 * CONST_PI));

   return c / sigma * std::exp(- (x - mu)*(x - mu) / (2.0 * sigma*sigma));

}


double normalDistribution(const double sigma, const double x)
{

   const double c = std::sqrt(1.0/(2.0 * CONST_PI));

   return c / sigma * std::exp(- x*x / (2.0 * sigma*sigma));

}




}}
