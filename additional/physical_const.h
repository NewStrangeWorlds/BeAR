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


#ifndef _PHYSICAL_CONSTANTS_H
#define	_PHYSICAL_CONSTANTS_H

#include <gsl/gsl_const_cgsm.h>

namespace helios{

//physical constants in cgs units from the GSL library
double const CONST_H = GSL_CONST_CGSM_PLANCKS_CONSTANT_H;                 //Planck constant
double const CONST_C = GSL_CONST_CGSM_SPEED_OF_LIGHT;                     //speed of light
double const CONST_K = GSL_CONST_CGSM_BOLTZMANN;                          //Boltzmann constant
double const CONST_SIGMA = GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT;      //Stefan-Boltzmann constant
double const CONST_MP = GSL_CONST_CGSM_MASS_PROTON;                       //proton mass
double const CONST_G = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT;             //Gravitational constant
double const CONST_AMU = 1.6605e-24;                                      //atomic mass unit
double const CONST_AW_C = 12.01;                                          //atomic weight of carbon
double const CONST_R = GSL_CONST_CGSM_MOLAR_GAS;                          //universal gas constant erg K^-1 mol^-1
double const CONST_E_MASS = GSL_CONST_CGSM_MASS_ELECTRON;                 //electron mass in g


//additional constants
double const CONST_EARTH_RADIUS = 6.37123e8;                              //mean Earth radius in cm
double const CONST_SOLAR_RADIUS = 696342.*100000.;                        //Solar radius in cm
double const CONST_JUPITER_RADIUS = 7.1492e9;                             //mean Jupiter radius in cm

double const CONST_PARSEC = 3.08567758135e18;                             //1 parsec in cm

double const CONST_EARTH_MASS = 5.97219e24 * 1000.;                       //Earth mass in g
double const CONST_JUPITER_MASS = 1.89813e27 * 1000.;                     //Jupiter mass in g

double const CONST_PI = 3.14159265358979323846;
double const DEG_TO_RAD = CONST_PI / 180.;
double const RAD_TO_DEG = 180./CONST_PI;


}

#endif
