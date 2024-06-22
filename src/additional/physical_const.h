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


#ifndef _PHYSICAL_CONSTANTS_H
#define	_PHYSICAL_CONSTANTS_H


namespace bear{  namespace constants{

//physical constants in cgs units from the GSL library
double const planck_h = 6.62606896e-27;                     //Planck constant in g cm^2 / s
double const light_c = 2.99792458e10;                       //speed of light in cm / s
double const boltzmann_k = 1.3806504e-16;                   //Boltzmann constant in g cm^2 / K s^2
double const stefan_boltzmann = 5.67040047374e-5;           //Stefan-Boltzmann constant in g / K^4 s^3
double const mass_proton = 1.67262158e-24;                  //proton mass in g
double const mass_electron = 9.10938188e-28;                //electron mass in g
double const gravitation_const = 6.673e-8;                  //Gravitational constant in cm^3 / g s^2
double const amu = 1.6605e-24;                              //atomic mass unit
double const atomic_weight_c = 12.01;                       //atomic weight of carbon
double const gas_constant = 8.314472e7;                     //universal gas constant erg K^-1 mol^-1

//additional constants
double const radius_earth = 6.37123e8;                      //mean Earth radius in cm
double const radius_sun = 696342.*100000.;                  //Solar radius in cm
double const radius_jupiter = 7.1492e9;                     //mean Jupiter radius in cm
double const mass_earth = 5.97219e24 * 1000.;               //Earth mass in g
double const mass_jupiter = 1.89813e27 * 1000.;             //Jupiter mass in g

double const parsec = 3.08567758135e18;                     //1 parsec in cm

double const pi = 3.14159265358979323846;
double const deg_to_rad = pi / 180.;
double const rad_to_deg = 180./pi;


}}

#endif
