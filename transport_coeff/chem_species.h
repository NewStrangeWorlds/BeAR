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


#ifndef _chem_species_h
#define _chem_species_h

#include <vector>
#include <string>

unsigned const int _TOTAL  = 0;
unsigned const int _H      = 1;
unsigned const int _H2     = 2;
unsigned const int _He     = 3;
unsigned const int _H2O    = 4;
unsigned const int _CO2    = 5;
unsigned const int _CO     = 6;
unsigned const int _CH4    = 7;
unsigned const int _HCN    = 8;
unsigned const int _NH3    = 9;
unsigned const int _C2H2   = 10;
unsigned const int _N2     = 11;
unsigned const int _Na     = 12;
unsigned const int _K      = 13;
unsigned const int _O3      = 14;
unsigned const int _N2O      = 15;
unsigned const int _H2S      = 16;


const std::vector<std::string> chemical_symbols{"Total",
                                                "H",
                                                "H2",
                                                "He",
                                                "H2O1",
                                                "C1O2",
                                                "C1O1",
                                                "C1H4",
                                                "C1H1N1",
                                                "H3N1",
                                                "C2H2",
                                                "N2",
                                                "Na",
                                                "K",
                                                "O3",
                                                "N2O1",
                                                "H2S1"};


#endif
