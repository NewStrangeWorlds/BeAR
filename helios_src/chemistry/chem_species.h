/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
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


namespace helios{

enum chemical_species_id {_TOTAL, _H, _He, _C, _O, _H2, _H2O, _CO2, _CO, _CH4, _HCN, _NH3, _C2H2, _N2, _Na, _K, _H2S};


struct chemistry_data{
  chemical_species_id id;
  std::string symbol;
  std::string fastchem_symbol;
  double molecular_weight;
};


//the chemical data is put into the constants namespace 
namespace constants{

const std::vector<chemistry_data> species_data{ {_TOTAL, "Total", "Total",  0.0},
                                                {_H,     "H",     "H",      1.00784},
                                                {_He,    "He",    "He",     4.002602},
                                                {_C,     "C",     "C",      12.0107},
                                                {_O,     "O",     "O",      15.999},
                                                {_H2,    "H2",    "H2",     2.01588},
                                                {_H2O,   "H2O",   "H2O1",   18.01528},
                                                {_CO2,   "CO2",   "C1O2",   44.01},
                                                {_CO,    "CO",    "C1O1",   28.0101},
                                                {_CH4,   "CH4",   "C1H4",   16.04246},
                                                {_HCN,   "HCN",   "C1H1N1", 27.0253},
                                                {_NH3,   "NH3",   "H3N1",   17.03052},
                                                {_C2H2,  "C2H2",  "C2H2",   26.04},
                                                {_N2,    "N2",    "N2",     28.0134},
                                                {_Na,    "Na",    "Na",     22.98977},
                                                {_K,     "K",     "K",      39.0983},
                                                {_H2S,   "H2S",   "H2S1",   34.09099}
                                              };
}
}

#endif
