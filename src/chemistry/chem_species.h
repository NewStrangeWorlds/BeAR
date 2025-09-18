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


#ifndef _chem_species_h
#define _chem_species_h

#include <vector>
#include <string>


namespace bear{

enum chemical_species_id {_TOTAL, _H, _He, _C, _O, _Fe, _Fep, _Ca, _Ti, _Tip, _H2, _H2O, _CO2, _CO, _CH4, _HCN, _NH3, _C2H2, _N2, _Na, _K, _H2S, _Hm, _TiO, _VO, _FeH, _SH, _MgO, _AlO, _CaO, _CrH, _MgH, _CaH, _TiH, _OH, _e, _V, _Vp, _Mn, _Si, _Cr, _Crp, _SiO, _SiO2, _SO2, _CS2};


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
                                                {_Fe,    "Fe",    "Fe",     55.845},
                                                {_Fep,   "Fe+",   "Fe+",    55.845},
                                                {_Ca,    "Ca",    "Ca",     40.078},
                                                {_Ti,    "Ti",    "Ti",     47.867},
                                                {_Tip,   "Ti+",   "Ti+",    47.867},
                                                {_H2,    "H2",    "H2",     2.01588},
                                                {_H2O,   "H2O",   "H2O1",   18.01528},
                                                {_CO2,   "CO2",   "C1O2",   44.01},
                                                {_CO,    "CO",    "C1O1",   28.0101},
                                                {_CH4,   "CH4",   "C1H4",   16.04246},
                                                {_HCN,   "HCN",   "C1H1N1_1", 27.0253},
                                                {_NH3,   "NH3",   "H3N1",   17.03052},
                                                {_C2H2,  "C2H2",  "C2H2",   26.04},
                                                {_N2,    "N2",    "N2",     28.0134},
                                                {_Na,    "Na",    "Na",     22.98977},
                                                {_K,     "K",     "K",      39.0983},
                                                {_H2S,   "H2S",   "H2S1",   34.09099},
                                                {_Hm,    "H-",    "H1-",    1.00784},
                                                {_TiO,   "TiO",   "O1Ti1",  63.8664},
                                                {_VO,    "VO",    "O1V1",   66.9409},
                                                {_FeH,   "FeH",   "H1Fe1",  56.853},
                                                {_SH,    "SH",    "H1S1",   34.08},
                                                {_MgO,   "MgO",   "Mg1O1",  40.3044}, 
                                                {_AlO,   "AlO",   "Al1O1",  42.981}, 
                                                {_CaO,   "CaO",   "Ca1O1",  56.0774}, 
                                                {_CrH,   "CrH",   "Cr1H1",  54.0040},
                                                {_MgH,   "MgH",   "H1Mg1",  26.3209},
                                                {_CaH,   "CaH",   "Ca1H1",  41.0859},
                                                {_TiH,   "TiH",   "H1Ti1",  48.87484},
                                                {_OH,    "OH",    "H1O1",   17.008},
                                                {_e,     "e-",    "e-",     5.4857990907e-4},
                                                {_V,     "V",     "V",      50.9415},
                                                {_Vp,    "V+",    "V1+",    50.9415},
                                                {_Mn,    "Mn",    "Mn",     54.938044},
                                                {_Si,    "Si",    "Si",     28.085},
                                                {_Cr,    "Cr",    "Cr",     51.996},
                                                {_Crp,   "Cr+",   "Cr1+",   51.996},
                                                {_SiO,   "SiO",   "O1Si1",  44.08},
                                                {_SiO2,  "SiO2",  "O2Si1",  60.08},
                                                {_SO2,   "SO2",   "O2S1",   64.066},
                                                {_CS2,   "CS2",   "C1S2",   76.139}
                                              };
}
}

#endif
