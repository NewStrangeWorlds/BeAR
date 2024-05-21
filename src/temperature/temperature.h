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


#ifndef _temperature_h
#define _temperature_h

#include <vector>
#include <cstddef>


namespace helios {


class Temperature{
  public:
    virtual ~Temperature() {}
    virtual bool calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      const std::vector<double>& pressure,
      std::vector<double>& temperature_profile) = 0;
    size_t nbParameters() {return nb_parameters;}
  protected:
    size_t nb_parameters {};
};


}
#endif 
