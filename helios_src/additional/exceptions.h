/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2022 Daniel Kitzmann
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


#ifndef _EXCEPTIONS_H
#define	_EXCEPTIONS_H


#include <iostream>
#include <exception>
#include <stdexcept>


namespace helios{

class FileNotFound : public std::runtime_error {
  public:
    FileNotFound(
      const std::string where, 
      const std::string what) : std::runtime_error("Critical Error - Aborting!") {
        std::cout << "File " << what << " in " << where << " not found!\n";
      }
};


class InvalidInput : public std::runtime_error {
  public:
    InvalidInput(
      const std::string where, 
      const std::string what) : std::runtime_error("Critical Error - Aborting!") {
        std::cout << "Invalid input in " << where << ":\n" << what << "\n";
      }
};


}

#endif
 
