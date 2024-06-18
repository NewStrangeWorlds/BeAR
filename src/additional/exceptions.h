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
 
