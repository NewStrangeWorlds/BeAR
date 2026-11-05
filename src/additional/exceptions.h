// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _EXCEPTIONS_H
#define	_EXCEPTIONS_H


#include <iostream>
#include <exception>
#include <stdexcept>


namespace bear{

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
 
