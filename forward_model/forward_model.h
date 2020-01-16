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


#ifndef _forward_model_h
#define _forward_model_h

#include <vector>



namespace helios {

//forward declaration
class Retrieval;


//abstract class for the forward model
//a derived class *has to* implement the calcModel and calcModelGPU methods that will be called by the retrieval 
class ForwardModel{
  public:
    virtual ~ForwardModel() {}
    virtual bool calcModel(const std::vector<double>& parameter, std::vector<double>& spectrum) = 0;
    virtual bool calcModelGPU(const std::vector<double>& parameter, double* model_spectrum) = 0;
    
  protected:
    Retrieval* retrieval;

};


}


#endif

