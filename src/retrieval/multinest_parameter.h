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


#ifndef _multinest_parameter_h
#define _multinest_parameter_h

#include <string>
#include "../config/global_config.h"


namespace bear {


//A struct with all available Multinest parameters
struct MultinestParameter{
  MultinestParameter(GlobalConfig* config);
  void loadParameterFile(const std::string& file_name);  //not implemented yet

  int is = 1;             //do Nested Importance Sampling?
  int mmodal = 0;         //do mode separation?
  int ceff = 0;           //run in constant efficiency mode?
  int nlive = 400;        //number of live points
  double efr = 0.8;       //set the required efficiency
  double tol = 0.5;       //tol, defines the stopping criteria
  int ndims = 6;          //dimensionality (no. of free parameters)
  int nPar = 6;           //total no. of parameters including free & derived parameters
  int nClsPar = 6;        //no. of parameters to do mode separation on
  int updInt = 1000;      //after how many iterations feedback is required & the output files should be updated
  double Ztol = -1E90;    //all the modes with logZ < Ztol are ignored
  int maxModes = 100;     //expected max no. of modes (used only for memory allocation)
  int pWrap[100];         //which parameters to have periodic boundary conditions?
  char root[1000];         //root for output files
  int seed = -1;          //random no. generator seed, if < 0 then take the seed from system clock
  int fb = 1;             //need feedback on standard output?
  int resume = 0;         //resume from a previous job?
  int outfile = 1;        //write output files?
  int initMPI = 0;        //initialize MPI routines?, relevant only if compiling with MPI
  double logZero = -1E90; //points with loglike < logZero will be ignored by MultiNest
  int maxiter = 0;        //max no. of iterations, a non-positive value means infinity.
  void *context = 0;      //not required by MultiNest, any additional information user wants to pass
};




}


#endif
