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


#ifndef _forward_model_h
#define _forward_model_h

#include <vector>

#include "../retrieval/priors.h"


namespace bear {


//abstract class for the forward model
//a derived class *has to* implement all the various, virtual methods
class ForwardModel{
  public:
    virtual ~ForwardModel() {}
    //calculate a model on the CPU
    //the return value signals the retrieval to neglect this model
    virtual bool calcModel(
      const std::vector<double>& parameter, 
      std::vector<double>& spectrum, 
      std::vector<double>& spectrum_bands) = 0;
    //calculate a model on the GPU
    //the return value signals the retrieval to neglect this model
    virtual bool calcModelGPU(
      const std::vector<double>& parameter, 
      double* model_spectrum, 
      double* model_spectrum_bands) = 0;
    //model-specific post process
    virtual void postProcess(
      const std::vector< std::vector<double> >& model_parameter, 
      const std::vector< std::vector<double> >& model_spectrum_bands,
      const size_t best_fit_model) = 0;
    //converts a high-resolution spectrum to a model observation (e.g. a transit depth)
    virtual std::vector<double> convertSpectrumToModel(const std::vector<double>& spectrum) = 0;
    //model-specific tests
    virtual bool testModel(
      const std::vector<double>& parameter, 
      double* model_spectrum) = 0;
  protected:
    virtual void setPriors(Priors* priors) = 0;
};


}


#endif

