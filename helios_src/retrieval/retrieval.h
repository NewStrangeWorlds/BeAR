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


#ifndef _retrieval_h
#define _retrieval_h

#include <vector>
#include <string>
#include <iostream>


//#include "../forward_model/forward_model.h"
#include "../spectral_grid/spectral_grid.h"
#include "../observations/observations.h"
#include "../config/global_config.h"
#include "priors.h"


namespace helios {

extern bool stop_model;

void signalHandler(int sig);

//forward declaration
class ForwardModel;



//the main class that does the retrieval
//it contains the spectral grid, observations, priors, the forward model, and runs MulitNest
class Retrieval{
  public:
    Retrieval(GlobalConfig* global_config);
    ~Retrieval();
    Priors priors;
    
    GlobalConfig* config;
    SpectralGrid spectral_grid;
    std::vector<Observation> observations;

    virtual bool doRetrieval();

    size_t nbObservations() {return observations.size();}

    std::vector<double> observation_data;               //combined vector of all observational data
    std::vector<double> observation_error;              //combined array of the corresponding observational errors
    std::vector<double> observation_likelihood_weight; //combined vector of likelihood weights

    double* observation_data_gpu = nullptr;             //pointer to the corresponding data on the GPU
    double* observation_error_gpu = nullptr;
    double* observation_likelihood_weight_gpu = nullptr;

    size_t nb_observations = 0;
    size_t nb_observation_points = 0;

    double* model_spectrum_gpu = nullptr;            //pointer to the high-res spectrum on the GPU
  protected:
   
    ForwardModel* forward_model = nullptr;
    void setAdditionalPriors();
    void loadObservations(
      const std::string file_folder, 
      const std::vector<std::string>& file_list);
    void loadObservationFileList(
      const std::string file_folder, 
      std::vector<std::string>& file_list);
    
    ForwardModel* selectForwardModel(const std::string model_description);
  private:
    static void multinestLogLike(
      double *Cube, 
      int &ndim, 
      int &npars, 
      double &lnew, 
      void *context);
    static void multinestLogLikeGPU(
      double *cube, 
      int &nb_dim, 
      int &nb_param, 
      double &new_log_like, 
      void *context);
    double logLikeDev(
      double* model_spectrum,
      const double error_inflation_coefficient);

    static void multinestDumper(
      int &nSamples, 
      int &nlive, 
      int &nPar, 
      double **physLive, 
      double **posterior,
      double **paramConstr, 
      double &maxLogLike, 
      double &logZ, 
      double &INSlogZ, 
      double &logZerr, 
      void *context);
};


}


#endif
