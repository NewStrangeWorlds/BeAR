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


#ifndef _retrieval_h
#define _retrieval_h

#include <vector>
#include <string>
#include <iostream>


#include "../spectral_grid/spectral_grid.h"
#include "../observations/observations.h"
#include "../config/global_config.h"
#include "priors.h"


namespace bear {

extern bool stop_model;

void signalHandler(int sig);

//forward declaration
class ForwardModel;



//the main class that does the retrieval
class Retrieval{
  public:
    Retrieval(
      GlobalConfig* global_config, 
      const std::string additional_observation_file);
    Retrieval(GlobalConfig* global_config);
    ~Retrieval();
    
    GlobalConfig* config = nullptr;
    SpectralGrid spectral_grid;
    std::vector<Observation> observations;
    Priors priors;
    
    size_t nb_observations = 0;

    virtual bool run();
  protected:
    ForwardModel* forward_model = nullptr;
    
    ForwardModel* selectForwardModel(
       const std::string model_description);
    void setAdditionalPriors();
    void loadObservations(
      const std::string file_folder, 
      const std::vector<std::string>& file_list,
      const std::vector<std::string>& modifier_list);
    void loadObservationFileList(
      const std::string file_folder, 
      std::vector<std::string>& file_list,
      std::vector<std::string>& modifier_list);
  private:
    void convertHypercubeParameters(
      double *cube,
      const size_t nb_param,
      std::vector<double>& parameter,
      std::vector<double>& physical_parameter);
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
      std::vector<double*> model_spectrum,
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
