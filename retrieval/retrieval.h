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


#ifndef _retrieval_h
#define _retrieval_h

#include <vector>
#include <string>


#include "../forward_model/forward_model.h"
#include "../spectral_grid/spectral_grid.h"
#include "../observations/observations.h"
#include "../config/global_config.h"
#include "multinest_parameter.h"


namespace helios {

//forward declaration
class ForwardModel;
class BasicPrior;



//the main class that does the retrieval
//it contains the spectral grid, observations, priors, the forward model, and runs MulitNest
class Retrieval{
  public:
    Retrieval(GlobalConfig* global_config);
    ~Retrieval();

    GlobalConfig* config;
    SpectralGrid spectral_grid;
    std::vector<Observation> observations;           //object that holds the single observations

    virtual void doRetrieval();

    size_t nbObservations() {return nb_observations;}

  protected:
    std::vector<double> observation_data;            //combined vector of all observational data
    std::vector<double> observation_error;           //combined array of the corresponding observational errors

    double* observation_data_gpu = nullptr;          //pointer to the corresponding data on the GPU
    double* observation_error_gpu = nullptr;

    size_t nb_observations = 0;                      //number of observations
    size_t nb_observation_points = 0;                //total number of observational data points


    int* band_indices_gpu = nullptr;                 //spectral indices of the high-res points of all bands on the GPU
    int* band_sizes_gpu = nullptr;                   //number of spectral points for each band on the GPU
    int* band_start_index_gpu = nullptr;             //start index for each band on the GPU
    size_t nb_total_bands = 0;                       //total number of observational bands 


    double* model_spectrum_gpu = nullptr;            //pointer to the high-res spectrum on the GPU

    std::vector<double*> convolved_spectra;          //pointer to the convolved spectra on the GPU
    std::vector<int*> observation_spectral_indices;  //pointer to the spectral indices from the high-res spectrum for each observational bin on the GPU
    std::vector<double*> observation_wavelengths;    //pointer to the corresponding wavelengths on the GPU
    std::vector<double*> observation_profile_sigma;  //pointer to the instrument profiles on the GPU
    std::vector<int*> convolution_start_index;       //pointer to the start indices for the convolution on the GPU 
    std::vector<int*> convolution_end_index;         //pointer to the end indices for the convolution on the GPU


    std::vector<double*> band_spectrum_id;           //contains a pointer for each observation to either the high-res spectrum or one of its convolutions on the GPU

    bool loadObservations(const std::string file_folder, const std::vector<std::string>& file_list);
    bool loadObservationFileList(const std::string file_folder, std::vector<std::string>& file_list);

  private:
    ForwardModel* forward_model;
    std::vector<BasicPrior*> priors;
    
    void setPriors();
    void runMultinest();

    static void multinestLogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context);
    static void multinestLogLikeGPU(double *cube, int &nb_dim, int &nb_param, double &new_log_like, void *context);

    static void multinestDumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior,
                                double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context);
};


}


#endif
