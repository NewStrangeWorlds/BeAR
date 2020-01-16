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


#include "retrieval.h"


#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>


#include "multinest_parameter.h"
#include "prior.h"
#include "../forward_model/forward_model.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../CUDA_kernels/log_like_kernels.h"
#include "../CUDA_kernels/band_integration_kernels.h"
#include "../CUDA_kernels/convolution_kernels.h"
#include "../observations/observations.h"

#include "../additional/physical_const.h"



namespace helios{


// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
// context						void pointer, any additional information
//
// Input/Output arguments
// Cube[npars] 						= on entry has the ndim parameters in unit-hypercube
//	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//
// Output arguments
// lnew 						= loglikelihood
void Retrieval::multinestLogLike(double *cube, int &ndim, int &nb_param, double &lnew, void *context)
{

  //context contains a pointer to the retrieval object
  //we now recast it accordingly to use the retrieval object here
  Retrieval *retrieval_ptr = static_cast<Retrieval*>(context);


  //convert the normalised cube values into the real parameter values
  //write back the parameter values into the cube for MultiNest output
  std::vector<double> parameter(nb_param, 0.0);

  for (size_t i=0; i<parameter.size(); ++i)
  {
    parameter[i] = retrieval_ptr->priors[i]->priorParameterValue(cube[i]);
    cube[i] = parameter[i];
  }


  //run the forward model with the parameter set to obtain a high-res model spectrum
  std::vector<double> model_spectrum(retrieval_ptr->spectral_grid.nbSpectralPoints(), 0.0);
  bool neglect = retrieval_ptr->forward_model->calcModel(parameter, model_spectrum);



  //integrate the high-res spectrum to observational bands
  //and convolve if neccessary 
  std::vector<double> model_spectrum_bands(retrieval_ptr->nb_observation_points, 0.0);
  std::vector<double> model_spectrum_convolved(retrieval_ptr->nb_observation_points, 0.0);
  std::vector<double>::iterator it = model_spectrum_bands.begin();


  for (size_t i=0; i<retrieval_ptr->nb_observations; ++i)
  {
    std::vector<double> observation_bands;


    if (retrieval_ptr->observations[i].instrument_profile_fwhm.size() == 0)
      retrieval_ptr->observations[i].spectral_bands.bandIntegrateSpectrum(model_spectrum, observation_bands);
    else
    {
      retrieval_ptr->observations[i].spectral_bands.convolveSpectrum(model_spectrum, model_spectrum_convolved);
      retrieval_ptr->observations[i].spectral_bands.bandIntegrateSpectrum(model_spectrum_convolved, observation_bands);
    }
    

    std::copy(observation_bands.begin(), observation_bands.end(), it);
    it += observation_bands.size();
  }


  //print the current model parameter values to the terminal
  if (retrieval_ptr->config->multinest_print_iter_values)
  {
    std::cout << "model ";
    
    for (auto & i : parameter) std::cout << i << "   ";

    std::cout << "\n";
  }

  
  //calculate the new likelihood value with the error inflation
  const double error_exponent = parameter.back();
  const double error_inflation = std::pow(10, error_exponent);


  double log_like = 0;

  for (size_t i=0; i<retrieval_ptr->nb_observation_points; ++i)
  {
    //Eq. 22 from Paper I
    const double error_square = retrieval_ptr->observation_error[i]*retrieval_ptr->observation_error[i] + error_inflation;
    
    //Eq. 23 from Paper I
    log_like += - 0.5 * std::log(error_square* 2.0 * CONST_PI)
                - 0.5 * (retrieval_ptr->observation_data[i] - model_spectrum_bands[i])*(retrieval_ptr->observation_data[i] - model_spectrum_bands[i]) / error_square;
  }

  
  //if the forward model tells us to neglect the current set of parameters, set the likelihood to a low value
  if (neglect == true) log_like = -1e30;


  //print the calculated log-like value
  if (retrieval_ptr->config->multinest_print_iter_values)
    std::cout << log_like << "\n";
  

  lnew = log_like;
}





void Retrieval::multinestLogLikeGPU(double *cube, int &nb_dim, int &nb_param, double &lnew, void *context)
{

  //context contains a pointer to the retrieval object
  //we now recast it accordingly to use the retrieval object here
  Retrieval *retrieval_ptr = static_cast<Retrieval*>(context);


  //convert the normalised cube values into the real parameter values
  //write back the parameter values into the cube for multinest output
  std::vector<double> parameter(nb_param, 0.0);

  for (size_t i=0; i<parameter.size(); ++i)
  {
    parameter[i] = retrieval_ptr->priors[i]->priorParameterValue(cube[i]);
    cube[i] = parameter[i];
  }

  
  //print the parameter values
  if (retrieval_ptr->config->multinest_print_iter_values)
  {
    std::cout << "model ";

    for (auto & i : parameter) std::cout << i << "   ";

    std::cout << "\n";
  }
  
  
  //allocate the memory for the spectra on the GPU
  size_t nb_points = retrieval_ptr->spectral_grid.nbSpectralPoints();

  //pointer to the spectrum on the GPU
  double* model_spectrum_bands = nullptr;
  allocateOnDevice(model_spectrum_bands, retrieval_ptr->nb_total_bands);


  //intialise the high-res spectrum on the GPU (set it to 0)
  intializeOnDevice(retrieval_ptr->model_spectrum_gpu, nb_points);


  //call the forward model
  bool neglect = retrieval_ptr->forward_model->calcModelGPU(parameter, retrieval_ptr->model_spectrum_gpu);

 
  //convolve the spectra of neccesary 
  for (size_t i=0; i<retrieval_ptr->nb_observations; ++i)
  {
    size_t nb_points_observation = retrieval_ptr->observations[i].spectral_bands.wavenumbers.size();

    if (retrieval_ptr->observations[i].instrument_profile_fwhm.size() != 0) 
      convolveSpectrumGPU(retrieval_ptr->model_spectrum_gpu, 
                          retrieval_ptr->observation_wavelengths[i], 
                          retrieval_ptr->observation_profile_sigma[i], 
                          retrieval_ptr->observation_spectral_indices[i],
                          retrieval_ptr->convolution_start_index[i], 
                          retrieval_ptr->convolution_end_index[i], 
                          nb_points_observation, 
                          retrieval_ptr->convolved_spectra[i]);
  }

  
  //integrate the high-res spectrum to observational bands on the GPU
  bandIntegrationGPU(retrieval_ptr->band_spectrum_id,
                     retrieval_ptr->band_indices_gpu,
                     retrieval_ptr->band_sizes_gpu,
                     retrieval_ptr->band_start_index_gpu,
                     retrieval_ptr->nb_total_bands,
                     retrieval_ptr->spectral_grid.wavenumber_list_gpu,
                     retrieval_ptr->spectral_grid.wavelength_list_gpu,
                     model_spectrum_bands);


  //calculate the likelihood on the GPU with the error inflation
  const double error_exponent = parameter.back();
  const double error_inflation = std::pow(10, error_exponent);

  double new_log_like = logLikeHost(retrieval_ptr->observation_data_gpu,
                                    retrieval_ptr->observation_error_gpu,
                                    model_spectrum_bands,
                                    retrieval_ptr->nb_total_bands, 
                                    error_inflation);
  

  //if the forward model tells us to neglect the current set of parameters, set the likelihood to a low value
  if (neglect == true) new_log_like = -1e30;
  

  //delete the spectra from the GPU
  deleteFromDevice(model_spectrum_bands);

  
  //print the calculated likelihood
  if (retrieval_ptr->config->multinest_print_iter_values)
    std::cout << new_log_like << "\n";


  lnew = new_log_like;
}






//The dumper routine from MultiNest, will be called every updInt*10 iterations
//At the moment it's empty
//If you need to print/test stuff during a MultiNest run, add it here
void Retrieval::multinestDumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
  
}



}
