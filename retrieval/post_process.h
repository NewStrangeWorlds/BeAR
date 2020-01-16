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


#ifndef _post_process_h
#define _post_process_h

#include <vector>
#include <string>


#include "../forward_model/forward_model.h"
#include "../spectral_grid/spectral_grid.h"
#include "../observations/observations.h"
#include "../config/global_config.h"
#include "multinest_parameter.h"
#include "retrieval.h"



namespace helios {

//forward declaration
class ForwardModel;


//the class that is responsible for the postprocessing
//it's derived from the retrieval class and inherits most of its capabilities
class PostProcess : public Retrieval{
  public:
    PostProcess(GlobalConfig* global_config);
    ~PostProcess();
    
    virtual void doRetrieval();

  private:
    ForwardModel* forward_model;

    std::vector< std::vector<double> > model_parameter;      //the values of the posteriors
    std::vector<double> log_like;                            //and their likelihood values
 
    size_t best_fit_model = 0;                               //best-fit model
    double best_log_like = 0;

    bool readPosteriorData();  
    void postProcessSpectra(std::vector< std::vector<double> >& model_spectrum_bands);
    void calcSpectrum(const unsigned int model_id, std::vector<double>& model_spectrum_bands);
    void calcSpectrumGPU(const unsigned int model_id, std::vector<double>& model_spectrum_bands);
    void postProcessEffectiveTemperatures(std::vector< std::vector<double> >& model_spectrum_bands, std::vector<double>& effective_temperatures);
    void postProcessTemperatures();

    void saveOutput(const std::vector< std::vector<double> >& model_spectrum_bands, const std::vector<double>& effective_temperatures);
    void saveBestFitSpectrum(const std::vector<double>& spectrum);
};


}


#endif
