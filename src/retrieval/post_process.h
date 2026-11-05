// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _post_process_h
#define _post_process_h

#include <vector>
#include <string>


#include "../spectral_grid/spectral_grid.h"
#include "../observations/observations.h"
#include "../config/global_config.h"
#include "retrieval.h"



namespace bear {

//forward declaration
class ForwardModel;


//the class that is responsible for the postprocessing
//it's derived from the retrieval class and inherits most of its capabilities
class PostProcess : public Retrieval{
  public:
    PostProcess(
      GlobalConfig* global_config);
    PostProcess(
      GlobalConfig* global_config,
      GenericConfig* model_config,
      GenericConfig* model_postprocess_config_,
      const std::vector<ObservationInput>& observation_input,
      const std::vector<PriorConfig>& prior_config);
    ~PostProcess();
    
    virtual bool run();
    bool run(
      const std::string posterior_file_path);
    bool run(
      std::vector<std::vector<double>>& posteriors, 
      std::vector<double>& log_likelihoods);
  private:
    GenericConfig* model_postprocess_config = nullptr;
    
    std::vector<std::vector<double>> model_parameter;      //the values of the posteriors
    std::vector<double> log_like;                            //and their likelihood values
 
    size_t best_fit_model = 0;                               //best-fit model
    double best_log_like = 0;                                //likelihood of best-fit model

    void readPosteriorData(const std::string file_path);
    void processPosteriorData();

    void postProcessSpectra(std::vector< std::vector<double> >& model_spectrum_bands);

    void deleteSamplerFiles(const std::vector<std::string>& file_list);
};


}


#endif
