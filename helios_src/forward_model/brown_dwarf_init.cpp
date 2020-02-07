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


#include "brown_dwarf.h"


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>


#include "../additional/exceptions.h"
#include "../retrieval/retrieval.h"
#include "../chemistry/isoprofile_chemistry.h"
#include "../chemistry/fastchem_chemistry.h"
#include "../temperature/piecewise_poly_temperature.h"


namespace helios{


//initialise radiative transfer model
void BrownDwarfModel::initRadiativeTransfer(const BrownDwarfConfig& model_config)
{
  if (model_config.radiative_transfer_model == 0)
  {
    ShortCharacteristics* scm = new ShortCharacteristics(&retrieval->spectral_grid);
    radiative_transfer = scm; 
  }

  if (model_config.radiative_transfer_model == 1)
  {
    if (retrieval->config->use_gpu)
    {
      std::string error_message = "Radiative transfer model CDISORT cannot run on the GPU\n";
      throw ExceptionInvalidInput(std::string ("BrownDwarfModel::BrownDwarfModel"), error_message);
    }

    DiscreteOrdinates* disort = new DiscreteOrdinates(&retrieval->spectral_grid, 4, nb_grid_points); 
    radiative_transfer = disort;
  }

}


//select and initialise the chemistry models
void BrownDwarfModel::initChemistry(const BrownDwarfConfig& model_config)
{
  chemistry.assign(model_config.chemistry_model.size(), nullptr);


  for (size_t i=0; i<model_config.chemistry_model.size(); ++i)
  {
    if (model_config.chemistry_model[i] == 0)
    {
      IsoprofileChemistry* model = new IsoprofileChemistry(model_config.chemistry_parameters[i]);
      chemistry[i] = model;
    }


    if (model_config.chemistry_model[i] == 1)
    {
      FastChemChemistry* model = new FastChemChemistry(retrieval->config->retrieval_folder_path + model_config.chemistry_parameters[i][0], retrieval->config->nb_omp_processes);
      chemistry[i] = model;
    }
  }

  
  nb_total_chemistry_param = 0;

  for (auto & i : chemistry)
    nb_total_chemistry_param += i->nbParameters(); 
}



//select and initialise the chemistry models
void BrownDwarfModel::initTemperature(const BrownDwarfConfig& model_config)
{

  PiecewisePolynomialTemperature* temp = new PiecewisePolynomialTemperature(model_config.nb_temperature_elements, model_config.temperature_poly_degree, model_config.atmos_boundaries);
  temperature_profile = temp; 
  
}



}

