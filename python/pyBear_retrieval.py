from src.setup_bear_model import BeARTransmissionRetrieval
import numpy as np
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt
import pymultinest


#setting the basic properties of the model
spectral_discretisation = 'const_wavenumber'
resolution = 1.0
cross_section_file_path = "/media/data/opacity_data/helios-k/"
use_gpu = True
retrieval_folder = "../test_trans/"

#create the BeAR forward model
bear_transmission = BeARTransmissionRetrieval(
  use_gpu,
  retrieval_folder,
  spectral_discretisation,
  resolution,
  cross_section_file_path)


cube_parameters = np.zeros(bear_transmission.retrieval.nbParameters())
parameters = np.zeros(bear_transmission.retrieval.nbParameters())
physical_parameters = np.zeros(bear_transmission.retrieval.nbParameters())


def priors(cube, ndim, nparams):
  for i in range(nparams):
    cube_parameters[i] = cube[i]

  converted_values = bear_transmission.retrieval.convertCubeParameters(
    cube_parameters)
  
  global parameters
  parameters = converted_values[0]

  global physical_parameters
  physical_parameters = converted_values[1]

  for i in range(nparams):
    cube[i] = parameters[i]
  
  return None


def loglike(cube, ndim, nparams):
  
  new_log_like = bear_transmission.retrieval.computeLikelihood(physical_parameters)
  
  return new_log_like


print("Starting retrieval")
pymultinest.run(
  loglike, 
  priors, 
  bear_transmission.retrieval.nbParameters(), 
	resume = False, 
  verbose = True, 
  outputfiles_basename=retrieval_folder)


print("Starting post processing")
bear_transmission.postProcess()

