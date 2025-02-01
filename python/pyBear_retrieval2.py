import pybear
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


observation = np.loadtxt("../test_trans/WASP-12b_kreidberg.dat", skiprows=11)
obs_data = observation[:,2]
obs_error = observation[:,3]


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


def loglike2(cube, ndim, nparams):

  output = bear_transmission.retrieval.computeModel(physical_parameters)
  
  model_spectrum_obs = np.array(output.spectrum_obs[0])

  obs_delta = obs_data - model_spectrum_obs

  log_like = 0

  for i in range(obs_data.size):
    error_square = obs_error[i]**2
    log_like += (- 0.5 * np.log(error_square* 2.0 * np.pi)
                 - 0.5 * obs_delta[i]**2 / error_square)

  return log_like


print("Starting retrieval")
pymultinest.run(
  loglike2, 
  priors, 
  bear_transmission.retrieval.nbParameters(), 
	resume = False, 
  verbose = True, 
  outputfiles_basename=retrieval_folder)


print("Starting post processing")
bear_transmission.postProcess()

