import pybear
from src.setup_transmission_model import BeARTransmissionRetrieval
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
obs_wavelength_bins = observation[:,0:2]
obs_data = observation[:,2]
obs_error = observation[:,3]

wasp_12b_obs = pybear.ObservationInput()

wasp_12b_obs.name = "WASP-12b_Kreidberg"
wasp_12b_obs.type = "band-spectroscopy"
wasp_12b_obs.bin_wavelength_edges = obs_wavelength_bins
wasp_12b_obs.data = obs_data
wasp_12b_obs.data_error = obs_error

observations = list([wasp_12b_obs])

priors = list([
  pybear.PriorConfig("gaussian", "log_g", [2.99, 0.03]),
  pybear.PriorConfig("uniform", "r_planet", [1.6, 2.3], "Rj"),
  pybear.PriorConfig("gaussian", "r_star", [1.57, 0.07], "Rs"),
  pybear.PriorConfig("log_uniform", "H2O", [1e-12, 1e-1]),
  pybear.PriorConfig("log_uniform", "TiO", [1e-12, 1e-3]),
  pybear.PriorConfig("log_uniform", "VO", [1e-12, 1e-3]),
  pybear.PriorConfig("log_uniform", "K", [1e-12, 1e-3]),
  pybear.PriorConfig("uniform", "temperature", [500, 1600])])

#create the BeAR forward model
bear_transmission = BeARTransmissionRetrieval(
  use_gpu,
  retrieval_folder,
  spectral_discretisation,
  resolution,
  cross_section_file_path,
  observations,
  priors)



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
