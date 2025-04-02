import os
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
parent_directory = os.path.dirname(current_directory)
sys.path.append(parent_directory)
sys.path.append(os.path.dirname(parent_directory))

from lib import pybear
import numpy as np


valid_spectral_discretisations = {'const_wavenumber', 'const_wavelength', 'const_resolution'}


class BeAROccultationModel:

  def __init__(
      self,
      use_gpu,
      grid_points_number,
      spectral_discretisation,
      wavelength_min,
      wavelength_max,
      resolution,
      cross_section_file_path, 
      opacity_species_data,
      stellar_spectrum_wavelengths,
      stellar_spectrum_flux,
      wavenumber_file_path = None) :
    
    if spectral_discretisation not in valid_spectral_discretisations:
      raise ValueError("Spectral discretisation must be one of %r." % valid_spectral_discretisations)
    
    self.nb_grid_points = grid_points_number
    
    bear_config = pybear.Config()
    
    bear_config.use_gpu = np.bool_(use_gpu)
    bear_config.forward_model_type = "secondary_eclipse"
    bear_config.cross_section_file_path = cross_section_file_path

    if wavenumber_file_path is not None:
      bear_config.wavenumber_file_path = wavenumber_file_path
    else :
      bear_config.wavenumber_file_path = ""
    
    if spectral_discretisation == 'const_wavenumber':
      bear_config.spectral_disecretisation = 0

    if spectral_discretisation == 'const_wavelength':
      bear_config.spectral_disecretisation = 1

    if spectral_discretisation == 'const_resolution':
      bear_config.spectral_disecretisation = 2
    
    bear_config.spectral_resolution = resolution
    
   
    self.spectral_grid = pybear.SpectralGrid(
      bear_config,
      wavelength_min,
      wavelength_max)
    
    self.wavelengths = np.flip(np.array(self.spectral_grid.wavelength_list))
    self.wavenumbers = np.flip(np.array(self.spectral_grid.wavenumber_list))

    opacity_species = opacity_species_data[:, 0]
    opacity_folders = opacity_species_data[:, 1]

    self.forward_model = pybear.OccultationModel(
      bear_config, 
      self.spectral_grid, 
      self.nb_grid_points, 
      stellar_spectrum_wavelengths,
      stellar_spectrum_flux,
      opacity_species, 
      opacity_folders)
  

  def calcSpectrum(
    self, 
    surface_gravity, 
    radius_ratio, 
    pressure, 
    temperature, 
    chem_species, 
    mixing_ratios,
    cloud_optical_depth) :

    cloud_tau = np.copy(cloud_optical_depth)
    
    #reverse the cloud optical depth array because BeAR uses wavenumbers in increasing order
    for i in range(cloud_tau.shape[0]):
      cloud_tau[i] = np.flip(cloud_tau[i])

    spectrum = np.array(
      self.forward_model.calcSpectrum(
        surface_gravity, 
        radius_ratio, 
        pressure, 
        temperature, 
        chem_species, 
        mixing_ratios, 
        cloud_tau))
    
    spectrum = np.flip(spectrum)

    return spectrum
