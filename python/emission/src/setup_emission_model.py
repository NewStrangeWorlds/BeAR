import os
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
parent_directory = os.path.dirname(current_directory)
sys.path.append(parent_directory)
sys.path.append(os.path.dirname(parent_directory))

from lib import bear
import numpy as np


valid_spectral_discretisations = {'const_wavenumber', 'const_wavelength', 'const_resolution'}


class BeAREmissionModel:

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
      wavenumber_file_path = None,
      radiative_transfer_model = "scm",
      radiative_transfer_model_parameters = []) :
    
    if spectral_discretisation not in valid_spectral_discretisations:
      raise ValueError("Spectral discretisation must be one of %r." % valid_spectral_discretisations)
    
    self.nb_grid_points = grid_points_number
    
    bear_config = bear.Config()
    
    bear_config.use_gpu = np.bool_(use_gpu)
    bear_config.forward_model_type = "emission"
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
    
   
    self.spectral_grid = bear.SpectralGrid(
      bear_config,
      wavelength_min,
      wavelength_max)
    
    self.wavelengths = np.flip(np.array(self.spectral_grid.wavelength_list))
    self.wavenumbers = np.flip(np.array(self.spectral_grid.wavenumber_list))

    opacity_species = opacity_species_data[:, 0]
    opacity_folders = opacity_species_data[:, 1]

    self.forward_model = bear.EmissionModel(
      bear_config, 
      self.spectral_grid, 
      self.nb_grid_points, 
      radiative_transfer_model,
      radiative_transfer_model_parameters,
      opacity_species, 
      opacity_folders)
  

  def calcSpectrum(
    self, 
    surface_gravity, 
    radius,
    distance, 
    pressure, 
    temperature, 
    chem_species, 
    mixing_ratios,
    cloud_optical_depth,
    use_variable_gravity=False,
    cloud_single_scattering=None,
    cloud_asym_param=None) :

    #reverse the cloud arrays along the spectral axis because BeAR uses wavenumbers in increasing order
    #all cloud arrays are indexed [layer, wavelength]
    def reverse_spectral_axis(array):
      array = np.copy(array)
      for i in range(array.shape[0]):
        array[i] = np.flip(array[i])
      return array

    cloud_tau = reverse_spectral_axis(cloud_optical_depth)

    if cloud_single_scattering is None and cloud_asym_param is None:
      spectrum = np.array(
        self.forward_model.calcSpectrum(
          surface_gravity,
          radius,
          distance,
          pressure,
          temperature,
          chem_species,
          mixing_ratios,
          cloud_tau,
          use_variable_gravity))
    else:
      if cloud_single_scattering is None:
        cloud_single_scattering = np.zeros_like(cloud_tau)
      if cloud_asym_param is None:
        cloud_asym_param = np.zeros_like(cloud_tau)

      spectrum = np.array(
        self.forward_model.calcSpectrum(
          surface_gravity,
          radius,
          distance,
          pressure,
          temperature,
          chem_species,
          mixing_ratios,
          cloud_tau,
          reverse_spectral_axis(cloud_single_scattering),
          reverse_spectral_axis(cloud_asym_param),
          use_variable_gravity))
    
    spectrum = np.flip(spectrum)

    return spectrum
