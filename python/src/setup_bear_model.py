
import pybear
import numpy as np


valid_spectral_discretisations = {'const_wavenumber', 'const_wavelength', 'const_resolution'}


class BeARTransmissionRetrieval:

  def __init__(
      self,
      use_gpu,
      retrieval_folder,
      spectral_discretisation,
      resolution,
      cross_section_file_path, 
      wavenumber_file_path = None) :
    
    self.bear_config = pybear.GlobalConfig()
    
    self.bear_config.use_gpu = np.bool_(use_gpu)
    self.bear_config.forward_model_type = "transmission"
    self.bear_config.cross_section_file_path = cross_section_file_path

    self.bear_config.retrieval_folder_path = retrieval_folder

    self.bear_config.multinest_print_iter_values = False

    if wavenumber_file_path is not None:
      self.bear_config.wavenumber_file_path = wavenumber_file_path
    else :
      self.bear_config.wavenumber_file_path = ""
    
    if spectral_discretisation == 'const_wavenumber':
      self.bear_config.spectral_disecretisation = 0
      self.bear_config.const_wavenumber_step = resolution

    if spectral_discretisation == 'const_wavelength':
      self.bear_config.spectral_disecretisation = 1
      self.bear_config.const_wavelength_step = resolution

    if spectral_discretisation == 'const_resolution':
      self.bear_config.spectral_disecretisation = 2
      self.bear_config.const_spectral_resolution = resolution
   
    self.retrieval = pybear.Retrieval(self.bear_config)

  def postProcess(self):
     self.post_process = pybear.PostProcess(self.bear_config)
     self.post_process.run()


class BeARTransmissionModel:

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
      wavenumber_file_path = None) :
    
    if spectral_discretisation not in valid_spectral_discretisations:
      raise ValueError("Spectral discretisation must be one of %r." % valid_spectral_discretisations)
    
    self.nb_grid_points = grid_points_number
    
    bear_config = pybear.GlobalConfig()
    
    bear_config.use_gpu = np.bool_(use_gpu)
    bear_config.forward_model_type = "transmission"
    bear_config.cross_section_file_path = cross_section_file_path

    if wavenumber_file_path is not None:
      bear_config.wavenumber_file_path = wavenumber_file_path
    else :
      bear_config.wavenumber_file_path = ""
    
    if spectral_discretisation == 'const_wavenumber':
      bear_config.spectral_disecretisation = 0
      bear_config.const_wavenumber_step = resolution

    if spectral_discretisation == 'const_wavelength':
      bear_config.spectral_disecretisation = 1
      bear_config.const_wavelength_step = resolution

    if spectral_discretisation == 'const_resolution':
      bear_config.spectral_disecretisation = 2
      bear_config.const_spectral_resolution = resolution
   
    self.spectral_grid = pybear.SpectralGrid(
      bear_config,
      wavelength_min,
      wavelength_max)
    
    self.wavelengths = np.flip(np.array(self.spectral_grid.wavelength_list))
    self.wavenumbers = np.flip(np.array(self.spectral_grid.wavenumber_list))

    opacity_species = opacity_species_data[:, 0]
    opacity_folders = opacity_species_data[:, 1]

    self.forward_model = pybear.TransmissionModel(
      bear_config, 
      self.spectral_grid, 
      self.nb_grid_points, 
      opacity_species, 
      opacity_folders)
  

  def calcSpectrum(
    self, 
    surface_gravity, 
    planet_radius, 
    radius_ratio, 
    pressure, 
    temperature, 
    chem_species, 
    mixing_ratios,
    cloud_optical_depth,
    use_variable_gravity=False) :

    cloud_tau = np.copy(cloud_optical_depth)
    
    #reverse the cloud optical depth array because BeAR uses wavenumbers in increasing order
    for i in range(cloud_tau.shape[0]):
      cloud_tau[i] = np.flip(cloud_tau[i])

    spectrum = np.array(
      self.forward_model.calcSpectrum(
        surface_gravity, 
        planet_radius, 
        radius_ratio, 
        pressure, 
        temperature, 
        chem_species, 
        mixing_ratios, 
        cloud_tau,
        np.bool_(use_variable_gravity)))
    
    spectrum = np.flip(spectrum)

    return spectrum



class BeARSecondaryEclipseModel:

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
    
    bear_config = pybear.GlobalConfig()
    
    bear_config.use_gpu = np.bool_(use_gpu)
    bear_config.forward_model_type = "secondary_eclipse"
    bear_config.cross_section_file_path = cross_section_file_path

    if wavenumber_file_path is not None:
      bear_config.wavenumber_file_path = wavenumber_file_path
    else :
      bear_config.wavenumber_file_path = ""
    
    if spectral_discretisation == 'const_wavenumber':
      bear_config.spectral_disecretisation = 0
      bear_config.const_wavenumber_step = resolution

    if spectral_discretisation == 'const_wavelength':
      bear_config.spectral_disecretisation = 1
      bear_config.const_wavelength_step = resolution

    if spectral_discretisation == 'const_resolution':
      bear_config.spectral_disecretisation = 2
      bear_config.const_spectral_resolution = resolution
      
    
   
    self.spectral_grid = pybear.SpectralGrid(
      bear_config,
      wavelength_min,
      wavelength_max)
    
    self.wavelengths = np.flip(np.array(self.spectral_grid.wavelength_list))
    self.wavenumbers = np.flip(np.array(self.spectral_grid.wavenumber_list))

    opacity_species = opacity_species_data[:, 0]
    opacity_folders = opacity_species_data[:, 1]

    self.forward_model = pybear.SecondaryEclipseModel(
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
