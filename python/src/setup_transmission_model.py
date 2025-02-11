
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
      observations,
      priors,
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


    bottom_pressure = 10
    top_pressure = 1e-5

    opacity_species_data = np.array([
      ['CIA-H2-H2', 'CIA/H2-H2'], 
      ['CIA-H2-He', 'CIA/H2-He'],
      ['H2', 'Rayleigh'],
      ['He', 'Rayleigh'],
      ['H2O', 'Molecules/1H2-16O__POKAZATEL_e2b_n'],
      ['CO', 'Molecules/12C-16O__Li2015_e2b_n'],
      ['CO2', 'Molecules/12C-16O2__UCL-4000_e2b'],
      ['H2O', 'Rayleigh'],
      ['CO', 'Rayleigh'],
      ['CO2', 'Rayleigh']])
    
    fit_mean_molecular_weight = False
    fit_scale_height = False
    use_variable_gravity = False
    nb_grid_points = 100
    atmos_bottom_pressure = 10
    atmos_top_pressure = 1e-5
    temperature_profile_model = "const"
    temperature_profile_parameters = []
    chemistry_model = ["iso", "bg"]
    chemistry_parameters = [["H2O", "TiO", "VO", "K"], ["H2He"]]
    cloud_model = []
    cloud_model_parameters = []
    modules = []
    modules_parameters = []
    opacity_species_symbol = opacity_species_data[:, 0]
    opacity_species_folder = opacity_species_data[:, 1]

    self.model_config = pybear.TransmissionModelConfig(
      fit_mean_molecular_weight,
      fit_scale_height,
      use_variable_gravity,
      nb_grid_points,
      atmos_bottom_pressure,
      atmos_top_pressure,
      temperature_profile_model,
      temperature_profile_parameters,
      chemistry_model,
      chemistry_parameters,
      cloud_model,
      cloud_model_parameters,
      modules,
      modules_parameters,
      opacity_species_symbol,
      opacity_species_folder)
    
   
    #self.retrieval = pybear.Retrieval(self.bear_config)
    self.retrieval = pybear.Retrieval(
      self.bear_config, 
      self.model_config, 
      observations,
      priors)

  def postProcess(self):
     self.post_process = pybear.PostProcess(self.bear_config)
     self.post_process.run()
