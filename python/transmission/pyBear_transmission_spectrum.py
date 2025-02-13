from src.setup_transmission_model import BeARTransmissionModel
import numpy as np
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt


#setting the basic properties of the model
spectral_discretisation = 'const_resolution'
wavelength_min = 0.4
wavelength_max = 10.0
resolution = 1000.0

cross_section_file_path = "/media/data/opacity_data/helios-k/"

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


use_gpu = True
grid_points_number = 100

#create the BeAR forward model
transmission_model = BeARTransmissionModel(
  use_gpu,
  grid_points_number,
  spectral_discretisation,
  wavelength_min,
  wavelength_max,
  resolution,
  cross_section_file_path, 
  opacity_species_data)


#define the planet/atmosphere parameters
surface_gravity = 2000.0  #in cm/s^2
planet_radius = 1.0 * const.R_jup.cgs.value
radius_ratio = planet_radius / const.R_sun.cgs.value

#p-T structure
temperature = np.full((grid_points_number), 500.0)
pressure = np.logspace(0, -7, grid_points_number)

#chemical composition
chem_species = np.array(['H2', 'He', 'H2O', 'CO', 'CO2'])
chem_species_concentration = np.array([0.0, 0.1, 1e-3, 1e-3, 1e-5])

chem_species_concentration[0] = 1.0 - np.sum(chem_species_concentration[1:])

#set the atmospheric mixing ratios
mixing_ratios = np.zeros((grid_points_number, chem_species.size))

for i in range(chem_species.size):
  mixing_ratios[:, i] = chem_species_concentration[i]

#cloud properties
#note: cloud optical depths are layer quantities, so the number of layers is grid_points_number-1
cloud_optical_depth = np.zeros((grid_points_number-1, transmission_model.wavelengths.size))

spectrum = transmission_model.calcSpectrum(
  surface_gravity, 
  planet_radius, 
  radius_ratio, 
  pressure, 
  temperature, 
  chem_species, 
  mixing_ratios,
  cloud_optical_depth,
  use_variable_gravity=True)

#adding a grey cloud layer with a constant optical depth at layer 30
cloud_optical_depth[30,:] = 100.0

spectrum2 = transmission_model.calcSpectrum(
  surface_gravity, 
  planet_radius, 
  radius_ratio, 
  pressure, 
  temperature, 
  chem_species, 
  mixing_ratios,
  cloud_optical_depth)


cloud_optical_depth = np.zeros((grid_points_number-1, transmission_model.wavelengths.size))
chem_species_concentration = np.array([0.0, 0.0, 0.3, 0.3, 0.3])
mixing_ratios = np.zeros((grid_points_number, chem_species.size))

for i in range(chem_species.size):
  mixing_ratios[:, i] = chem_species_concentration[i]

#temperature = np.full((grid_points_number), 1500.0)

spectrum3 = transmission_model.calcSpectrum(
  surface_gravity, 
  planet_radius, 
  radius_ratio, 
  pressure, 
  temperature, 
  chem_species, 
  mixing_ratios,
  cloud_optical_depth)


fig, ax = plt.subplots()
ax.plot(transmission_model.wavelengths, spectrum)
ax.plot(transmission_model.wavelengths, spectrum2)
ax.plot(transmission_model.wavelengths, spectrum3)
plt.show()