from src.setup_bear_model import BeARSecondaryEclipseModel
import numpy as np
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from astropy.modeling import models
from astropy import units as u


#setting the basic properties of the model
spectral_discretisation = 'const_wavenumber'
wavelength_min = 1.0
wavelength_max = 10.0
resolution = 1.0

cross_section_file_path = "/media/data/opacity_data/helios-k/"

opacity_species_data = np.array([
  ['CIA-H2-H2', 'CIA/H2-H2'], 
  ['CIA-H2-He', 'CIA/H2-He'],
  ['H2O', 'Molecules/1H2-16O__POKAZATEL_e2b_n'],
  ['CO', 'Molecules/12C-16O__Li2015_e2b_n'],
  ['CO2', 'Molecules/12C-16O2__UCL-4000_e2b']
  ])


#define stellar spectrum, here we use a simple black body
stellar_spectrum_wavelengths = np.linspace(wavelength_min, wavelength_max, 1000)

bb_lam = models.BlackBody(temperature=4000*u.K, scale=1.0 * u.erg / (u.cm ** 2 * u.AA * u.s * u.sr))
stellar_spectrum_flux = bb_lam(stellar_spectrum_wavelengths * u.micron).to('W m-2 micron-1 sr-1').value * np.pi


use_gpu = True
grid_points_number = 100

#create the BeAR forward model
seconary_eclipse_model = BeARSecondaryEclipseModel(
  use_gpu,
  grid_points_number,
  spectral_discretisation,
  wavelength_min,
  wavelength_max,
  resolution,
  cross_section_file_path, 
  opacity_species_data,
  stellar_spectrum_wavelengths,
  stellar_spectrum_flux)


#define the planet/atmosphere parameters
surface_gravity = 2000.0  #in cm/s^2
planet_radius = 1.0 * const.R_jup.cgs.value
radius_ratio = planet_radius / const.R_sun.cgs.value


#create a fake T-p profile using a spline
temperature = np.full((grid_points_number), 0.0)
pressure = np.logspace(np.log10(200), -5, grid_points_number)

pressure_points = np.array([1e-5, 0.1, 1, 10, 200])
temperature_points = np.array([300, 800, 1000, 1500, 3000])

cs = CubicSpline(np.log10(pressure_points), temperature_points)
temperature = cs(np.log10(pressure))


#chemical composition
chem_species = np.array(['H2', 'He', 'H2O', 'CO', 'CO2'])
chem_species_concentration = np.array([0.0, 0.145, 1e-5, 1e-5, 1e-7])

chem_species_concentration[0] = 1.0 - np.sum(chem_species_concentration[1:])

#set the atmospheric mixing ratios
mixing_ratios = np.zeros((grid_points_number, chem_species.size))

for i in range(chem_species.size):
  mixing_ratios[:, i] = chem_species_concentration[i]

#cloud properties
#note: cloud optical depths are layer quantities, so the number of layers is grid_points_number-1
cloud_optical_depth = np.zeros((grid_points_number-1, seconary_eclipse_model.wavelengths.size))

spectrum = seconary_eclipse_model.calcSpectrum(
  surface_gravity, 
  radius_ratio, 
  pressure, 
  temperature, 
  chem_species, 
  mixing_ratios,
  cloud_optical_depth)

#adding a grey cloud layer with a constant optical depth at layer 30
cloud_optical_depth[30,:] = 100.0

spectrum2 = seconary_eclipse_model.calcSpectrum(
  surface_gravity, 
  radius_ratio, 
  pressure, 
  temperature, 
  chem_species, 
  mixing_ratios,
  cloud_optical_depth)


fig, ax = plt.subplots()
ax.plot(seconary_eclipse_model.wavelengths, spectrum)
ax.plot(seconary_eclipse_model.wavelengths, spectrum2)
plt.show()