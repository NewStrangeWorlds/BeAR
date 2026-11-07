#This Python script can be used to automatically generate the filelist.dat
#from DACE opacity files with the standard HELIOS-K naming convention for
#species and temperature & pressure
#The script will scan all folders in the folder main_directory

import numpy as np
import os
from pyteomics import mass
import re

main_directory = 'Molecules/'


print("Scanning directory", main_directory)

#scan the directory for all 'bin' files and extract pressures and temperatures
for directory in os.listdir(main_directory):
  if not os.path.isdir(main_directory + directory):
    continue
  
  molecular_formula = directory.split('__')[0]

  element_formula = molecular_formula.split('-')
  nb_elements = len(element_formula)
  
  print(directory, molecular_formula, element_formula)

  formula = ''
  skip = False

  for e in element_formula:
    #test = re.split('(\d+)',e)
    element_components = re.findall('\d+|\D+', e)
    
    if len(element_components) < 2:
      print("Error: could not parse element", e, "in formula", molecular_formula)
      skip = True
      break
    
    formula += element_components[1] + '[' + element_components[0] + ']'
    if (len(element_components) > 2):
      formula += element_components[2]
  
  if skip:
    continue
  formula = formula.replace("_p", "")

  print("Parsed formula:", formula)
  species_mass = mass.calculate_mass(formula = formula)
  
  print("Molecular weight for", molecular_formula, "is", species_mass)

  file_list = []

  for file_name in os.listdir(main_directory+directory):
    if file_name.endswith(".bin"):
        
      #print("found", file_name)
        
      pressure_string = file_name.split('.')[0].split("_")[-1]
      temperature_string = file_name.split('.')[0].split("_")[-2]
    
      pressure_sign = "+"
    
      if pressure_string[0] == 'n':
        pressure_sign = '-'

      temperature = float(temperature_string)
      log_pressure = float(pressure_sign + pressure_string[1] + '.' + pressure_string[2:] + pressure_string[2:] + pressure_string[2:] + pressure_string[2:])
      pressure = 10**log_pressure

      file_list.append((pressure, temperature, file_name))

  file_list.sort(key=lambda x: (x[0], x[1]))

  # for i in range(len(file_list)):
  #   print("Pressure:", '{0:>1.5e}'.format(file_list[i][0]), "Temperature:", '{0:>1.5e}'.format(file_list[i][1]), "File:", file_list[i][2])

  file = open(main_directory+directory+'/filelist.dat', 'w')

  file.write(str(species_mass)+'\n')

  for i in range(len(file_list)):
    file.write('{0:>1.5e}  {1:>1.5e}  {2:>s}\n'.format(file_list[i][0], file_list[i][1], file_list[i][2]))
  
  file.close()
