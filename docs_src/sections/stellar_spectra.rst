
.. _sec:stellar_spectra:

Stellar spectra
===============

Some forward models or optional modules in BeAR require a stellar spectrum as input. BeAR can handle
stellar spectra in the three different ways:

  - :ref:`Blackbody spectra <sec:stellar_spectrum_bb>`

  - :ref:`A single, tabulated stellar spectrum <sec:stellar_spectrum_file>`
  
  -  :ref:`A grid of stellar atmosphere spectra <sec:stellar_spectrum_bb>`


.. _sec:stellar_spectrum_bb: 

Blackbody spectra
-----------------

The simplest way to provide a stellar spectrum is by using a blackbody spectrum. It is chosen by using
the keyword :code:`blackbody` in the corresponding section in the ``forward_model.config`` file:

.. code:: 

   #Stellar spectrum
   blackbody

This option requires the user to provide the temperature of the blackbody spectrum in Kelvin as a free
parameter in the priors configuration file.


.. _sec:stellar_spectrum_file:

Tabulated stellar spectrum
--------------------------

If a tabulated stellar spectrum is used, the keyword :code:`file` has to be used in the corresponding

`forward_model.config`` file:

.. code:: 

   #Stellar spectrum
   file WASP-43.dat

with the file name of the stellar spectrum following the keyword. The file location hast to be given relative 
to the path of the BeAR executable. It has the following structure:

.. include:: ../examples/WASP-43.dat
   :literal:

The first column contains the wavelengths in micrometers, while the second column contains the stellar photospheric
flux in :math:`\mathrm{W/m^2/\mu m}`. The spectrum will be interpolated to the wavenumber grid of the retrieval.
However, it needs to cover the entire wavelength range of the observations that are used in a retrieval. 
BeAR will not extrapolate the spectrum but replace missing values outside of the tabulated range with zeros.
This stellar spectrum does not require a free parameter in the priors configuration file.


.. _sec:stellar_spectrum_grid:

Grid of tabulated stellar spectra
---------------------------------

In some cases, it is necessary to use a grid of stellar atmosphere models. This is especially the case
when one wants to sample stellar spectra for a range of different effective temperatures, surface gravities, or metallicities.
Examples of grids of stellar spectra are the `PHOENIX stellar atmosphere library <https://phoenix.astro.physik.uni-goettingen.de/>`_ 
or the `SPHINX library <https://zenodo.org/records/7416042>`_

To use a grid of stellar spectra, the keyword :code:`grid` has to be used in the corresponding ``forward_model.config`` file:

.. code:: 

   #Stellar spectrum
   grid /data/phoenix_grid/

followed by a path to the folder that contains the grid of stellar atmosphere models. 
BeAR currently assumes that the grid is given in terms of:

  - stellar effective temperature in Kelvin

  - logarithm of the surface gravity in :math:`\mathrm{cm/s^2}`

  - metallicity in logarithm of solar units (typically [Fe/H])


The folder with the stellar grid should contain a file called ``grid_parameters.dat`` that 
lists the parameters of the grid. It needs to have the following structure:

.. include:: ../examples/stellar_grid_parameters.dat
   :literal:

The first line contains the list of stellar effective temperatures, the second line the list of surface gravities, 
and the third line the list of metallicities. BeAR assumes that the grid rectangular and regular, which means that
every possible combinations of the parameters should have a tabulated stellar spectrum.
If the original grid is not complete, the user has to extrapolate the grid to cover all possible combinations or 
create dummy data sets if those spectra will not get sampled during a retrieval.

The fourth line contains the name of the file that contains the wavelength grid for the stellar spectra. Each spectrum
has to be tabulated at these wavelengths. These wavelengths have to be in micrometers and saved in a binary file that
contains single precision floating point numbers.

Below the wavelength file, the actual grid of stellar spectra is listed. 
Each line contains the effective temperature, the surface gravity, the metallicity, and the file name of the tabulated spectrum.
Like the wavelength file, the spectra have to be saved in binary format and have to be tabulated at the wavelengths listed 
in the wavelength file. The spectra need to cover the entire wavelength range of the observations that are used in a retrieval. 
BeAR will not extrapolate the spectra but replace missing values outside of the tabulated range with zeros.

The binary files should contain the stellar photospheric flux in :math:`\mathrm{W/m^2/\mu m}` as single-precision 
floating point numbers. As mentioned above, BeAR expects to find a stellar spectrum entry for every possible combination
of stellar parameters. If this is not the case, an error message will be displayed and the code will terminate.

The use of the stellar grid requires the user to provide the effective temperature, the surface gravity, and the 
metallicity as free parameters in the priors configuration file.
