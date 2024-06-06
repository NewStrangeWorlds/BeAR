
Opacity data
============

.. _sec:opacity_data:

BeAR uses mainly tabulated absorption cross-sections for its radiative transfer calculations. 
The cross-sections have to be given in units of :math:`\mathrm{cm^{2}}` or in cross-sections
per species mass :math:`\mathrm{cm^{2}/g}`. BeAR is primarily designed to be used 
with the output files from the HELIOS-k opacity calculator. HELIOS-k is an open-source code 
that can be found `here <https://github.com/exoclime/HELIOS-K/>`_.

Pre-calculated opacity data for a broad range of species that have been curated by Simon Grimm 
can be found on the `DACE platform <https://dace.unige.ch//>`_ under ``Opacity``. 
This data is already in a format that BeAR can use directly.

While BeAR is designed to work primarily with the opacity data format of the HELIOS-k 
opacity calculator, one can of course also use other opacities. In that case, they have to 
be converted to the format described below.


Global wavenumber grid
......................

Cross-sections are required to be tabulated as a function of wavenumber. All cross-sections
have to be given on the same wavenumber grid to avoid having to constantly interpolating them
in wavenumber space when reading them in. By default, BeAR will use the HELIOS-k opacities. 
In that case, the opacity data is calculated on a wavenumber grid with a constant step of 
0.01 :math:`\mathrm{cm^{-1}}`, starting from 0 :math:`\mathrm{cm^{-1}}`.

BeAR expects to find this global wavenumber grid in a special file. The location of the opacity 
data and the wavenumber file are handled by the ``retrieval.config`` configuration file.

This file has the following, simple structure.

.. include:: ../examples/wavenumbers_full.dat_example
   :literal:
   
The file starts with an integer value that is equal to the total number of wavenumbers in
the file. After that, the single wavenumbers are tabulated in ascending order. 

   
Opacity data for a specific species
...................................

Opacity data for a specific species have to be located in a folder that is usually
specified in the opacity species list within the ``forward_model.config``. Inside
this folder, BeAR expects to find a summary file ``filelist.dat`` that contains
a list of pressures and temperatures as well as filenames the opacities have
been tabulated at.

In general it has the following structure:

.. include:: ../examples/filelist_water.dat
   :literal:
   
The first line is a header with important information. Internally, BeAR uses
cross-sections in units of :math:`\mathrm{cm^{2}}`, while the the standard
HELIOS-k data is usually given in :math:`\mathrm{cm^{2}/g}`, that means 
cross-sections per species mass.

For the conversion between the two, the molecular weight of the species is
required, which has to be placed within the header as the first value. In
the above example, this corresponds to the molecular weight of water.
In case the actual opacity data is already given in :math:`\mathrm{cm^{2}}`,
a value of ``0`` needs to be used here, instead.

Below the header is a list of all available opacity files together with
the pressure in bar and temperature in Kelvin. The first column refers to
the pressure, the second column to the temperature, while the third one
is the corresponding file name where the corresponding, tabulated opacity
data can be found.

All individual opacity files are saved in binary format to decrease the time 
it requires to read them in. They have to be tabulated at exactly the same
wavenumbers that are listed in ``wavenumber_full.dat`` file.

However, not all chemical species have absorption lines that cover the entire, 
global wavenumber range. Some, for example, might only absorb in the infrared but 
not in the UV. Therefore, the cross-sections don’t need to be tabulated over the 
entire global wavenumber range. Instead, they may stop once no absorption
lines are present any more. Internally, the missing cross-sections are set to zero. 
Note, however, that there must be no holes in the data files. 
If a molecule has holes in its absorption data, the corresponding cross-sections 
have to be stored as zeros in the binary file.

Inside the binary file itself, the data has to be stored as single-precision values,
not the usual double-precision ones. This is already the case for the standard
HELIOS-k opacity files that can be found on DACE, for example.


Special cases
-------------

Sometimes, opacity data can be too small to be represented by a single-precision
float value. This is often the case with collision-induced absorption (CIA). In this case,
the opacities can also be saved and read-in in log10 space.

In order to signal BeAR that the opacities are stored in log space, the header of the
summary file ``filelist.dat`` has to be adjusted slightly as shown in the example below.

.. include:: ../examples/filelist_h2h2.dat
   :literal:
   
In this example, cross-sections in :math:`\mathrm{cm^{2}}` are stored in log10 space.
Thus, the molecular weight is set to 0, and the parameter ``log`` is used afterwards.
If the data is stored in linear space, ``lin`` can be used. Since, however, this is
standard behaviour for HELIOS-k data, it can be omitted.

Another special case for opacities is pressure broadening by a special species. Usually, BeAR
treats the pressure the opacities are tabulated at as the total gas pressure. However,
sometimes absorption cross-sections have been calculated for a specific background perturber
species. This can also be taken into account in the header of the ``filelist.dat`` file as
shown in the example below.

.. include:: ../examples/filelist_K.dat
   :literal:

Here, cross-sections per mass for potassium are stored in linear space with molecular hydrogen as the
perturber species. Thus, instead of the total gas pressure, BeAR will use the partial
pressure of H2 for these opacities. The perturber species is listed in the third
column of the header. Note, that when this option is used, the second column with
the storage format of the opacities has to be present as well.





