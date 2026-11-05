
.. _sec:opacity_data:

Opacity data
============

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

Besides the tabulated opacities, BeAR offers a variety of built-in opacities, namely Rayleigh
scattering and continuum absorption. The available data is 
discussed :ref:`here <sec:built-in-opacity_data>`.


Global wavenumber grid
......................

Cross-sections are required to be tabulated as a function of wavenumber. All cross-sections
have to be given on the same wavenumber grid to avoid having to constantly interpolating them
in wavenumber space when reading them in. By default, BeAR will use the HELIOS-k opacities. 
In that case, the opacity data is calculated on a wavenumber grid with a constant step of 
0.01 :math:`\mathrm{cm^{-1}}`, starting from 0 :math:`\mathrm{cm^{-1}}`. Unless indicated otherwise,
BeAR will implicitely assume that this standard HELIOS-k grid is used and no further configuration is needed.

However, if other opacity data is used, the wavenumber grid has to be provided as well. In that case,
BeAR expects to find this global wavenumber grid in a special file. The location of the opacity
data and the wavenumber file are handled by the ``retrieval.toml`` configuration file. The root
opacity folder is set via the ``opacity_data_folder`` key in the ``[retrieval]`` section of
``retrieval.toml``. Within this root opacity folder, BeAR expects to find a file
called ``wavenumber_full.dat``.

This file has the following, simple structure.

.. include:: ../examples/wavenumber_full.dat_example
   :literal:
   
The file starts with an integer value that is equal to the total number of wavenumbers in
the file. After that, the single wavenumbers are tabulated in ascending order. 

   
Opacity data for a specific species
...................................

Opacity data for a specific species have to be located in a folder that is usually
specified in the opacity species list within the ``forward_model.toml``. Inside
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
is the file name where the corresponding, tabulated opacity
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
the opacities can also be saved and read in in log10 space.

In order to signal BeAR that the opacities are stored in log space, the header of the
summary file ``filelist.dat`` has to be adjusted slightly as shown in the example below.

.. include:: ../examples/filelist_h2h2.dat
   :literal:
   
In this example, cross-sections in :math:`\mathrm{cm^{2}}` are stored in log10 space.
Thus, the molecular weight is set to 0, and the parameter ``log`` is used afterwards.
If the data is stored in linear space, ``lin`` can be used. Since, however, this is
standard behaviour for HELIOS-k data, it can be omitted as shown in the first example.

Another special case for opacities is pressure broadening by a specific perturber species.
BeAR always treats the pressure the opacities are tabulated at as the total gas pressure.
Cross-sections that have been computed for a specific background perturber (e.g. H2) are
therefore interpolated at the total pressure rather than at the partial pressure of the
perturber. A header line such as

.. include:: ../examples/filelist_K.dat
   :literal:

with the perturber species in a third column is accepted (the additional column is ignored by
the parser), but the partial-pressure treatment is currently not implemented.

.. _sec:built-in-opacity_data:

Opacities included in BeAR
..........................

BeAR includes a selection of different opacity sources that can be calculated on-the-fly and do not
need to be read in as tabulated data. This currently includes Rayleigh scattering by

  - molecular hydrogen, H2

  - atomic hydrogen, H

  - atomic helium, He

  - water, H2O

  - carbon monoxide, CO

  - carbon dioxide, CO2

  - methane, CH4

Additionally, BeAR can calculate the bound-free and free-free continuum absorption for
the hydrogen anion H-. This, however, is only relevant for very hot scenarios as 
otherwise the abundance of H- would not be high enough to have a strong impact on
the overall opacity.

If additional Rayleigh scattering or more continuum absorbers are needed, they can be
either added to the code or simply be pre-tabulated and then used as the regular opacity
data discussed above.


.. _sec:opacity_config:

Choosing opacity sources in forward models
..........................................

Most forward models in BeAR require the user to specify the opacity sources that should be used.
These are given in the ``forward_model.toml`` file as a TOML array of tables assigned to the
``opacity`` key. The general format for this is as follows:

.. include:: ../examples/opacity_sources1.toml
   :literal:

Each entry has a ``species`` field with the species' chemical formula and a ``folder`` field.
The species names need to be listed in the internal
chemical species list of BeAR as discussed :ref:`here <sec:chemistry_models>`. The ``folder`` field
indicates the folder where the opacity data, or more specifically the summary file ``filelist.dat``,
for the species can be found. The folders here are relative to the root opacity folder that is set
via the ``opacity_data_folder`` key in the ``[retrieval]`` section of ``retrieval.toml``.
The species do not have to appear in any specific order.

For the Rayleigh-scattering species mentioned above, the folder name is replaced by the keyword
``Rayleigh`` as shown in the following example, where Rayleigh scattering for H2 and He is added.

.. include:: ../examples/opacity_sources2.toml
   :literal:

Special cases are the bound-free and free-free continuum absorptions for H- and collision-induced
absorption (CIA). These are added by using the keyword ``H-`` and the ``CIA-`` species names,
respectively, as shown below:

.. include:: ../examples/opacity_sources3.toml
   :literal:

BeAR currently can add the following CIA sources:

  - ``CIA-H2-H2`` - CIA by H2 and H2

  - ``CIA-H2-He`` - CIA by H2 and He

  - ``CIA-H-He`` - CIA by atomic H and He

The opacity for the hydrogen anion H- is calculated by BeAR on-the-fly and does not need a folder
to be specified. Thus, the keyword ``none`` is used here instead.

.. _sec:opacity_highres:

Separate opacity sources for the high-resolution model
......................................................

The phase-curve (high-resolution) forward model can read a *separate* opacity table that is used
only for its high-resolution part. This is done via an optional ``opacity_highres`` key in
``forward_model.toml``, which uses exactly the same array-of-tables format as the normal
``opacity`` block described above:

.. code::

   opacity_highres = [
     { species = "CO", folder = "CO/high_res" },
     { species = "H2O", folder = "H2O/high_res" }
   ]

Each entry again has a ``species`` field with the chemical formula and a ``folder`` field pointing to
the folder that contains the corresponding ``filelist.dat``, relative to the root opacity folder set
via the ``opacity_data_folder`` key in the ``[retrieval]`` section of ``retrieval.toml``. All of the
special keywords (``Rayleigh``, ``H-``, the ``CIA-`` sources) are supported here as well.

The purpose of this separate list is to let the high-resolution part of the model use different
opacity data than the low-resolution part, typically cross-sections computed at a higher spectral
resolution or from a different line list, while keeping the (usually coarser) low-resolution
opacities in the normal ``opacity`` block.

If ``opacity_highres`` is absent or empty, the high-resolution model automatically falls back to the
species and folders given in the normal ``opacity`` list, so specifying it is only necessary when the
two parts of the model should use genuinely different opacity data.

See :ref:`sec:high_resolution` for how the high-resolution opacities enter the re-injection and
cross-correlation calculation.