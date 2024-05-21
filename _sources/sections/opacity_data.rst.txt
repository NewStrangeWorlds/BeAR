Opacity data
============

``BeAR`` uses tabulated absorption cross-sections for its radiative transfer calculations. The
cross-sections have to be given in units of cm2 .
You can of course use your own opacity data. If you choose to do that, you need either
convert your data to the format described below or change the corresponding code within
``BeAR``. More in-depth details about the implementation in the code can be found in Part
III.


Required format of the tabulated data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cross-sections are required to be tabulated as a function of wavenumber. All cross-sections
have to be given on the same wavenumber grid. ``BeAR`` expects to find this global
wavenumber grid in a special file. The location of the opacity data and the wavenumber
file are handled by the specific forward model. The configuration of the forward models will
usually have an option to specify specific folders where the data for the different species can
be found by ``BeAR``.
