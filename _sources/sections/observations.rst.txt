Observations
============

Supported observational types
-----------------------------

BeAR currently supports three different types of observations:

- :ref:`Spectroscopy <sec:spectroscopy>`

- :ref:`Band-Spectroscopy <sec:band_spectroscopy>`

- :ref:`Photometry <sec:photometry>`

Based on the type of observation, the required format of the data files
differs slightly. In the following sections, the three basic types and
their required input formats are described. Observations have to be
provided in wavelength space. Internally, however, BeAR performs the
calculations in wavenumber space.

Typically, the computions are done with a higher spectral resolution than
the observations and then integrated to the observational wavelength structure.
The resolution is determined by the corresponding configuration parameter in
the retrieval.config file.

For spectroscopy and band spectroscopy obervations, BeAR has the option to
use :ref:`instrument line profiles <sec:instrument_line_profile>`, which typically
spreads the flux at a given wavelength over several adjacent pixels.

Additionally, BeAR can use :ref:`filter transmission functions <sec:filter_response>`. 
While this is typically only used in photometric observations, BeAR also supports 
them for both spectroscopy and band spectroscopy.


.. _sec:spectroscopy:

Spectroscopy
------------

In spectroscopy mode, an observational spectrum is given at specific
wavelengths :math:`\lambda_i`, from :math:`\lambda_1` to :math:`\lambda_N`.
The internal respresentation in wavenumber space :math:`\nu`
is shown in the figure below.

.. image:: ../images/spectroscopy.png
  :width: 700
  :align: center

BeAR sets up a high-resolution wavenumber grid, with a step size determined 
by the corresponding option in the retrieval.config file. This grid is
symbolised by the red lines in the figure.

To simulate the observed flux at the given wavelengths, BeAR creates a 
structure composed of spectral bands, one band for each observational wavelength. 
The boundaries of these bands in wavenumber space are halfway between adjacent
wavelengths. For example, the boundaries :math:`\nu_{i,1}` and :math:`\nu_{i,2}` 
for the i-th band, corresponding to the wavelength :math:`\lambda_i`, are
determined by the adjacent wavelengths :math:`\lambda_{i-1}` and :math:`\lambda_{i+1}`.
BeAR calculates its model spectrum on the high-resolution wavenumber grid (the red lines). 
It then obtains the mean flux in each of the bands i via integration and identifies
the result with the flux at the observational wavelengths :math:`\lambda_i`. 


Input file structure
....................

A basic example for an input file for a spectroscopic observation is shown below.

.. include:: ../examples/gj570d_spex_min.dat
   :literal:
   
The file consists of a header that contains some basic information. The name of the 
observation/instrument is not used during the calculation but will determine the
file name of the posterior spectra file. 

For spectroscopy, the ``#type`` needs to be set to ``spectroscopy``. This is followed by
an optional :ref:`filter bandpass transmission function <sec:filter_response>`.
If no filter transmission is used, this should be set to ``none`` as in the example above.

The actual spectroscopic data is given in three columns. The first column is the wavelength
in units of :math:`\mathrm{\mu m}`, the second the observational data. The units of the data depend on the
chosen forward model. For example, the ``emission`` forward model expects a radiation flux
in units of :math:`\mathrm{W} \mathrm{m}^{-2}  \mathrm{\mu m}^{-1}`, while the ``transmission`` spectroscopy model requires
the transit depth in ppm. The third column contains the error of the observational data in
the same units as the previous column.


As mentionend above, an optional Gaussian :ref:`instrument line profiles <sec:instrument_line_profile>`,
characterised by its FWHM, can be used in BeAR. This information is added in an optional fourth column 
as shown below.

.. include:: ../examples/gj570d_spex.dat
   :literal:
   
The fourth column contains the FWHM of the Gaussian profile in :math:`\mathrm{\mu m}`. Setting the FWHM to
0 will result in the instrument line profile being neglected. Another optional fifth column contains a 
weighting factor for each observational point. This allows to give unreliable data points a lower impact during the
computation of the likelihood or to neglect certain points entirely.


.. _sec:band_spectroscopy:

Band Spectroscopy
-----------------

Band-spectroscopy is a degraded form of spectroscopy, where individual wavelengths have
been summed up into bands to e.g. increase the signal-to-noise of a low-signal observation.
This is, for example, commonly done for exoplanet observations with the WFC3 instrument
on the Hubble Space Telescope. The band structure itself does not need to be regular.


As depicted in the figure above, the observational data is assumed to consist of
math:`i = 1 ... N` spectral bands, each with given wavelength boundaries :math:`\lambda_{i,1}` 
and :math:`\lambda_{i,2}` . 

.. image:: ../images/band_spectroscopy.png
  :width: 700
  :align: center

BeAR will create the same band structure in wavenumber space. Just like for spectroscopy calculations,
the high-resolution spectrum of BeAR will be integrated over each band math:`i` to obtain the
mean flux of the corresponding observation. Optionally, before integrating the spectrum, the
high-resolution spectrum can be convolved with a given instrument line profile to simulate the
flux received by the detector.


Input file structure
....................

A basic example for an input file for a band spectroscopy observation is shown below.

.. include:: ../examples/WASP-12b_kreidberg_min.dat
   :literal:

The file consists of a header that contains some basic information. The name of the 
observation/instrument is not used during the calculation but will determine the
file name of the posterior spectra file. 
For band spectroscopy, the ``#type`` needs to be set to ``band-spectroscopy``.
This is followed by an optional :ref:`filter bandpass transmission function <sec:filter_response>`.
If no filter transmission is used, this should be set to ``none`` as in the example above.

The actual spectroscopic data is given in four columns. The first two columns describe the boundaries 
of the wavelength bins in units of :math:`\mathrm{\mu m}`. 
The third column refers to the observational data. The units of the data depend on the
chosen forward model. For example, the ``emission`` forward model expects a radiation flux
in units of :math:`\mathrm{W} \mathrm{m}^{-2}  \mathrm{\mu m}^{-1}`, while the ``transmission`` spectroscopy model requires
the transit depth in ppm as shown in the example above.

The fourth column contains the error of the observational data in the same units as the previous column.

Just like spectroscopic data, an optional Gaussian :ref:`instrument line profiles <sec:instrument_line_profile>`, 
characterised by its FWHM, can be used for band spectroscopy as well. This information is added 
in an optional fifth column as shown below.

.. include:: ../examples/WASP-12b_kreidberg.dat
   :literal:

The additional column contains the FWHM of the Gaussian profile in :math:`\mathrm{\mu m}`. Setting the FWHM to
0 will result in the instrument line profile being neglected as shown in the example above. 
Another optional sixth column contains a  weighting factor for each observational band. 
This allows to give unreliable data points a lower impact during thecomputation of the likelihood 
or to neglect certain points entirely.


.. _sec:photometry:

Photometry
----------

Photometry is essentially band-spectroscopy with just one broad band between two wavelengths 
:math:`\lambda_{1}` and :math:`\lambda_{2}` as depicted in the figure below.

.. image:: ../images/photometry.png
  :width: 700
  :align: center

The high-resolution spectrum calculated by BeAR will
be integrated over the bandpass in wavenumber space to obtain the mean flux in the filter.
The conceptual difference between band-spectroscopy and photometry within BeAR
is that unlike the former, a photometry observation does not have an instrument line profile
because it’s supposed to cover a broader wavelength range. Instead, it can be processed
through a :ref:`filter transmission function <sec:filter_response>` to simulate the observation 
through a specific filter.


Input file structure
....................

A basic example for an input file for a photometric observation is shown below.

.. include:: ../examples/wasp-43b_spitzer_2_min.dat
   :literal:
   
The file consists of a header that contains some basic information. The name of the 
observation/instrument is not used during the calculation but will determine the
file name of the posterior spectra file. 
For photometry, the ``#type`` needs to be set to ``photometry``.

This is followed by the location of the file with the :ref:`bandpass transmission function <sec:filter_response>`. 
When setting this to ``none``, BeAR will use a transmission function of unity within 
the wavelength boundaries given below.

The observational data is given in at least four columns. The first two columns respresent
the wavelength boundaries over which the photometric data should be integrated. 
The third column reprents the the observational photometry data. The units of the data depend on the
chosen forward model. For example, the ``emission`` forward model expects a radiation flux
in units of :math:`\mathrm{W} \mathrm{m}^{-2}  \mathrm{\mu m}^{-1}`, while the ``transmission`` 
spectroscopy model requires the transit depth in ppm. The foruth column contains the error of 
the observational data in the same units as the previous column.

Just like for the previous observational types, an optional weight for the photometric
data point can be included in a fifth column. This is shown in the example below.

.. include:: ../examples/wasp-43b_spitzer_2.dat
   :literal:
   
Unlike the input for spectroscopic data, no instrument line profile is used here. Since photometry data
is integrated over a wider bandpass anyway, the impact of a Gaussian line profile would be 
negligible.


.. _sec:instrument_line_profile:

Instrument Line Profile
-----------------------

Optionally, before integrating the spectrum, a high-resolution spectrum can be convolved 
with a given instrument line profile to simulate the flux received by the detector. The
instrument line profile can be used for spectroscopy and band spectroscopy observations.
It is not required for photometry since in this case the flux is already integrated over
a wider filter bandpass.

.. image:: ../images/line_profile.png
  :width: 400
  :align: center

With an ideal spectrograph, the flux at a given wavelength will only be received by
a single pixel, following the dispersion relation of the instrument. 
This would correspond to the blue line in the above figure.
In the real world, due to the finite slit width of a spectrograph,
flux at a given wavelengths, however, will be spread out across several pixels. This is
depicted by the red curve in the above plot.

This spread is usually described by a Gaussian profile with a given, instrument-dependent width.
In order to compute the flux at a given pixel, BeAR needs to take into account the spread of
flux at a discrete wavelength over multiple pixels.

The Gaussian is described by its corresponding full width at half maximum (FWHM) in wavelength units. 
This, usually wavelength-dependent, FWHM needs to be supplied as an input to BeAR. In order to save 
computation time, BeAR will limit the contributions of the profile to a distance of five standard deviations 
from the profile centre.



.. _sec:filter_response:

Filter transmission function
----------------------------

BeAR has the option to use filter transmission functions to simulate the 
passing of light through a specific filter before it reaches the detector.

This is typically used for photometric observations. However, for the unlikely 
case that for a given spectroscopic observation a filter has been placed
before the spectrograph or grism, both spectroscopy and band spectroscopy also
allow for the use of a filter transmission function.

For a given wavelength-dependent filter transmission function :math:`T(\lambda)`,
the flux :math:`F(\lambda)` the integrated, photometric flux after the filter
is given by 

.. math::
  F_\mathrm{phot} = \frac{\int F(\lambda) T(\lambda) \mathrm{d} \lambda} 
    {\int T(\lambda) \mathrm{d} \lambda}

for an energy counting detector and 

.. math::
  F_\mathrm{phot} = \frac{\int F(\lambda) \lambda T(\lambda) \mathrm{d} \lambda} 
    {\int T(\lambda) \lambda \mathrm{d} \lambda}

for a photon counter. The additional factor :math:`\lambda` in the latter case
converts the energy flux :math:`F(\lambda)` into a photon flux. In case of the two
spectroscopic observational modes, the computed fluxes are simply multiplied by the
filter transmission function.

Input file structure
....................

A basic example for an filter transmission input file is shown below.

.. include:: ../examples/k_filter_transmission.dat
   :literal:

It has to contain the definition for the detector, which is either an energy counter 
(``Energy counter``) or a photon counter (``Photon counter``).

This is followed by two columns, with the wavelength in micrometers in the first
and the filter response function in the second column. The provided filter transmission
curve will be interpolated onto the internal high-resolution spectral grid used
by BeAR for a given retrieval calculation.

BeAR already comes with a set of selected filter response function files. They can
be found in the folder ``telescope_data``.