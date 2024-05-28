Observations
============

Supported observational types
-----------------------------

BeAR currently supports three different types of observations:

- spectroscopy

- band spectroscopy

- photometry

Based on the type of observation, the required format of the data files
differs slightly. In the following sections, the three basic types and
their required input formats are described. Observations have to be
provided in wavelength space. Internally, however, BeAR performs the
calculations in wavenumber space.

Typically, the computions are done with a higher spectral resolution than
the observations and then integrated to the observational wavelength structure.
The resolution is determined by the corresponding configuration parameter in
the retrieval.config file.


Spectroscopy
............

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

Optionally, before integrating the spectrum, the high-resolution spectrum can be convolved 
with a given instrument line profile to simulate the flux received by the detector.

.. image:: ../images/line_profile.png
  :width: 400
  :align: center

With an ideal spectrograph, the flux at a given wavelength will only be received by
a single pixel. This would correspond to the blue line in the above figure.
In the real world, due to the finite slit width of a spectrograph,
flux at a given wavelengths, however, will be spread out across several pixels. This is
depicted by the red line.

This spread is usually described by a Gaussian profile with a given, instrument-dependent width.
In order to compute the flux at a given wavelength, BeAR needs to integrate over the entire
profile to simulate the flux at the different wavelengths received by a pixel.

This, usually wavelength-dependent, FWHM of the instrument line profile needs to be supplied 
as an input to BeAR. In order to save computation time, BeAR will limit the integration
over the profile to a distance of five standard deviations from the profile centre.
