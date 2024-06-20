
.. _sec:forward_models:

Forward Models
==============

BeAR currently includes the following forward models:

  - :ref:`Transmission spectrum <sec:forward_model_transmission>`

  - Emission spectrum

  - Secondary eclipse spectrum

  - Secondary eclipse spectra with a planetary blackbody

  - :ref:`Flat line <sec:forward_model_flat>`

The latter two models are usually used to check if the observational
data is more likely explained by more simpler models, such as a flat line.

Each model requires its own config file ``forward_model.config`` that are 
discussed below.

.. _sec:forward_model_flat:

Flat line
---------

This model simply fits a flat line through the observational data. Usually, this
model is used to test if the use of a more complex model is warranted to explain
the observation. The test is normally done by comparing their 
Bayesian evidences.

In the main retrieval config file ``retrieval.config`` it is selected by choosing:

.. code:: 

   #Forward model
   flat_line

The flat line model does not need a ``forward_model.config`` file since it
has no configurable parameters. In the prior distribution file, this model 
requires a single free parameter that determines the flat line. 

This parameter needs to have the same units as the observational data. For example,
in case of a transmission spectrum, this parameter refers to the transit depth in ppm, while
for an emission spectrum, it needs to have units of 
:math:`\mathrm{W} \mathrm{m^{-2}} \mathrm{\mu m^{-1}}`.


.. _sec:forward_model_transmission:

Transmission spectrum
---------------------

This forward model computes the wavelength-dependent transit depth 
:math:`D(\lambda)` of an exoplanet atmosphere, given by

.. math::
  D(\lambda) = \frac{R_p(\lambda)}{R_*} \ ,

where :math:`R_p(\lambda)` is the wavelength-dependent planetary radius
and :math:`R_*` the radius of the host star. In BeAR, :math:`D(\lambda)`
has units of ppm.