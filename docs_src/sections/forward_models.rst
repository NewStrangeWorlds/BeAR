
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
  D(\lambda) = \left(\frac{R_p(\lambda)}{R_*}\right)^2 \ ,

where :math:`R_p(\lambda)` is the wavelength-dependent planetary radius
and :math:`R_*` the radius of the host star. In BeAR, :math:`D(\lambda)`
has units of ppm.

In the retrieval config file ``retrieval.config`` it is selected by choosing:

.. code:: 

   #Forward model
   transmission

The transmission spectrum forward model has three general free parameters that
have to be added to the priors configuration file in the following order:

  - logarithm of surface gravity :math:`\log g` in cgs units

  - planet radius in Jupiter radii

  - stellar radius in solar radii

The ``forward_model.config`` file for the transmission spectrum model has
the following structure:

.. include:: ../examples/forward_model_transmission.dat
   :literal:

The first three entries refer to the vertical discretisation of the atmosphere.
This includes the number of atmospheric layers as well as the bottom and 
top-of-atmosphere pressures in units of bars.

With the next entry, the retrieval of the mean molecular or the scale height
can be chosen. Usually, BeAR will determine the mean molecular based on
the chemical abundances of the species included in the retrieval.
Sometimes, however, the background species in an atmosphere might not be
known. For such a case, BeAR can use the mean molecular weight as a free 
parameter. 

Furthermore, for transmission spectra, the surface gravity, the mean molecular
weight, and the temperature might become degenerate in a retrieval. This is
often the case when no constraints on the surface gravity can be provided
and the dominating background specise is not known. In such a scenario, BeAR
can use the (constant) atmospheric scale height in units of km as a free 
parameter.

For the standard case, the keyword ``no`` needs to be used here. If BeAR should
use the mean molecular weight as a free parameter ``mmw`` is used, while
``sh`` is used when the scale should be used as a free parameter instead. If either
the mean molecular weight or the scale height is chosen as a free parameter,
a corresponding prior needs to be added as a fourth model parameter in the
prior configuration file.

With the next two entries, the parametrisation for the 
:ref:`temperature profile <sec:temperature_profiles>` and
the :ref:`cloud models <sec:cloud_models>` are set. After than, optional
modules can be added to the forward model.

Finally, the different :ref:`chemistry models <sec:chemistry_models>`, chemical species, and the
:ref:`opacity sources <sec:opacity_data>` are selected. It is important to note that
chemical species that are used as part of the chemistry models, should normally also have
an associated opacity source. Otherwise, the impact of that species on the 
resulting spectrum might be negligible and, therefore, its abundance will be 
likely be unconstrained.

On the other hand, it is theoretically also possible to add opacity species
without including this species in any of the chemistry models. In this case, the
abundance of this species will be zero and, thus, won't show up in the spectrum.

