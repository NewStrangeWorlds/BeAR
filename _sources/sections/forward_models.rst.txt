
.. _sec:forward_models:

Forward Models
==============

BeAR currently includes the following forward models:

  - :ref:`Transmission spectrum <sec:forward_model_transmission>`

  - :ref:`Emission spectrum <sec:forward_model_em>`

  - :ref:`Secondary eclipse spectrum <sec:forward_model_se>`

  - :ref:`Secondary eclipse spectra with a planetary blackbody <sec:forward_model_se_bb>`

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
the :ref:`cloud models <sec:cloud_models>` are set. After that, optional
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



.. _sec:forward_model_se:

Secondary eclipse spectrum
--------------------------

This forward model computes the wavelength-dependent secondary eclipse (or occulation) 
depth :math:`D(\lambda)` of an exoplanet atmosphere, given by

.. math::
  D(\lambda) = \frac{F_p(\lambda)}{F_*(\lambda)} \left(\frac{R_p}{R_*}\right)^2 \ ,

where :math:`F_p(\lambda)` is the outgoing flux at top of the planet's atmosphere,
:math:`F_*(\lambda)` is the stellar photospheric flux, :math:`R_p` 
is the planetary radius and :math:`R_*` the radius of the host star. 
In BeAR, :math:`D(\lambda)` has units of ppm.

In the retrieval config file ``retrieval.config`` it is selected by choosing:

.. code:: 

   #Forward model
   secondary_eclipse

The secondary eclipse spectrum forward model has two general free parameters that
have to be added to the priors configuration file in the following order:

  - logarithm of surface gravity :math:`\log g` in cgs units

  - ratio of the planet's and stellar radius :math:`\mathrm{R_p/R_*}`

The ``forward_model.config`` file for the secondary eclipse spectrum model has
the following structure:

.. include:: ../examples/forward_model_se.dat
   :literal:

The first three entries refer to the vertical discretisation of the atmosphere.
This includes the number of atmospheric layers as well as the bottom and 
top-of-atmosphere pressures in units of bars.

With the next two entries, the parametrisation for the :ref:`temperature profile <sec:temperature_profiles>`, 
the :ref:`stellar spectrum <sec:stellar_spectra>`, and the :ref:`cloud models <sec:cloud_models>` are set. 

After that, the radiative transfer scheme is chosen. There are currently two 
different options:

  - ``scm`` - the short characteristic method, available for CPU and GPU

  - ``disort`` - the discrete-ordinate solver DISORT, only available for CPU

Finally, the different :ref:`chemistry models <sec:chemistry_models>`, chemical species, and the
:ref:`opacity sources <sec:opacity_data>` are selected. It is important to note that
chemical species that are used as part of the chemistry models, should normally also have
an associated opacity source. Otherwise, the impact of that species on the 
resulting spectrum might be negligible and, therefore, its abundance will be 
likely be unconstrained.

On the other hand, it is theoretically also possible to add opacity species
without including this species in any of the chemistry models. In this case, the
abundance of this species will be zero and, thus, won't show up in the spectrum.



.. _sec:forward_model_se_bb:

Secondary eclipse spectrum with planetary blackbody
---------------------------------------------------

This forward model computes the wavelength-dependent secondary eclipse (or occulation) 
depth :math:`D(\lambda)` of an exoplanet atmosphere, given by

.. math::
  D(\lambda) = \frac{F_p(\lambda)}{F_*(\lambda)} \left(\frac{R_p}{R_*}\right)^2 \ ,

where :math:`F_p(\lambda)` is the outgoing flux at top of the planet's atmosphere,
:math:`F_*(\lambda)` is the stellar photospheric flux, :math:`R_p` 
is the planetary radius and :math:`R_*` the radius of the host star. 
In BeAR, :math:`D(\lambda)` has units of ppm.

This model is a special case of the secondary eclipse spectrum model, where the planet's flux
is assumed to be blackbody radiation. This is often used to test if the observational data warrants 
a more complex model or when only a few, usually photometric, data points are available.

In the retrieval config file ``retrieval.config`` it is selected by choosing:

.. code:: 

   #Forward model
   secondary_eclipse_bb

This forward model has two general free parameters that
have to be added to the priors configuration file in the following order:

  - the planet's effective temperature in Kelvin

  - ratio of the planet's and stellar radius :math:`\mathrm{R_p/R_*}`

The ``forward_model.config`` file for the secondary eclipse spectrum model has
the following structure:

.. include:: ../examples/forward_model_se_bb.dat
   :literal:

It contains a single option related to description of the stellar spectrum. Information
on the available options can be found in the :ref:`section <sec:stellar_spectra>` on stellar spectra.


.. _sec:forward_model_em:

Emission spectrum
-----------------

This forward model computes the emission spectrum :math:`F(\lambda)` of an exoplanet or brown dwarf atmosphere.
In BeAR, :math:`F(\lambda)` has units of :math:`\mathrm{W} \mathrm{m^{-2}} \mathrm{\mu m^{-1}}`.

In the retrieval config file ``retrieval.config`` it is selected by choosing:

.. code:: 

   #Forward model
   emission

The emission spectrum forward model has two general free parameters that
have to be added to the priors configuration file in the following order:

  - logarithm of surface gravity :math:`\log g` in cgs units

  - a scaling factor :math:`f` for the radius/distance relationship, 
    where the radius is internally set to 1 Jupiter radius

  - distance to the object in parsecs

The ``forward_model.config`` file for the emission spectrum model has
the following structure:

.. include:: ../examples/forward_model_em.dat
   :literal:

The first three entries refer to the vertical discretisation of the atmosphere.
This includes the number of atmospheric layers as well as the bottom and 
top-of-atmosphere pressures in units of bars.

With the next two entries, the parametrisation for the 
:ref:`temperature profile <sec:temperature_profiles>` and
the :ref:`cloud models <sec:cloud_models>` are set. 

After that, the radiative transfer scheme is set. There are currently two 
different options:

  - ``scm`` - the short characteristic method, available for CPU and GPU

  - ``disort`` - the discrete-ordinate solver DISORT, only available for CPU

Finally, the different :ref:`chemistry models <sec:chemistry_models>`, chemical species, and the
:ref:`opacity sources <sec:opacity_data>` are selected. It is important to note that
chemical species that are used as part of the chemistry models, should normally also have
an associated opacity source. Otherwise, the impact of that species on the 
resulting spectrum might be negligible and, therefore, its abundance will be 
likely be unconstrained.

On the other hand, it is theoretically also possible to add opacity species
without including this species in any of the chemistry models. In this case, the
abundance of this species will be zero and, thus, won't show up in the spectrum.

