
.. _sec:temperature_profiles:

Temperature profiles
====================

BeAR currently includes the following parametrisations for atmospheric
temperature-pressure profiles

  - Milne's solution
  
  - Guillot temperature profile
  
  - Isoprofile

  - Piecewise polynomials

  - Cubic b splines
  
Details on the implementation of the latter three parametrisations are 
discussed :ref:`here <sec:profile_parametrisations>`. 


Milne's solution
----------------

Milne's solution is a simple, analytic temperature profile that has
been derived for a self-luminous object with a grey atmosphere. 
In the corresponding forward_model.config of the chosen ``forward_model.config``, 
this parametrisation is chosen by using the keyword :code:`milne`:

.. code:: 

   #Temperature profile
   milne

Milne's profile is given by

.. math:: T^4(\tau) = \frac{3}{4} T_\mathrm{eff}^4 \left(\tau + q(\tau) \right) ,

where :math:`T_\mathrm{eff}` is the effective temperature of the self-luminous object, 
:math:`\tau = \int \kappa(z) \mathrm d z` the grey optical depth for the wavelength-dependent
opacity :math:`\kappa` (usually taken to be the Rosseland mean opacity), 
and :math:`q` the Hopf function. In the Eddington approximation, :math:`q` would be 
:math:`1/3`. Internally, BeAR has a parametrised form of :math:`q(\tau)`, taken from
Mihalas (1979).

The Milne temperature profile parametrisation has two free parameters that have to
be listed in the prior config file in the following order:

  - :math:`\kappa`, the constant grey opacity in units of :math:`\mathrm{cm}^2/\mathrm{g}`

  - the effective temperature :math:`T_\mathrm{eff}`
  

Guillot profile
---------------

The analytical profiles by Tristan Guillot 
(`Guillot (2010) <https://ui.adsabs.harvard.edu/abs/2010A%26A...520A..27G/>`_) have been
derived for irridiated planetary atmospheres. BeAR contains two different implementations
of the Guillot profile, one for an incident stellar beam (Eq. 27 from Guillot, 2010) and 
for isotropic stellar irridiation (Eq. 29 from Guillot, 2010). They are chosen by setting 
the temperature profile to :code:`guillot beam`:

.. code:: 

   #Temperature profile
   guillot beam

and :code:`guillot iso`, respectively:

.. code:: 

   #Temperature profile
   guillot iso

The stellar beam case has five free parameters that have to appear ordered as shown below
in the prior config file:

  - :math:`\kappa`, the (constant) grey thermal opacity in units of :math:`\mathrm{cm}^2/\mathrm{g}`
  
  - :math:`T_\mathrm{irr}`, the temperature of incident stellar radiation (assumed to be 
    black body radiation)

  - :math:`T_\mathrm{int}`, the planet's internal temperature. For self-luminous objects,
    this would be the effective temperature

  - :math:`\gamma`, the ratio of the short-wave to the thermal opacity

  - :math:`\mu_*`, the cosine of the zenith angle of the incident stellar beam


When using isotropic stellar insolation, the following five free parameters that have to appear 
ordered as shown below in the prior config file:

  - :math:`\kappa`, the (constant) grey thermal opacity in units of :math:`\mathrm{cm}^2/\mathrm{g}`
  
  - :math:`T_\mathrm{irr}`, the temperature of incident stellar radiation (assumed to be 
    black body radiation)

  - :math:`T_\mathrm{int}`, the planet's internal temperature. For self-luminous objects,
    this would be the effective temperature

  - :math:`\gamma`, the ratio of the short-wave to the thermal opacity

  - :math:`f`, energy distribution factor. Typical values are :math:`f=1` for the sub-stellar point,
    :math:`f=1/2` for the day-side average, and :math:`f=1/4` for an average over the whole
    planetary surface.


Isoprofiles
-----------

When using an isoprofile, a constant temperature as a function of pressure is set internally.
This is the most common temperature profile employed for transmission spectra of exoplanets.
In the ``forward_model.config`` file it is chosen by using the keyword :code:`const`:

.. code:: 

   #Temperature profile
   const

The isoprofile parametrisation has a single free parameter: the constant temperature.


Piecewise polynomials
---------------------

If the temperature profile should be described by using
the parametrisation with piecewise polynomials, the keyword :code:`poly`, followed by 
the number of elements :code:`k`, and the polynomial degree :code:`q` need to be added to 
the ``forward_model.config`` file:

.. code:: 

   #Temperature profile
   poly k q

As stated :ref:`here <sec:profile_parametrisations>`, the number of free parameters for this
parametrisation is :math:`k q + 1`.
The first free parameter in the prior distrubution file refers to the bottom temperature. 
All subsequent :math:`k q` parameters are factors :math:`b_i`, such that the temperature at the
control point :math:`i` is determined by the one from the previous control point :math:`i-1` 
following :math:`T_i = T_{i-1} b_i`. The control points are distributed 
equidistantly in logarithmic pressure space.

The use of the :math:`b_i` factors allow some control over the general form of the temperature
profile. By not allowing them to exceed unity, for example, temperature inversions can be 
prevented.


Cubic b splines
---------------

For a description of the temperature profile by cubic b spline the keyword 
:code:`cubicbspline`, followed by the number of control 
points :code:`k` need to be added to the ``forward_model.config`` file:

.. code:: 

   #Temperature profile
   cubicbspline k

As stated :ref:`here <sec:profile_parametrisations>`, the number of free parameters for this
parametrisation is equal to the number of control points :math:`k` and needs to be at least 5.

The first free parameter in the prior distrubution file refers to the bottom temperature. 
All subsequent :math:`k-1` parameters are factors :math:`b_i`, such that the temperature at the
control point :math:`i` is determined by the one from the previous control point :math:`i-1` 
following :math:`T_i = T_{i-1} b_i`. The control points are distributed 
equidistantly in logarithmic pressure space.

The use of the :math:`b_i` factors allow some control over the general form of the temperature
profile. By not allowing them to exceed unity, for example, temperature inversions can be 
prevented.
