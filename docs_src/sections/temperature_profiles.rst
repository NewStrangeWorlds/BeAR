
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

  - PCHIP (monotone piecewise cubic Hermite interpolation)

  - Madhusudhan & Seager (2009) analytic profile

  - Adiabatic spline (radiative--convective profile with an adiabatic lower part)

Details on the implementation of the piecewise polynomial, cubic b spline, and
PCHIP parametrisations are discussed :ref:`here <sec:profile_parametrisations>`.


Milne's solution
----------------

Milne's solution is a simple, analytic temperature profile that has
been derived for a self-luminous object with a grey atmosphere.
In the ``forward_model.toml`` file, this parametrisation is chosen by setting
the ``temperature`` model to :code:`milne`:

.. code:: toml

   temperature = { model = "milne", params = [] }

Milne's profile is given by

.. math:: T^4(\tau) = \frac{3}{4} T_\mathrm{eff}^4 \left(\tau + q(\tau) \right) ,

where :math:`T_\mathrm{eff}` is the effective temperature of the self-luminous object, 
:math:`\tau = \int \kappa(z) \mathrm d z` the grey optical depth for the wavelength-dependent
opacity :math:`\kappa` (usually taken to be the Rosseland mean opacity), 
and :math:`q` the Hopf function. In the Eddington approximation, :math:`q` would be 
:math:`1/3`. Internally, BeAR has a parametrised form of :math:`q(\tau)`, taken from
Mihalas (1979).

The Milne temperature profile parametrisation has two free parameters that have to
be listed in the prior config file (identified by their prior names):

  - :code:`temp_kappa_ross`: :math:`\kappa`, the constant grey opacity in units of
    :math:`\mathrm{cm}^2/\mathrm{g}`

  - :code:`temp_eff`: the effective temperature :math:`T_\mathrm{eff}`
  

Guillot profile
---------------

The analytical profiles by Tristan Guillot 
(`Guillot (2010) <https://ui.adsabs.harvard.edu/abs/2010A%26A...520A..27G/>`_) have been
derived for irradiated planetary atmospheres. BeAR contains two different implementations
of the Guillot profile, one for an incident stellar beam (Eq. 27 from Guillot, 2010) and 
for isotropic stellar irradiation (Eq. 29 from Guillot, 2010). They are chosen by setting
the ``params`` of the ``guillot`` model to :code:`beam`:

.. code:: toml

   temperature = { model = "guillot", params = ["beam"] }

and :code:`isotropic`, respectively:

.. code:: toml

   temperature = { model = "guillot", params = ["isotropic"] }

The stellar beam case has five free parameters, identified in the prior config file
by their prior names:

  - :code:`temp_kappa_ir`: :math:`\kappa`, the (constant) grey thermal opacity in units of
    :math:`\mathrm{cm}^2/\mathrm{g}`

  - :code:`temp_irr`: :math:`T_\mathrm{irr}`, the temperature of incident stellar radiation
    (assumed to be black body radiation)

  - :code:`temp_int`: :math:`T_\mathrm{int}`, the planet's internal temperature. For
    self-luminous objects, this would be the effective temperature

  - :code:`temp_gamma`: :math:`\gamma`, the ratio of the short-wave to the thermal opacity

  - :code:`temp_mu`: :math:`\mu_*`, the cosine of the zenith angle of the incident stellar beam


When using isotropic stellar insolation, the five free parameters, identified in the
prior config file by their prior names, are:

  - :code:`temp_kappa_ir`: :math:`\kappa`, the (constant) grey thermal opacity in units of
    :math:`\mathrm{cm}^2/\mathrm{g}`

  - :code:`temp_irr`: :math:`T_\mathrm{irr}`, the temperature of incident stellar radiation
    (assumed to be black body radiation)

  - :code:`temp_int`: :math:`T_\mathrm{int}`, the planet's internal temperature. For
    self-luminous objects, this would be the effective temperature

  - :code:`temp_gamma`: :math:`\gamma`, the ratio of the short-wave to the thermal opacity

  - :code:`temp_f`: :math:`f`, energy distribution factor. Typical values are :math:`f=1`
    for the sub-stellar point, :math:`f=1/2` for the day-side average, and :math:`f=1/4`
    for an average over the whole planetary surface.


Isoprofiles
-----------

When using an isoprofile, a constant temperature as a function of pressure is set internally.
This is the most common temperature profile employed for transmission spectra of exoplanets.
In the ``forward_model.toml`` file it is chosen by using the model keyword :code:`const`:

.. code:: toml

   temperature = { model = "const", params = [] }

The isoprofile parametrisation has a single free parameter, the constant temperature,
identified in the prior config file by the prior name :code:`temperature`.


Piecewise polynomials
---------------------

If the temperature profile should be described by using
the parametrisation with piecewise polynomials, the ``poly`` model is used, with the number
of elements :code:`k` and the polynomial degree :code:`q` supplied as ``params`` in the
``forward_model.toml`` file:

.. code:: toml

   temperature = { model = "poly", params = ["k", "q"] }

As stated :ref:`here <sec:profile_parametrisations>`, the number of free parameters for this
parametrisation is :math:`k q + 1`.

By default, this model uses the *relative* parametrisation. The first free parameter,
identified by the prior name :code:`temp_t0`, refers to the bottom temperature. All
subsequent :math:`k q` parameters, named :code:`temp_b1`, :code:`temp_b2`, ..., are factors
:math:`b_i`, such that the temperature at the control point :math:`i` is determined by the
one from the previous control point :math:`i-1` following :math:`T_i = T_{i-1} b_i`. The
control points are distributed equidistantly in logarithmic pressure space.

The use of the :math:`b_i` factors allow some control over the general form of the temperature
profile. By not allowing them to exceed unity, for example, temperature inversions can be
prevented.

Alternatively, the *absolute* parametrisation can be selected by appending the keyword
:code:`absolute` to the ``params`` array, e.g. ``params = ["k", "q", "absolute"]``. In this
mode every free parameter is the temperature at its control point directly, using the prior
names :code:`temp_t0`, :code:`temp_t1`, ..., :code:`temp_t(kq)`. See
:ref:`here <sec:profile_parametrisations>` for a description of both modes.


Cubic b splines
---------------

For a description of the temperature profile by cubic b spline the ``cubicbspline`` model
is used, with the number of control points :code:`k` supplied as ``params`` in the
``forward_model.toml`` file:

.. code:: toml

   temperature = { model = "cubicbspline", params = ["k"] }

As stated :ref:`here <sec:profile_parametrisations>`, the number of free parameters for this
parametrisation is equal to the number of control points :math:`k` and needs to be at least 5.

By default, this model uses the *relative* parametrisation. The first free parameter,
identified by the prior name :code:`temp_t0`, refers to the bottom temperature. All
subsequent :math:`k-1` parameters, named :code:`temp_b1`, :code:`temp_b2`, ..., are factors
:math:`b_i`, such that the temperature at the control point :math:`i` is determined by the
one from the previous control point :math:`i-1` following :math:`T_i = T_{i-1} b_i`. The
control points are distributed equidistantly in logarithmic pressure space.

The use of the :math:`b_i` factors allow some control over the general form of the temperature
profile. By not allowing them to exceed unity, for example, temperature inversions can be
prevented.

Alternatively, the *absolute* parametrisation can be selected by appending the keyword
:code:`absolute` to the ``params`` array, e.g. ``params = ["k", "absolute"]``. In this mode
every free parameter is the temperature at its control point directly, using the prior names
:code:`temp_t0`, :code:`temp_t1`, ..., :code:`temp_t(k-1)`. See
:ref:`here <sec:profile_parametrisations>` for a description of both modes.


.. _sec:pchip_profile:

PCHIP profile
-------------

The ``pchip`` model describes the temperature profile with a *monotone piecewise
cubic Hermite interpolating polynomial* (PCHIP). Like the cubic b spline, the
profile is anchored at :math:`k` control points that are distributed equidistantly
in logarithmic pressure space between the bottom and the top of the atmosphere.
Between two neighbouring control points the temperature is represented by a cubic
Hermite polynomial whose derivatives at the knots are chosen following the
Fritsch--Carlson construction. This guarantees that the interpolant is *shape
preserving*: it never overshoots the control-point values and stays monotone on
every interval on which the control points are monotone. Compared to the cubic
b spline, this suppresses the spurious oscillations ("ringing") that global
splines can introduce between widely spaced control points. BeAR uses the
``pchip`` implementation from the Boost.Math interpolators library.

In the ``forward_model.toml`` file this parametrisation is chosen with the
``pchip`` model, with the number of control points :code:`k` supplied as
``params``:

.. code:: toml

   temperature = { model = "pchip", params = ["k"] }

The number of free parameters is equal to the number of control points :math:`k`,
which must be at least 4. Temperatures below 50 K are clipped to 50 K and cause the
model to be flagged as neglected.

Unlike the piecewise-polynomial and cubic-b-spline profiles, the PCHIP profile
uses the *absolute* parametrisation **by default**. In this mode every free
parameter is the temperature at its control point directly, identified by the
prior names :code:`temp_t0`, :code:`temp_t1`, ..., :code:`temp_t(k-1)`, where
:code:`temp_t0` is the temperature at the bottom (deepest, highest-pressure)
control point.

The *relative* parametrisation can be selected by appending the keyword
:code:`relative` to the ``params`` array, e.g. ``params = ["k", "relative"]``. In
that case the first free parameter :code:`temp_t0` is again the bottom temperature,
while the remaining :math:`k-1` parameters, named :code:`temp_b1`, :code:`temp_b2`,
..., are multiplicative factors :math:`b_i` with :math:`T_i = T_{i-1} b_i`. See
:ref:`here <sec:relative_absolute>` for a description of both modes.


Madhusudhan & Seager profile
----------------------------

This is the analytic three-layer pressure--temperature parametrisation introduced
by `Madhusudhan & Seager (2009)
<https://ui.adsabs.harvard.edu/abs/2009ApJ...707...24M/>`_. It is widely used in
exoplanet atmosphere retrievals because it can represent both monotone profiles
and thermal inversions with a small number of physically motivated parameters. In
the ``forward_model.toml`` file it is chosen with the ``madhusudhan_seager`` model
and takes no additional ``params``:

.. code:: toml

   temperature = { model = "madhusudhan_seager", params = [] }

The atmosphere is divided into three layers, separated by the pressure levels
:math:`P_1` and :math:`P_3`. With :math:`P_0` the pressure at the top of the
atmosphere and :math:`T_0` the temperature there, the profile is given by

.. math::

   T(P) = \begin{cases}
     \left[ \dfrac{\ln(P/P_0)}{\alpha_1} \right]^2 + T_0
       & P < P_1  \quad \text{(layer 1)} \\[2ex]
     \left[ \dfrac{\ln(P/P_2)}{\alpha_2} \right]^2 + T_2
       & P_1 \le P < P_3 \quad \text{(layer 2)} \\[2ex]
     T_3 & P \ge P_3 \quad \text{(layer 3)} .
   \end{cases}

Layer 1 is the upper radiative region above :math:`P_1`, layer 2 the intermediate
region between :math:`P_1` and :math:`P_3`, and layer 3 the deep, isothermal region
below :math:`P_3`. The reference temperature :math:`T_2` of layer 2 is fixed by
enforcing continuity of the profile at :math:`P_1`,

.. math::

   T_2 = \left[ \frac{\ln(P_1/P_0)}{\alpha_1} \right]^2 + T_0
       - \left[ \frac{\ln(P_1/P_2)}{\alpha_2} \right]^2 ,

and the isothermal deep-layer temperature is the value of the layer-2 expression at
:math:`P_3`, i.e. :math:`T_3 = \left[ \ln(P_3/P_2)/\alpha_2 \right]^2 + T_2`. The
pressure :math:`P_2` sets the location of the temperature extremum of layer 2; a
thermal inversion arises when :math:`P_2` lies inside the layer, allowing the model
to describe non-monotone profiles.

The parametrisation has six free parameters, identified in the prior config file by
their prior names:

  - :code:`temp_t0`: :math:`T_0`, the temperature at the top of the atmosphere in K

  - :code:`temp_log_p1`: :math:`\log_{10}(P_1)`, the pressure of the boundary
    between layers 1 and 2 (in log\ :sub:`10` bar)

  - :code:`temp_log_p2`: :math:`\log_{10}(P_2)`, the reference pressure of layer 2
    that sets the location of the temperature extremum (in log\ :sub:`10` bar)

  - :code:`temp_log_p3`: :math:`\log_{10}(P_3)`, the pressure below which the
    atmosphere is isothermal (in log\ :sub:`10` bar)

  - :code:`temp_alpha1`: :math:`\alpha_1`, the gradient parameter of layer 1
    (:math:`> 0`); smaller values give steeper profiles

  - :code:`temp_alpha2`: :math:`\alpha_2`, the gradient parameter of layer 2
    (:math:`> 0`)


Adiabatic spline
----------------

The ``adiabate_spline`` model is a radiative--convective profile that joins an
adiabatic deep region to a spline-interpolated upper region at a
radiative--convective boundary (RCB). Below the RCB (i.e. at pressures higher than
the boundary pressure), where the atmosphere is assumed to be convective, the
temperature follows a dry adiabat

.. math:: T(P) = T_\mathrm{bot} \left( \frac{P}{P_\mathrm{bot}} \right)^{1 - 1/\gamma} ,

anchored at the bottom temperature :math:`T_\mathrm{bot}` and the bottom (highest)
pressure :math:`P_\mathrm{bot}`. Above the RCB (lower pressures), the radiative
region is described by a monotone piecewise cubic Hermite (PCHIP) spline, using the
same Boost.Math interpolator as the :ref:`PCHIP profile <sec:pchip_profile>`
above. The spline is anchored at the RCB, where its temperature is set to the
adiabat value at the boundary so that the profile is continuous, and reaches up to
the top of the atmosphere through :math:`k` control points distributed equidistantly
in logarithmic pressure between the RCB and the top of the atmosphere.

In the ``forward_model.toml`` file this parametrisation is chosen with the
``adiabate_spline`` model, with the number of spline control points :code:`k`
supplied as ``params``:

.. code:: toml

   temperature = { model = "adiabate_spline", params = ["k"] }

The number of control points :math:`k` must be at least 4. The total number of free
parameters is :math:`k + 2`: the three parameters describing the RCB and the
adiabat, followed by the :math:`k-1` control-point temperatures of the radiative
spline (the deepest spline point is fixed by continuity with the adiabat and is not
a free parameter). The prior names, in order, are:

  - :code:`temp_pressure_rcb`: the pressure of the radiative--convective boundary
    (in bar). Levels at higher pressure follow the adiabat, levels at lower pressure
    follow the spline

  - :code:`temp_conv_gamma`: :math:`\gamma`, the adiabatic index that sets the
    adiabatic gradient of the deep, convective region through the exponent
    :math:`1 - 1/\gamma`

  - :code:`temp_bottom`: :math:`T_\mathrm{bot}`, the temperature at the bottom
    (deepest, highest-pressure) level of the atmosphere in K

  - :code:`temp_t1`, :code:`temp_t2`, ..., :code:`temp_t(k-1)`: the absolute
    temperatures at the :math:`k-1` control points of the radiative spline above the
    radiative--convective boundary in K
