.. _sec:high_resolution:

============================
High-resolution retrievals
============================

In addition to the classic low-resolution forward models, BeAR can perform
retrievals on high-resolution cross-correlation spectroscopy. This mode targets
ground-based, high-resolving-power time series (e.g. CRIRES+, IGRINS, SPIRou)
in which the planetary signal is not measured as an absolute flux but is instead
extracted from the Doppler shift of thousands of individual spectral lines as
the planet moves along its orbit.

BeAR implements this following the now-standard cross-correlation-to-likelihood
mapping of `Brogi & Line (2019) <https://ui.adsabs.harvard.edu/abs/2019AJ....157..114B/abstract>`_,
with the per-pixel-uncertainty extension of
`Gibson et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.4618G/abstract>`_
and an optional CHIMERA-style signal re-injection scheme.

.. note::

   The high-resolution workflow is under active development. The interface
   described here corresponds to the current ``dev_highres`` branch.


Overview
========

High-resolution retrievals in BeAR are driven by two ingredients:

* one or more observations flagged with the ``high_resolution`` modifier in
  ``observations.list`` (see :ref:`sec:highres_observations` for the data-file
  format), and
* a forward model that provides a high-resolution spectrum. This is primarily the
  :ref:`sec:forward_model_phase_curve` model, which produces both a
  low-resolution occultation spectrum and a high-resolution planet-to-star flux
  ratio :math:`F_\mathrm{p}/F_\mathrm{s}`. The transmission model can also be
  used with high-resolution observations; in that case the transit depth is
  converted to a signal :math:`-\delta` (in ppm) that is cross-correlated with the
  data. The separate ``opacity_highres`` block is, however, only available for the
  phase-curve model.

When at least one high-resolution observation is present, BeAR builds a **dual
spectral grid**: the normal low-resolution grid (set by ``spectral_resolution``)
and a separate high-resolution grid (set by ``spectral_resolution_highres`` in
``retrieval.toml``). The high-resolution grid is used for the line-by-line model
that enters the cross-correlation likelihood, and it may use its own opacity
tables via the ``opacity_highres`` block (see :ref:`sec:opacity_highres`).

The presence of high-resolution observations also automatically appends a small
set of retrieval parameters — the *high-resolution tail* — to the forward
model's own parameters. These are described below.


The orbital-velocity model
==========================

The planetary lines are shifted, exposure by exposure, by the line-of-sight
velocity of the planet. BeAR models the radial velocity at orbital phase
:math:`\phi` as

.. math::

   v_\mathrm{rad}(\phi) = \big(K_{p,\mathrm{ref}} + K_p\big)\,
   \sin\!\big(2\pi(\phi + \Delta\phi)\big)
   + \big(V_{\mathrm{sys,ref}} + V_\mathrm{sys}\big) + v_\mathrm{bary},

where :math:`v_\mathrm{bary}` is the per-exposure barycentric velocity
taken from the data file (``#barycentric_velocities``; zero if the user has
already applied the correction in pre-processing). Note that this value is
*added* to the radial velocity, so the sign convention of the values in the
data file has to be chosen accordingly. The high-resolution wavelength grid is
automatically extended by :math:`\pm 300\,\mathrm{km\,s^{-1}}` beyond the range
of the observations to accommodate the Doppler shifts. Only the planetary flux :math:`F_\mathrm{p}` is Doppler-shifted
by :math:`v_\mathrm{rad}`; the stellar flux :math:`F_\mathrm{s}` is kept at rest,
with a per-pixel correction factor :math:`F_\mathrm{s}(\lambda\,v_\mathrm{dop})
/F_\mathrm{s}(\lambda_\mathrm{rest})` applied to undo the shift that the
:math:`F_\mathrm{p}/F_\mathrm{s}` ratio would otherwise impose on the stellar
component.

The three retrieval parameters that control this model are the high-resolution
tail priors:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Prior name
     - Meaning
   * - ``kp``
     - Radial-velocity semi-amplitude :math:`K_p`, added as an **offset** to the
       reference value ``kp_ref`` from the data file (km/s).
   * - ``vsys``
     - Systemic velocity :math:`V_\mathrm{sys}`, added as an **offset** to the
       reference value ``vsys_ref`` (km/s).
   * - ``dphi``
     - Orbital-phase offset :math:`\Delta\phi` (in units of orbital phase).

Because ``kp`` and ``vsys`` are offsets from the reference values, a
well-centred prior is a narrow interval around zero. If ``kp_ref`` and
``vsys_ref`` are absent from the data file they default to zero, recovering the
older behaviour in which ``kp`` and ``vsys`` are retrieved directly.


.. _sec:exposure_blurring:

Exposure (velocity) blurring
============================

During a finite exposure the planet's line-of-sight velocity drifts, so the
recorded spectrum is the time-average of the instantaneous, Doppler-shifted
spectrum. BeAR models this smearing as a top-hat (boxcar) in velocity of full
width

.. math::

   \Delta v(\phi) = \big(K_{p,\mathrm{ref}} + K_p\big)\,
   \cos\!\big(2\pi(\phi + \Delta\phi)\big)\,
   \frac{2\pi}{P}\,t_\mathrm{exp} ,

which is maximal at conjunction (transit/eclipse) and vanishes at quadrature.
Here :math:`P` is the orbital period and :math:`t_\mathrm{exp}` the exposure
integration time. The running mean of the model over this velocity interval is
evaluated from a precomputed cumulative trapezoidal integral of the model, so
the cost is :math:`O(1)` per pixel and independent of the box width.

Exposure blurring is enabled automatically when **both** ``#orbital_period``
(in days) and ``#exposure_times`` (per-exposure, in seconds) are present in the
high-resolution data file (see :ref:`sec:highres_observations`). If either is
absent, no blurring is applied. If the number of values in ``#exposure_times``
does not match ``#nb_exposures``, a warning is printed and blurring is disabled.
Exposures with a negligible velocity drift (:math:`|\Delta v| < 10^{-3}`
km/s) are not blurred.


.. _sec:highres_filtering:

Model filtering and re-injection
================================

Ground-based high-resolution spectra are dominated by telluric and stellar
features that vary slowly in time. These are typically removed with a
data-detrending step (e.g. SysRem or an SVD/PCA decomposition), which projects
out the strongest time-correlated modes. To compare model and data on an equal
footing, the same operator must act on the model.

BeAR loads a per-order basis :math:`U` of the dominant temporal modes
(``#filtering_basis`` / ``#nb_basis_vectors``), augments it with a constant
column (mean removal), and constructs the projection matrix

.. math::

   P = U \left( U^\mathrm{T} \Lambda^2 U \right)^{-1} U^\mathrm{T} \Lambda^2 ,

where :math:`\Lambda` is a diagonal matrix with the inverse of the per-exposure
mean uncertainties (``#uncertainties``; the identity matrix if absent). The
operator :math:`(I - P)` removes the temporal modes. Two behaviours are
available:

* **Filtered model** (``#filter_model 1``, the default): the :math:`(I - P)`
  operator is applied to both the data (once, when the file is loaded) and the
  model (for every likelihood evaluation). Afterwards, both the filtered data and
  the filtered model are additionally detrended per exposure by removing a
  quadratic polynomial in pixel index. The likelihood is evaluated per order and
  exposure. This is the standard Gibson et al. (2022) fast model-filtering
  approach.
* **Data-only filtering** (``#filter_model 0``): :math:`(I - P)` is applied only
  to the data (without the subsequent detrending), while the model is centred in
  time by subtracting its temporal mean per pixel. The cross-correlation sums
  are accumulated over all exposures of an order before the logarithm is taken,
  i.e. one likelihood term per order. Since the model does not undergo the
  distortion that the filter imposes on the data, this mode should be regarded
  as a diagnostic rather than an unbiased likelihood.

.. note::

   On the GPU, the filtered mode currently supports at most 512 exposures per
   observation. Larger data sets have to be split into several observation
   files.

**CHIMERA-style re-injection** (``#reinject_model 1``) embeds the model planet
signal into the detector-unit background before filtering: the model
:math:`F_\mathrm{p}/F_\mathrm{s}` is multiplied by the SVD-captured background
:math:`S = P\,D`, computed from the raw (not detrended) data :math:`D` when the
file is loaded, and then passed through :math:`(I - P)`. This
reproduces the way a real planet signal is attenuated by the detrending. Since
:math:`(I-P)S = 0` by construction, only the planetary lines survive the
projection. Re-injection requires a ``#filtering_basis`` section and
``#filter_model 1``; BeAR aborts with an error otherwise. It is also strongly
recommended to retrieve the line-contrast parameter ``alpha`` (see below), since
otherwise the absolute line amplitude is not constrained; BeAR only prints a
warning in this case. If a ``#model_scale`` section is present as well, it is
ignored in favour of :math:`S = P\,D`. For re-injection, the stellar spectrum
that enters :math:`F_\mathrm{p}/F_\mathrm{s}` should be smooth, which is achieved
with the ``highres_stellar_smooth_sigma`` option of the phase-curve model (see
:ref:`sec:stellar_spectra`).

Finally, a dayside **phase function**

.. math::

   w(\phi) = \left[ 0.5\,\big(1 + \cos(2\pi\,\phi - \pi)\big) \right]^2 ,

following `Herman et al. (2022) <https://ui.adsabs.harvard.edu/abs/2022AJ....163..248H/abstract>`_
and Pelletier et al. (2025), can be applied per exposure to the
model by setting ``#phase_function 1`` in the data file. Here, :math:`\phi = 0`
corresponds to the transit, such that :math:`w = 1` at secondary eclipse,
:math:`1/4` at quadrature, and :math:`0` at transit. This weights the dayside
emission contribution with orbital phase, as appropriate for a phase-curve
observation. The weight is applied in all filtering modes, before the
:math:`(I - P)` projection where applicable.


.. _sec:highres_likelihood:

Likelihood
==========

BeAR maps the cross-correlation between model and data to a likelihood using one
of three modes. The mode is **selected automatically** from whether the line
contrast parameter ``alpha`` is a free parameter and whether the data file
provides per-pixel flux uncertainties:

.. list-table::
   :header-rows: 1
   :widths: 22 20 20 38

   * - Mode
     - ``alpha`` free?
     - Per-pixel :math:`\sigma`?
     - Description
   * - ``marginalized_alpha``
     - no
     - not used
     - Brogi & Line (2019); both the line-contrast scaling and the noise scale
       are analytically marginalized. Default when no ``alpha`` prior is given.
   * - ``free_alpha``
     - yes
     - no
     - Brogi & Line (2019) with an explicit line-contrast factor ``alpha``; the
       noise scale is marginalized.
   * - ``gibson``
     - yes
     - yes
     - Gibson et al. (2022), Eq. 4; per-pixel uncertainties with the noise
       scaling :math:`\beta` marginalized.

For the Brogi & Line modes, the per-exposure log-likelihood over :math:`N`
pixels is

.. math::

   \log L = -\frac{N}{2}\,
   \log\!\left(s_f^2 + \alpha^2 s_g^2 - 2\,\alpha\,R\right),

where :math:`s_f^2` and :math:`s_g^2` are the data and model variances and
:math:`R` is their cross-covariance. Internally, BeAR evaluates the equivalent,
normalised form :math:`\log L = -\tfrac{N}{2}\log\!\big(1 + (\alpha^2 s_g^2 -
2\alpha R)/s_f^2\big)`, which differs from the expression above only by a
constant and yields :math:`\log L = 0` for a model that is uncorrelated with the
data. In the ``marginalized_alpha`` mode this
reduces to :math:`\log L = -\tfrac{N}{2}\log\!\big(1 - R^2/(s_f^2 s_g^2)\big)`,
i.e. a monotonic function of the cross-correlation coefficient. In the
``free_alpha`` mode :math:`\alpha` is retrieved and directly measures the
absolute strength (line contrast) of the planetary lines relative to the model.

The ``gibson`` mode instead evaluates a :math:`\chi^2` built from the per-pixel
uncertainties,

.. math::

   \log L = -\frac{N}{2}\,\log\!\left(\frac{\chi^2}{\chi^2_\mathrm{data}}\right),
   \qquad
   \chi^2 = \chi^2_\mathrm{data}
   + \alpha^2\!\left(S_{mm} - \tfrac{S_m^2}{S_1}\right)
   - 2\,\alpha\!\left(S_{fm} - \tfrac{S_f S_m}{S_1}\right),

with the weighted sums :math:`S_1 = \sum_i 1/\sigma_i^2`,
:math:`S_f = \sum_i f_i/\sigma_i^2`, and analogously for the model, precomputed
once per observation. The sums are computed from the filtered fluxes but with the
original per-pixel uncertainties. Note that the ``gibson`` mode is only selected
if ``alpha`` is a free parameter; ``#flux_uncertainties`` present in a data file
are ignored otherwise.

.. _sec:highres_alpha:

The ``alpha`` parameter
-----------------------

The line-contrast parameter ``alpha`` becomes a free retrieval parameter **iff**
a prior named ``alpha`` is present in ``priors.config``. This single switch
determines the likelihood mode (see the table above) and, together with
``filter_model``, enables re-injection. Retrieving ``alpha`` allows the absolute
line contrast to be measured rather than assumed, at the cost of one extra
parameter; leaving it out selects the fully marginalized Brogi & Line
likelihood, which is scale-free.


Modules
=======

*Modules* are optional forward-model components that post-process the computed
spectrum. They are configured with a ``modules`` array in ``forward_model.toml``
(available for the transmission and phase-curve models), e.g.

.. code:: toml

   modules = [
     { model = "velocity_broadening", params = [] },
   ]

The ``velocity_broadening`` and ``phase_resolved_broadening`` modules take no
configuration ``params`` of their own, while ``stellar_contamination`` requires
the description of a stellar spectrum model in its ``params`` (see below). Each
module exposes a set of free retrieval parameters (name-keyed priors, listed
below).

Velocity broadening
--------------------

``velocity_broadening`` convolves the spectrum with a combined rotational and
instrumental kernel. Its four priors are:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Prior name
     - Meaning
   * - ``vsini``
     - Projected rotational velocity :math:`v \sin i` (km/s).
   * - ``sigma_inst``
     - Instrumental Gaussian broadening width (km/s).
   * - ``epsilon``
     - Linear limb-darkening coefficient of the rotational kernel.
   * - ``v_wind``
     - Additional (super-rotation) wind broadening velocity (km/s), added in
       quadrature to the instrumental Gaussian width.

The rotational kernel follows Gray (2005) with a linear limb-darkening law. The
convolution is carried out in velocity space on the (non-uniform) wavenumber grid,
i.e. the kernel is evaluated with the actual velocity step of each grid point.
Kernels narrower than the local grid spacing are skipped; the Gaussian is
truncated at :math:`5\sigma`.

Phase-resolved broadening
-------------------------

``phase_resolved_broadening`` applies the limb-darkened transit kernel of
`Brogi et al. (2016) <https://ui.adsabs.harvard.edu/abs/2016ApJ...817..106B/abstract>`_
(their Sect. 4.2), which accounts for planetary rotation and super-rotation
winds (within :math:`\pm 25^\circ` latitude) across the annulus of the planet
that is probed during transit, including the quadratic limb darkening of the
star. Despite its name, the kernel is currently static: it is constructed once
from the transit geometry and applied to the spectrum as a whole, not separately
for each exposure. Its seven priors are:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Prior name
     - Meaning
   * - ``v_eq``
     - Equatorial rotation velocity (km/s).
   * - ``v_wind``
     - Equatorial super-rotation wind velocity (km/s).
   * - ``sigma_inst``
     - Instrumental Gaussian broadening width (km/s).
   * - ``u1``
     - First quadratic stellar limb-darkening coefficient.
   * - ``u2``
     - Second quadratic stellar limb-darkening coefficient.
   * - ``rp_rs``
     - Planet-to-star radius ratio.
   * - ``impact_b``
     - Transit impact parameter (in stellar radii).

Stellar contamination
----------------------

``stellar_contamination`` models the transit light source effect of unocculted
spots and faculae. Its ``params`` array has to contain the stellar spectrum model
and its parameters, e.g. ``params = ["grid", "/data/phoenix_grid/"]`` (see
:ref:`sec:stellar_spectra`). In addition to the underlying stellar-model parameters
(prefixed ``tls_``), it exposes four priors:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Prior name
     - Meaning
   * - ``tls_delta_t_faculae``
     - Faculae temperature offset, added to the stellar effective temperature.
   * - ``tls_delta_t_spot``
     - Spot temperature offset, subtracted from the stellar effective
       temperature.
   * - ``tls_fraction_faculae``
     - Faculae covering fraction.
   * - ``tls_fraction_spot``
     - Spot covering fraction.


Prior summary
=============

For a high-resolution retrieval the ``priors.config`` file must contain, in
addition to the forward model's own parameters (see
:ref:`sec:forward_model_phase_curve`):

* the high-resolution tail ``kp``, ``vsys``, ``dphi`` (always), and
* optionally ``alpha`` (only if you want to retrieve the line contrast), plus
* any priors required by the modules listed in ``forward_model.toml``.

As everywhere in BeAR, these are matched by name and may be listed in any order
(see :ref:`sec:prior_distributions`).
