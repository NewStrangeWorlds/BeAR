.. _sec:radiative_transfer:

Radiative transfer schemes
==========================

The forward models that compute the thermal emission of an atmosphere, i.e. the
:ref:`emission <sec:forward_model_em>`, the :ref:`secondary eclipse <sec:forward_model_se>`,
and the :ref:`phase-curve <sec:forward_model_phase_curve>` models, solve the radiative
transfer equation along the vertical direction to obtain the outgoing flux at the top of
the atmosphere. The :ref:`transmission spectrum model <sec:forward_model_transmission>`
does not require a radiative transfer solver because it only needs the slant optical
depths through the atmosphere; it includes the Rayleigh scattering and cloud opacities
directly as extinction.

The radiative transfer scheme is chosen in the ``forward_model.toml`` file via

.. code::

   radiative_transfer = { model = "...", params = [...] }

BeAR currently offers the following schemes (``model`` keyword in brackets):

  - :ref:`Short characteristics <sec:rt_scm>` (``scm``) - pure absorption, no scattering

  - :ref:`Toon source function method <sec:rt_toon>` (``toon``) - two-stream scattering
    source function with a formal solution for the intensity

  - :ref:`Four-stream discrete ordinate adding <sec:rt_dda>` (``dda``) - analytic
    four-stream solution with multiple scattering

  - :ref:`Adding-doubling <sec:rt_adding_doubling>` (``adding_doubling``) - multi-stream
    reference solution with multiple scattering

  - :ref:`DisORT++ <sec:rt_disort>` (``disort``) - external discrete-ordinate solver (CPU only)

The following table summarises their main properties:

======================= ========================= ============== ===== ================================================
Keyword                 Scattering                Streams        GPU   ``params``
======================= ========================= ============== ===== ================================================
``scm``                 none                      2 angles       yes   ``[]``
``toon``                two-stream, multiple      1 - 8 angles   yes   ``[]`` or ``["n_angles"]``
``dda``                 four-stream, multiple     4              yes   ``[]``, ``["delta_m"]``, or ``["no_delta_m"]``
``adding_doubling``     multi-stream, multiple    2N             yes   ``["N"]`` or ``["N", "no_delta_m"]``
``disort``              multi-stream, multiple    4, 8, 16, 32   no    ``["n_streams"]``
======================= ========================= ============== ===== ================================================

All schemes return the upward flux at the top of the atmosphere, integrated over the
hemisphere, in units of :math:`\mathrm{W\,m^{-2}\,\mu m^{-1}}` on the spectral grid of
the retrieval. The computation is done independently for every spectral point and, on the
CPU, is parallelised over the spectral grid with OpenMP. On the GPU, every scheme (except
DisORT++) is implemented as a CUDA kernel that evaluates all spectral points in parallel;
the GPU versions use single-precision arithmetic.


.. _sec:rt_optical_properties:

Optical properties of a layer
-----------------------------

All schemes work on the layer optical properties that BeAR assembles from the different
opacity sources. For each layer between two atmospheric levels and each spectral point,
these are

  - the total optical depth :math:`\tau`, which is the sum of the gas absorption
    (tabulated cross-sections, continuum, and CIA), the gas Rayleigh scattering, and the
    cloud extinction,

  - the single scattering albedo :math:`\omega`, the ratio of the scattering optical depth
    (gas Rayleigh scattering plus the scattering part of the cloud extinction) to the
    total optical depth, and

  - the phase function of the scattering. The gas is assumed to scatter with the Rayleigh
    phase function, while clouds scatter with a Henyey-Greenstein phase function with the
    asymmetry parameter :math:`g` supplied by the cloud model. Each layer's phase function
    is the mixture of the two, weighted by their respective scattering optical depths.

The gas opacities are computed per atmospheric level and are converted to layer optical
depths with the trapezoidal rule. The cloud models, on the other hand, directly provide
the cloud optical depth, the single scattering albedo, and the asymmetry parameter per
layer, see the :ref:`section on cloud models <sec:cloud_models>`. Without a scattering
cloud model or a Rayleigh scattering opacity source, all layers are purely absorbing and
every scheme reduces to the same, scattering-free solution.

Which of these properties are actually used depends on the scheme. The short
characteristics method ignores scattering altogether. The Toon and DDA solvers use the
total single scattering albedo and the effective asymmetry parameter of the mixed phase
function (the Rayleigh contribution is treated as isotropic in the Henyey-Greenstein
sense, i.e. with :math:`g = 0`). The adding-doubling and DisORT++ solvers expand the
mixed phase function in Legendre moments and use the Rayleigh moments
(:math:`\chi_0 = 1`, :math:`\chi_2 = 0.1`) for the gas and :math:`\chi_l = g^l` for the
clouds.


.. _sec:rt_scm:

Short characteristics
---------------------

The short characteristics method is the simplest and fastest scheme in BeAR. It is chosen
with

.. code:: toml

   radiative_transfer = { model = "scm", params = [] }

It integrates the formal solution of the radiative transfer equation along two
Gauss-Legendre angles per hemisphere, assuming that the Planck function varies linearly
with optical depth inside each layer. The intensities at the two angles are integrated
over the hemisphere with the corresponding Gauss-Legendre weights to obtain the flux.

This scheme treats the atmosphere as purely absorbing. Rayleigh scattering opacities are
ignored entirely, and clouds are treated as pure absorbers, i.e. the single scattering
albedo and the asymmetry parameter provided by the cloud models are not used and the
entire cloud optical depth acts as absorption. This is appropriate for most emission and
secondary eclipse retrievals at infrared wavelengths of objects without strongly
scattering clouds.

The scheme takes no parameters. It is available on the CPU and the GPU.


.. _sec:rt_toon:

Toon source function method
---------------------------

This scheme implements the source function technique for thermal emission with multiple
scattering by `Toon et al. (1989) <https://ui.adsabs.harvard.edu/abs/1989JGR....9416287T/>`_.
It is selected with

.. code:: toml

   radiative_transfer = { model = "toon", params = [] }

A hemispheric-mean two-stream solution of the scattering problem provides the scattering
source function in every layer. The emergent intensity is then obtained from the formal
solution of the radiative transfer equation along a set of Gauss-Legendre angles, using
this source function, and integrated over the hemisphere to the flux. In the absence of
scattering, the result is identical to the short characteristics solver with the same
angular quadrature.

The optional, single parameter is the number of Gauss-Legendre angles used for the
angular integration of the intensity, between 1 and 8. By default, two angles are used,
i.e. the same quadrature as in the short characteristics scheme:

.. code:: toml

   radiative_transfer = { model = "toon", params = ["4"] }

Clouds enter with their single scattering albedo and asymmetry parameter. The two-stream
scattering source function is an approximation and is less accurate than the
:ref:`DDA <sec:rt_dda>` or :ref:`adding-doubling <sec:rt_adding_doubling>` solvers for
optically thick, strongly scattering clouds. For layers with a single scattering albedo
very close to unity, the albedo is internally limited to :math:`1 - 10^{-6}` to keep the
two-stream system non-degenerate.

The scheme is available on the CPU and the GPU.


.. _sec:rt_dda:

Four-stream discrete ordinate adding (DDA)
------------------------------------------

The DDA scheme solves the thermal radiative transfer with multiple scattering by using
the four-stream discrete ordinate adding method of
`Zhang et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013JQSRT.122...64Z/>`_, as
implemented in the Exo-FMS code (``lw_DDA``). It is chosen with

.. code:: toml

   radiative_transfer = { model = "dda", params = [] }

Each homogeneous layer is solved analytically for its :math:`2 \times 2` reflection and
transmission matrices and its thermal source terms. The layers are then combined with the
adding method from the bottom of the atmosphere upwards. The scheme uses a double-Gauss
quadrature with two Gauss-Legendre angles per hemisphere, the same angles that are used
by the short characteristics solver. Without scattering, the result is therefore
identical to the short characteristics scheme. The gas scatters with the Rayleigh phase
function, clouds with a Henyey-Greenstein phase function.

Strongly forward-peaked phase functions of large cloud particles are by default treated
with the delta-M scaling, which truncates the forward-scattering peak and rescales the
layer optical depth, single scattering albedo, and asymmetry parameter accordingly. This
can be switched off with the optional parameter ``no_delta_m`` (``delta_m`` explicitly
enables it):

.. code:: toml

   radiative_transfer = { model = "dda", params = ["no_delta_m"] }

In tests with scattering clouds, the DDA solver agrees with the adding-doubling solver with
eight quadrature points to within about one per cent, at a computational cost that is
similar to the Toon scheme. It is, therefore, the recommended scheme for retrievals with
scattering clouds. The scheme is available on the CPU and the GPU. The GPU version uses
single precision, with the layer algebra arranged such that it stays accurate for optically
thin layers and for (nearly) conservative scattering.


.. _sec:rt_adding_doubling:

Adding-doubling
---------------

The adding-doubling scheme uses the external ``AddingDoublingRT`` library that is
downloaded and compiled together with BeAR. It is chosen with

.. code:: toml

   radiative_transfer = { model = "adding_doubling", params = ["4"] }

The first, mandatory parameter is the number of quadrature points :math:`N` per
hemisphere, which must be at least 2. This corresponds to :math:`2N` streams in total.
The phase function is expanded in :math:`2N` Legendre moments (plus one additional
moment when the delta-M scaling is used). Larger values of :math:`N` increase the
accuracy of the angular integration at the expense of computation time. For reference
calculations, :math:`N = 8` is usually sufficient, while :math:`N = 2` yields a
four-stream solution that is comparable to the :ref:`DDA scheme <sec:rt_dda>`.

The delta-M scaling of forward-peaked phase functions is enabled by default. It can be
switched off with the optional second parameter ``no_delta_m`` (``delta_m`` explicitly
enables it):

.. code:: toml

   radiative_transfer = { model = "adding_doubling", params = ["8", "no_delta_m"] }

The scheme is available on the CPU and the GPU. On the GPU, the batched solver supports
scattering clouds only for up to 16 quadrature points; a larger number of quadrature points
can only be used for purely absorbing clouds and will otherwise result in an error
message.

The adding-doubling scheme is the most accurate but also the slowest of the GPU-capable
schemes and is best suited as a reference for the faster approximate solvers.


.. _sec:rt_disort:

DisORT++
--------

The discrete-ordinate solver DisORT++ is an external library that is downloaded and
compiled together with BeAR. It is chosen with

.. code:: toml

   radiative_transfer = { model = "disort", params = ["8"] }

The single, mandatory parameter is the number of streams, which must be 4, 8, 16, or 32.

DisORT++ is only available on the CPU. If ``use_gpu`` is set to ``true`` in the
``retrieval.toml`` file, BeAR will stop with an error message when this scheme is
selected. Since all other parts of the forward model would still run on the GPU, this
scheme is mostly useful for verifying the other solvers on the CPU rather than for
production retrievals.


Choosing a scheme
-----------------

For retrievals without scattering clouds and without Rayleigh scattering opacity sources,
all schemes produce the same result and the short characteristics method should be used
because it is the fastest. The DDA scheme is the recommended default when the retrieval
includes scattering clouds (for example the ``grey`` cloud model with the ``scattering``
option or the ``mie`` cloud model, see :ref:`sec:cloud_models`) or when Rayleigh
scattering is expected to affect the emission spectrum. The Toon scheme is of comparable
speed but, being based on a two-stream approximation of the scattering source function,
is less accurate for optically thick scattering clouds. The adding-doubling scheme with a
larger number of quadrature points serves as a reference solution.

Note that the single-precision GPU implementations of the scattering solvers can deviate
slightly from their double-precision CPU counterparts, in particular for atmospheres with
many optically very thin layers. If in doubt, a CPU run with the same configuration
provides the reference.
