
.. _sec:python_interface:

####################
The Python interface
####################

BeAR is used through its Python module, simply called ``bear``. It allows the user to perform retrieval calculations
from within Python or call forward models directly to calculate synthetic spectra for a given temperature
structure and chemical composition. The interface is built on top of ``PyBind11`` and the compiled
extension module is named ``bear``. The installation of
the ``bear`` module is described in the :ref:`installation section <sec:install_config>`. The compiled Python
module (``bear``) is placed in the ``python/lib`` directory of the BeAR root folder. The MultiNest shared
library ``libmultinest.so`` stays in the CMake build directory; its location is recorded by the build in the
small helper module ``python/lib/bear_multinest_path.py`` (see below).

The simplest way to start is one of the retrieval example folders in the BeAR root directory, e.g.
``EmissionExample``, each of which contains a driver script ``run_multinest.py`` that loads the configuration
files of the example and runs the retrieval and the post process with the ``bear`` module and PyMultiNest:

.. code:: bash

   python EmissionExample/run_multinest.py

The ``python`` directory also contains subdirectories with further examples that show how to use the ``bear`` module with
the different forward models available in BeAR. The provided Python examples in these subdirectories are
designed to be run from the root directory of BeAR, rather than from their own directory. Thus, one of the
transmission spectroscopy examples is run via:

.. code:: bash

   python python/transmission/bear_retrieval_1.py

Note that the scripts named ``bear_retrieval_0.py`` in these subdirectories call a ``model.run()`` method
that is not provided by the current ``bear`` module and therefore do not work; use ``run_multinest.py`` of the
example folders or the ``bear_retrieval_1.py`` to ``bear_retrieval_3.py`` scripts instead.

The ``bear`` module is able to use the standard configuration files of BeAR
to set up the retrieval calculations. In addition, the entire configuration of the forward models and the
retrieval setup can be done directly in Python.

Retrievals under Python are driven by an external sampler. The ``bear`` module provides the prior transformation and
likelihood calculations of BeAR, which are handed to a Python sampler such as MultiNest (via PyMultiNest)
or dynesty. The user can equally write their own corresponding functions and use any other sampler that is
written in Python.

Below, a simple example of how to start a retrieval calculation is given. A detailed
description of the different classes and functions of the ``bear`` module is provided :ref:`here <sec:python_interface_details>`. The provided Python examples
can be used as a reference for the usage of the different classes and functions.

***************************
Starting a simple retrieval
***************************

| In the most simple case, a retrieval calculation can be started by loading the configuration files of BeAR directly.
  Examples of this can be found in the driver scripts named ``run_multinest.py`` provided in every retrieval example
  folder. The high-resolution example ``WASP-77Ab_highres`` additionally contains a ``run_nautilus.py`` that drives
  the same retrieval with the nautilus sampler, and ``python/transmission/bear_retrieval_1_dynesty.py`` shows
  the equivalent setup with dynesty.

First, the compiled ``bear`` module has to be imported. In the example scripts this is done via:

.. code:: python

   from lib import bear

To use the MultiNest sampler, its shared library also has to be loaded before ``pymultinest`` is imported.
The build writes a helper module ``bear_multinest_path`` that records the location of ``libmultinest.so``:

.. code:: python

   from lib.bear_multinest_path import MULTINEST_LIB_DIR
   import ctypes as _ctypes
   _ctypes.CDLL(MULTINEST_LIB_DIR + '/libmultinest.so', mode=_ctypes.RTLD_GLOBAL)
   import pymultinest

The following command will then load the ``retrieval.toml`` file from the folder ``retrieval_folder``:

.. code:: python

   model_config = bear.Config(retrieval_folder)

The structure of this file is described in the :ref:`configuration section <sec:config_files>`.

With the retrieval configuration loaded, the retrieval object of ``bear`` can now be created:

.. code:: python

   model = bear.Retrieval(model_config)

During the setup, the object class will read in the other configuration files for the chosen
forward model, the observational data, and the prior information. If successful, the retrieval object
is ready to be used for the retrieval calculations.

The ``Retrieval`` object does not run a sampler by itself. Instead, it provides the prior transformation
and the likelihood, which are handed to an external sampler. For MultiNest, the prior callback calls
``model.convertCubeParameters()`` and the likelihood callback calls ``model.computeLikelihood()``. Both are
passed to ``pymultinest.run(...)`` together with the number of free parameters returned by
``model.nbParameters()``:

.. code:: python

   pymultinest.run(
     loglike,
     priors,
     model.nbParameters(),
     outputfiles_basename=retrieval_folder)

This uses the likelihood functions of BeAR to perform the retrieval calculations, and the sampler output
files are written into the folder ``retrieval_folder``.

After a retrieval calculation is finished, the posterior samples can be used to perform the post process, including
the calculation of spectra or effective temperatures. This can be done with the PostProcess object class.
In the most simplest case, the object is created with the already existing model configuration:

.. code:: python

   post_process = bear.PostProcess(model_config)

During its setup, the PostProcess object will read in the postprocess configuration file in the ``retrieval_folder``.
The structure of this file depends on the chosen forward model is discussed in the :ref:`section <sec:forward_models>`
on forward models. Once the object has been created, the post process can be started by calling the ``run()`` method:

.. code:: python

   post_process.run()


*******************************
Advanced retrieval calculations
*******************************

Additional examples of how to set up more advanced retrieval calculations are provided in the Python scripts
named ``bear_retrieval_1.py``, ``bear_retrieval_2.py``, and ``bear_retrieval_3.py`` in the respective forward model directories. These examples
show how to set up the forward models, the observational data, and the prior information directly in Python without
the need to use configuration files. This allows the user to have more flexibility in setting up the retrieval calculations.
The provided examples also show how to write custom likelihood functions and use alternative samplers that are written in Python.


*******************************
Calling forward models directly
*******************************

Besides running retrievals, the ``bear`` module can also be used to call the forward models directly from Python. This allows the user to calculate synthetic spectra
for a given set of temperature profiles and chemical compositions without the need to set up a full retrieval. Most forward models have their own example scripts located
in the respective forward model directories.


Transmission Spectrum Forward Model
###################################

This forward model calculates transmission spectra. An exampe of how to use this forward model
is provided in the script ``python/transmission/bear_transmission_spectrum.py``. The source code for the Python interface of this forward model is located in 
``python/transmission/src`` folder.


Secondary Eclipse Spectrum Forward Model
########################################

This forward model calculates secondary eclipse spectra. An exampe of how to use this forward model
is provided in the script ``python/secondary_eclipse/bear_secondary_eclipse_spectrum.py``. The source code for the Python interface of this forward model is located in 
``python/secondary_eclipse/src`` folder.


Emission Spectrum Forward Model
###############################

This forward model calculates emission spectra. An exampe of how to use this forward model
is provided in the script ``python/emission/bear_emission_spectrum.py``. The source code for the Python interface of this forward model is located in 
``python/emission/src`` folder.
