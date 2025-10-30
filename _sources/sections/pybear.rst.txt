
.. _sec:pybear:

##################
The pyBeAR package
##################

The Python interface to BeAR is called ``pyBeAR``. It allows the user to perform retrieval calculations
from within Python or call forward models directly to calculate synthetic spectra for a given temperature
structure and chemical composition. The interface is built on top of the ``PyBind11``. The installation of 
the pyBeAR package is described in the :ref:`installation section <sec:install_config>`. The compiled Python
package of pyBeAR is located in the ``python/lib`` directory of the BeAR root folder.

The ``python`` directories also contains subdirectories with the examples that show how to use pyBeAR with
the different forward models available in BeAR. The provided Python examples are designed to be run from the 
root directory of BeAR, rather than from their own directory. Thus, one of transmission spectroscopy examples
is run via:

.. code:: bash

   python python/transmission/pyBear_retrieval_0.py

pyBeAR is able to use the standard configuration files of BeAR
to set up the retrieval calculations. In addition, the entire configuration of the forward models and the 
retrieval setup can be done directly in Python. 

While the retrievals under Python can use the internal MultiNest sampler and the likelihood calculations of
BeAR, the user can also write their own corresponding functions to do this and use alternative samplers that 
are written in Python.

Below, a simple example of how to start a retrieval calculation with pyBeAR is given. A detailed
description of the different classes and functions of pyBeAR is provided :ref:`here <sec:pybear_details>`. The provided Python examples
can be used as a reference for the usage of the different classes and functions.

**********************************
Starting a simple pyBeAR retrieval
**********************************

| In the most simple case, a retrieval calculation can be started by loading the configuration files of BeAR directly. 
  Examples of this can be found in the provided Python scripts name ``pyBear_retrieval_0.py`` that are available for all forward models.

First, the pyBeAR package has to be imported:

.. code:: python

   from lib import pybear

The following command will then load the retrieval.config file from the folder ``retrieval_folder``:

.. code:: python

   model_config = pybear.Config(retrieval_folder)

The structure of this file is described in the :ref:`configuration section <sec:config_files>`.

With the retrieval configuration loaded, the retrieval object of pyBeAR can now be created:

.. code:: python

   model = pybear.Retrieval(model_config)

During the setup, the object class will read in the other configuration files for the chosen
forward model, the observational data, and the prior information. If successful, the retrieval object
is ready to be used for the retrieval calculations. The retrieval can be started by calling the
``run()`` method of the retrieval object:

.. code:: python

   model.run()

This will use the internal functions of BeAR and the MultiNest library to perform the retrieval calculations.
The MultiNest output files will be put into the folder ``retrieval_folder``.

After a retrieval calculation is finished, the posterior samples can be used to perform the post process, including
the calculation of spectra or effective temperatures. In pyBear, this can be done with the PostProcess object class.
In the most simplest case, the object is created with the already existing model configuration:

.. code:: python

   post_process = pybear.PostProcess(model_config)

During its setup, the PostProcess object will read in the postprocess configuration file in the ``retrieval_folder``.
The structure of this file depends on the chosen forward model is discussed in the :ref:`section <sec:forward_models>`
on forward models. Once the object has been created, the post process can be started by calling the ``run()`` method:

.. code:: python

   post_process.run()


**************************************
Advanced pyBeAR retrieval calculations
**************************************

Additional examples of how to set up more advanced retrieval calculations with pyBeAR are provided in the Python scripts
named ``pyBear_retrieval_1.py``, ``pyBear_retrieval_2.py``, and ``pyBear_retrieval_3.py`` in the respective forward model directories. These examples
show how to set up the forward models, the observational data, and the prior information directly in Python without
the need to use configuration files. This allows the user to have more flexibility in setting up the retrieval calculations.
The provided examples also show how to write custom likelihood functions and use alternative samplers that are written in Python.


*******************************
Calling forward models directly
*******************************

Besides running retrievals, pyBeAR can also be used to call the forward models directly from Python. This allows the user to calculate synthetic spectra
for a given set of temperature profiles and chemical compositions without the need to set up a full retrieval. Most forward models have their own example scripts located
in the respective forward model directories.


Transmission Spectrum Forward Model
###################################

This forward model calculates transmission spectra. An exampe of how to use this forward model
is provided in the script ``python/transmission/pyBear_transmission_spectrum.py``. The source code for the Phyton interface of this forward model is located in 
``python/transmission/src`` folder.


Secondary Eclipse Spectrum Forward Model
########################################

This forward model calculates secondary eclipse spectra. An exampe of how to use this forward model
is provided in the script ``python/secondary_eclipse/pyBear_secondary_eclipse_spectrum.py``. The source code for the Phyton interface of this forward model is located in 
``python/secondary_eclipse/src`` folder.


Emission Spectrum Forward Model
###############################

This forward model calculates emission spectra. An exampe of how to use this forward model
is provided in the script ``python/emission/pyBear_emission_spectrum.py``. The source code for the Phyton interface of this forward model is located in 
``python/emission/src`` folder.
