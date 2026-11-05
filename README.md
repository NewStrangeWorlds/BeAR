# BeAR - The Bern Atmospheric Retrieval code

<p align="left">
  <picture>
    <img alt="BeAR logo" src="bear_logo/bear.png" width="100">
  </picture>
</p>

#### Authors: Daniel Kitzmann, Esther Wang ####


# Overview #

BeAR is an open-source code that performs atmospheric retrievals of exoplanet and brown dwarf spectra. It is the successor of the retrieval code HELIOS-r2 introduced in Kitzmann et al. (2020). The core of BeAR is written in C++ and CUDA; it is used exclusively through its Python module `bear`.

BeAR follows a Bayesian approach: a nested-sampling (or other) sampler explores the parameter space, generates posterior distributions, and computes the Bayesian evidence. The forward models and likelihood calculations are provided by BeAR, while the sampling is driven from Python by an external sampler such as MultiNest (via PyMultiNest), nautilus, dynesty, or UltraNest. MultiNest is compiled together with BeAR and can optionally be built with MPI support to distribute the likelihood evaluations over several GPUs.

The computationally demanding parts of the code are written in CUDA. BeAR can run on a CPU alone, but this is not recommended since runtimes are then larger by factors of 10 to 100 compared to a GPU.

Supported forward models:

- emission spectra (brown dwarfs, directly imaged planets)
- transmission spectra
- secondary eclipse spectra, including a black-body planet variant
- flat-line models for model comparison
- high-resolution phase-curve spectra for cross-correlation spectroscopy (Brogi & Line 2019 and Gibson et al. 2022 likelihoods, model filtering, CHIMERA-style re-injection, exposure blurring)

Observations can be spectroscopic, band-spectroscopic, photometric, or high-resolution time series. Atmospheric structure is described by a range of temperature-profile parametrisations (piecewise polynomials, splines, PCHIP, Guillot, Madhusudhan & Seager, Milne, adiabat plus spline), chemistry models (equilibrium chemistry via FastChem, isoprofiles, centred-log-ratio priors, isotopologue ratios, vertically varying profiles), and cloud models (grey, non-grey, power-law), all configurable in TOML files.


# User Guide #

The user guide is available at https://newstrangeworlds.github.io/BeAR/ . It describes the installation and usage of BeAR, the configuration files, and the different forward models, modules, and the Python interface.


# Requirements #

- a C++17 compiler (e.g. `g++`), NVIDIA's CUDA toolkit and drivers, CMake 3.15 or newer, OpenMP
- Python 3 with development headers (the `bear` module is built with pybind11)
- a Fortran compiler and BLAS/LAPACK for MultiNest
- optionally an MPI installation with Fortran support for the MPI version of MultiNest

All other dependencies (FastChem, MultiNest, DisORT++, AddingDoublingRT, Eigen, Boost.Math, toml++, pybind11) are downloaded and compiled automatically by CMake.


# Installation #

```bash
git clone https://github.com/newstrangeworlds/bear
cd bear
mkdir build && cd build
cmake ..
make
```

Useful CMake options are `-DSM=<arch>` to set the CUDA compute capability (e.g. `sm_86`), `-DUSE_MULTINEST_MPI=ON` to additionally build `libmultinest_mpi.so`, and `-DDEBUG_CUDA=ON` for full kernel synchronisation during debugging.

After compilation, the Python module `bear` is placed in `python/lib`, together with a small helper module that records the location of the MultiNest library in the build folder.


# Running a retrieval #

Each retrieval is configured by a set of files in a retrieval folder: `retrieval.toml`, `forward_model.toml`, `post_process.toml`, `priors.config`, and `observations.list`. The retrieval is started with a small Python driver script, e.g.

```bash
python EmissionExample/run_multinest.py
```

Example folders for the different forward models are provided in the repository (`EmissionExample`, `TransmissionExample`, `SecondaryEclipseExample`, `SecondaryEclipseExampleBB`, `FlatLineExample`, `WASP-77Ab_highres`). Configuration files from older BeAR versions can be converted to the TOML format with `tools/config_to_toml.py`.


# Python Interface #

The `bear` module exposes the retrieval setup, prior transformation, likelihood, forward-model, and post-processing functions of BeAR to Python. Besides running retrievals with any Python sampler, the forward models can be called directly to compute spectra for user-defined temperature profiles and chemical compositions. Examples are provided in the `python` folder, with scripts for each forward model. A detailed description of the interface is given in the user guide.


# Licence #

BeAR is released under the GNU General Public License v3.0 (see `LICENSE`).

If you use BeAR in your work, please cite Kitzmann et al. (2020), ApJ 890, 174 and Wang et al. (2026)
