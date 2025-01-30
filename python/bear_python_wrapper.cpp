
#ifdef _SETUP_PY
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#else
#include "../_deps/pybind11-src/include/pybind11/pybind11.h"
#include "../_deps/pybind11-src/include/pybind11/stl.h"
#endif

#include "../src/config/global_config.h"
#include "../src/spectral_grid/spectral_grid.h"
#include "../src/retrieval/retrieval.h"
#include "../src/retrieval/post_process.h"
#include "../src/forward_model/transmission/transmission.h"
#include "../src/forward_model/secondary_eclipse/secondary_eclipse.h"


namespace py = pybind11;

PYBIND11_MODULE(pybear, m) {
    py::class_<bear::GlobalConfig>(m, "GlobalConfig")
        .def(py::init<>())
        .def("loadConfigFile", &bear::GlobalConfig::loadConfigFile)
        .def_readwrite("forward_model_type", &bear::GlobalConfig::forward_model_type)
        .def_readwrite("retrieval_folder_path", &bear::GlobalConfig::retrieval_folder_path)
        .def_readwrite("wavenumber_file_path", &bear::GlobalConfig::wavenumber_file_path)
        .def_readwrite("cross_section_file_path", &bear::GlobalConfig::cross_section_file_path)
        .def_readwrite("spectral_disecretisation", &bear::GlobalConfig::spectral_disecretisation)
        .def_readwrite("const_wavenumber_step", &bear::GlobalConfig::const_wavenumber_step)
        .def_readwrite("const_wavelength_step", &bear::GlobalConfig::const_wavelength_step)
        .def_readwrite("const_spectral_resolution", &bear::GlobalConfig::const_spectral_resolution)
        .def_readwrite("multinest_print_iter_values", &bear::GlobalConfig::multinest_print_iter_values)
        .def_readwrite("use_gpu", &bear::GlobalConfig::use_gpu);

    py::class_<bear::Retrieval>(m, "Retrieval")
        .def(py::init<bear::GlobalConfig*>())
        .def("nbParameters", &bear::Retrieval::nbParameters)
        .def("convertCubeParameters", &bear::Retrieval::convertCubeParameters)
        .def("computeLikelihood", &bear::Retrieval::computeLikelihood);

    py::class_<bear::PostProcess>(m, "PostProcess")
        .def(py::init<bear::GlobalConfig*>())
        .def("run", &bear::PostProcess::run);

    py::class_<bear::SpectralGrid>(m, "SpectralGrid")
        .def(py::init<bear::GlobalConfig*>())
        .def(py::init<bear::GlobalConfig*, double, double>())
        .def_readwrite("wavenumber_list", &bear::SpectralGrid::wavenumber_list)
        .def_readwrite("wavelength_list", &bear::SpectralGrid::wavelength_list)
        .def("wavelengthToWavenumber", static_cast<double (bear::SpectralGrid::*)(const double)>(&bear::SpectralGrid::wavelengthToWavenumber));

    py::class_<bear::TransmissionModel>(m, "TransmissionModel")
        .def(py::init<
          bear::GlobalConfig*, 
          bear::SpectralGrid*, 
          const size_t, 
          const std::vector<std::string>, 
          const std::vector<std::string>>())
        .def("calcSpectrum", &bear::TransmissionModel::calcSpectrum);

    py::class_<bear::SecondaryEclipseModel>(m, "SecondaryEclipseModel")
        .def(py::init<
          bear::GlobalConfig*, 
          bear::SpectralGrid*, 
          const size_t, 
          const std::vector<double>, 
          const std::vector<double>, 
          const std::vector<std::string>, 
          const std::vector<std::string>>())
        .def("calcSpectrum", &bear::SecondaryEclipseModel::calcSpectrum);
}