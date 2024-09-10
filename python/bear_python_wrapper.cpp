
#ifdef _SETUP_PY
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#else
#include "../_deps/pybind11-src/include/pybind11/pybind11.h"
#include "../_deps/pybind11-src/include/pybind11/stl.h"
#endif

#include "../src/config/global_config.h"
#include "../src/spectral_grid/spectral_grid.h"
#include "../src/forward_model/transmission/transmission.h"


namespace py = pybind11;

PYBIND11_MODULE(pybear, m) {
    py::class_<bear::GlobalConfig>(m, "GlobalConfig")
        .def(py::init<>())
        .def("loadConfigFile", &bear::GlobalConfig::loadConfigFile)
        .def_readwrite("forward_model_type", &bear::GlobalConfig::forward_model_type)
        .def_readwrite("wavenumber_file_path", &bear::GlobalConfig::wavenumber_file_path)
        .def_readwrite("cross_section_file_path", &bear::GlobalConfig::cross_section_file_path)
        .def_readwrite("spectral_disecretisation", &bear::GlobalConfig::spectral_disecretisation)
        .def_readwrite("const_wavenumber_step", &bear::GlobalConfig::const_wavenumber_step)
        .def_readwrite("const_wavelength_step", &bear::GlobalConfig::const_wavelength_step)
        .def_readwrite("const_spectral_resolution", &bear::GlobalConfig::const_spectral_resolution)
        .def_readwrite("use_gpu", &bear::GlobalConfig::use_gpu);

    py::class_<bear::SpectralGrid>(m, "SpectralGrid")
        .def(py::init<bear::GlobalConfig*>())
        .def(py::init<bear::GlobalConfig*, double, double>())
        .def_readwrite("wavenumber_list", &bear::SpectralGrid::wavenumber_list)
        .def_readwrite("wavelength_list", &bear::SpectralGrid::wavelength_list)
        .def("wavelengthToWavenumber", static_cast<double (bear::SpectralGrid::*)(const double)>(&bear::SpectralGrid::wavelengthToWavenumber));

    py::class_<bear::TransmissionModel>(m, "TransmissionModel")
        .def(py::init<bear::GlobalConfig*, bear::SpectralGrid*, const size_t, const std::vector<std::string>, const std::vector<std::string>>())
        .def("calcSpectrum", &bear::TransmissionModel::calcSpectrum);
}