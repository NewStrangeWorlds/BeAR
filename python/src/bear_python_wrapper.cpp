
#ifdef _SETUP_PY
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#else
#include "../../_deps/pybind11-src/include/pybind11/pybind11.h"
#include "../../_deps/pybind11-src/include/pybind11/stl.h"
#endif

#include "../../src/config/global_config.h"
#include "../../src/spectral_grid/spectral_grid.h"
#include "../../src/retrieval/retrieval.h"
#include "../../../src/retrieval/priors.h"
#include "../../src/retrieval/post_process.h"
#include "../../src/observations/observations.h"
#include "../../src/forward_model/forward_model.h"
#include "../../../src/forward_model/generic_config.h"
#include "../../src/forward_model/transmission/transmission.h"
#include "../../src/forward_model/secondary_eclipse/secondary_eclipse.h"
#include "../../src/forward_model/emission/emission.h"


namespace py = pybind11;

PYBIND11_MODULE(pybear, m) {
    py::class_<bear::GlobalConfig>(m, "Config")
        .def(py::init<>())
        .def(py::init<const std::string>())
        .def(py::init<
            const bool, 
            const std::string, 
            const std::string, 
            const std::string, 
            const double,
            const std::string, 
            const std::string>())
        .def("loadConfigFile", &bear::GlobalConfig::loadConfigFile)
        .def_readwrite("forward_model_type", &bear::GlobalConfig::forward_model_type)
        .def_readwrite("retrieval_folder_path", &bear::GlobalConfig::retrieval_folder_path)
        .def_readwrite("multinest_output_path", &bear::GlobalConfig::multinest_output_path)
        .def_readwrite("post_output_path", &bear::GlobalConfig::post_output_path)
        .def_readwrite("wavenumber_file_path", &bear::GlobalConfig::wavenumber_file_path)
        .def_readwrite("cross_section_file_path", &bear::GlobalConfig::cross_section_file_path)
        .def_readwrite("spectral_disecretisation", &bear::GlobalConfig::spectral_disecretisation)
        .def_readwrite("const_wavenumber_step", &bear::GlobalConfig::const_wavenumber_step)
        .def_readwrite("const_wavelength_step", &bear::GlobalConfig::const_wavelength_step)
        .def_readwrite("const_spectral_resolution", &bear::GlobalConfig::const_spectral_resolution)
        .def_readwrite("multinest_print_iter_values", &bear::GlobalConfig::multinest_print_iter_values)
        .def_readwrite("multinest_efficiency", &bear::GlobalConfig::multinest_efficiency)
        .def_readwrite("multinest_nb_living_points", &bear::GlobalConfig::multinest_nb_living_points)
        .def_readwrite("multinest_nb_iterations", &bear::GlobalConfig::multinest_nb_iterations)
        .def_readwrite("multinest_feedback", &bear::GlobalConfig::multinest_feedback)
        .def_readwrite("use_error_inflation", &bear::GlobalConfig::use_error_inflation)
        .def_readwrite("use_gpu", &bear::GlobalConfig::use_gpu)
        .def_readwrite("nb_omp_processes", &bear::GlobalConfig::nb_omp_processes);

    py::class_<bear::ObservationInput>(m, "Observation")
        .def(py::init<>())
        .def(py::init<std::string, std::string>())
        .def(py::init<
          std::string, 
          std::string, 
          const std::vector<double>& ,
          const std::vector<double>& ,
          const std::vector<double>&>())
        .def(py::init<
          std::string, 
          std::string, 
          const std::vector<std::vector<double>>& ,
          const std::vector<double>& ,
          const std::vector<double>&>())
        .def_readwrite("name", &bear::ObservationInput::name)
        .def_readwrite("type", &bear::ObservationInput::type)
        .def_readwrite("wavelengths", &bear::ObservationInput::wavelengths)
        .def_readwrite("bin_wavelength_edges", &bear::ObservationInput::bin_wavelength_edges)
        .def_readwrite("data", &bear::ObservationInput::data)
        .def_readwrite("data_error", &bear::ObservationInput::data_error)
        .def_readwrite("likelihood_weight", &bear::ObservationInput::likelihood_weight)
        .def_readwrite("instrument_profile_fwhm", &bear::ObservationInput::instrument_profile_fwhm)
        .def_readwrite("filter_response", &bear::ObservationInput::filter_response)
        .def_readwrite("filter_detector_type", &bear::ObservationInput::filter_detector_type)
        .def_readwrite("spectrum_modifier_id", &bear::ObservationInput::spectrum_modifier_id);

    py::class_<bear::PriorConfig>(m, "Prior")
        .def(py::init<
           std::string, 
           std::string, 
           std::vector<double>>())
        .def(py::init<
           std::string, 
           std::string, 
           std::vector<double>,
           std::string>())
        .def_readwrite("type", &bear::PriorConfig::type)
        .def_readwrite("description", &bear::PriorConfig::description)
        .def_readwrite("parameter", &bear::PriorConfig::parameter)
        .def_readwrite("unit", &bear::PriorConfig::unit);

    py::class_<bear::ForwardModelOutput>(m, "ForwardModelOutput")
        .def(py::init<>())
        .def_readwrite("neglect_model", &bear::ForwardModelOutput::neglect_model)
        .def_readwrite("spectrum", &bear::ForwardModelOutput::spectrum)
        .def_readwrite("spectrum_obs", &bear::ForwardModelOutput::spectrum_obs);

    py::class_<bear::AtmosphereOutput>(m, "AtmosphereOutput")
        .def(py::init<>())
        .def_readwrite("neglect_model", &bear::AtmosphereOutput::neglect_model)
        .def_readwrite("pressure", &bear::AtmosphereOutput::pressure)
        .def_readwrite("altitude", &bear::AtmosphereOutput::altitude)
        .def_readwrite("temperature", &bear::AtmosphereOutput::temperature)
        .def_readwrite("species_symbols", &bear::AtmosphereOutput::species_symbols)
        .def_readwrite("mixing_ratios", &bear::AtmosphereOutput::mixing_ratios);

    py::class_<bear::Retrieval>(m, "Retrieval")
        .def(py::init<
          bear::GlobalConfig*>())
        .def(py::init<
          bear::GlobalConfig*,
          bear::GenericConfig*,
          const std::vector<bear::ObservationInput>&,
          const std::vector<bear::PriorConfig>&>())
        .def_readwrite("spectral_grid", &bear::Retrieval::spectral_grid)
        .def("run", &bear::Retrieval::run)
        .def("nbParameters", &bear::Retrieval::nbParameters)
        .def("convertCubeParameters", &bear::Retrieval::convertCubeParameters)
        .def("convertToPhysicalParameters", &bear::Retrieval::convertToPhysicalParameters)
        .def("computeLikelihood", &bear::Retrieval::computeLikelihood)
        .def("computeModel", &bear::Retrieval::computeModel);

    py::class_<bear::PostProcess>(m, "PostProcess")
        .def(py::init<bear::GlobalConfig*>())
        .def(py::init<
          bear::GlobalConfig*,
          bear::GenericConfig*,
          bear::GenericConfig*,
          const std::vector<bear::ObservationInput>&,
          const std::vector<bear::PriorConfig>&>())
        .def_readwrite("spectral_grid", &bear::PostProcess::spectral_grid)
        .def("convertToPhysicalParameters", &bear::PostProcess::convertToPhysicalParameters)
        .def("computeModel", &bear::PostProcess::computeModel)
        .def("computeAtmosphereStructure", &bear::PostProcess::computeAtmosphereStructure)
        .def("run", static_cast<bool (bear::PostProcess::*)()>(&bear::PostProcess::run))
        .def("run", static_cast<bool (bear::PostProcess::*)(
            const std::string)>(&bear::PostProcess::run))
        .def("run", static_cast<bool (bear::PostProcess::*)(
            std::vector<std::vector<double>>&, 
            std::vector<double>&)>(&bear::PostProcess::run));

    py::class_<bear::SpectralGrid>(m, "SpectralGrid")
        .def(py::init<bear::GlobalConfig*>())
        .def(py::init<bear::GlobalConfig*, double, double>())
        .def_readwrite("wavenumber_list", &bear::SpectralGrid::wavenumber_list)
        .def_readwrite("wavelength_list", &bear::SpectralGrid::wavelength_list)
        .def("wavelengthToWavenumber", static_cast<double (bear::SpectralGrid::*)(const double)>(&bear::SpectralGrid::wavelengthToWavenumber));
    
    py::class_<bear::GenericConfig> genericConfig(m, "GenericConfig");

    py::class_<bear::TransmissionModelConfig>(m, "TransmissionModelConfig", genericConfig)
        .def(py::init<const std::string&>())
        .def(py::init<
          const int,
          const double,
          const double,
          const std::string&,
          const std::vector<std::string>& ,
          const std::vector<std::string>& ,
          const std::vector<std::vector<std::string>>&,
          const std::vector<std::string>&,
          const std::vector<std::string>&>())
        .def(py::init<
          const bool, 
          const bool, 
          const bool,
          const int,
          const double,
          const double,
          const std::string&,
          const std::vector<std::string>&,
          const std::vector<std::string>&,
          const std::vector<std::vector<std::string>>&,
          const std::vector<std::string>&,
          const std::vector<std::string>&,
          const std::vector<std::string>&,
          const std::vector<std::vector<std::string>>&,
          const std::vector<std::string>&,
          const std::vector<std::vector<std::string>>&>())
        .def_readwrite("nb_grid_points", &bear::TransmissionModelConfig::nb_grid_points)
        .def_readwrite("atmos_boundaries", &bear::TransmissionModelConfig::atmos_boundaries)
        .def_readwrite("fit_mean_molecular_weight", &bear::TransmissionModelConfig::fit_mean_molecular_weight)
        .def_readwrite("fit_scale_height", &bear::TransmissionModelConfig::fit_scale_height)
        .def_readwrite("temperature_profile_model", &bear::TransmissionModelConfig::temperature_profile_model)
        .def_readwrite("temperature_profile_parameters", &bear::TransmissionModelConfig::temperature_profile_parameters)
        .def_readwrite("chemistry_model", &bear::TransmissionModelConfig::chemistry_model)
        .def_readwrite("chemistry_parameters", &bear::TransmissionModelConfig::chemistry_parameters)
        .def_readwrite("cloud_model", &bear::TransmissionModelConfig::cloud_model)
        .def_readwrite("cloud_model_parameters", &bear::TransmissionModelConfig::cloud_model_parameters)
        .def_readwrite("modules", &bear::TransmissionModelConfig::modules)
        .def_readwrite("modules_parameters", &bear::TransmissionModelConfig::modules_parameters)
        .def_readwrite("opacity_species_symbol", &bear::TransmissionModelConfig::modules)
        .def_readwrite("opacity_species_folder", &bear::TransmissionModelConfig::opacity_species_folder);

    py::class_<bear::SecondaryEclipseConfig>(m, "SecondaryEclipseConfig", genericConfig)
        .def(py::init<
          const std::string&>())
        .def(py::init<
          const int,
          const double,
          const double,
          const std::string,
          const std::vector<std::string>&,
          const std::string,
          const std::vector<std::string>&,
          const std::vector<std::string>&,
          const std::vector<std::vector<std::string>>&,
          const std::vector<std::string>&,
          const std::vector<std::string>&,
          const std::string,
          const std::vector<std::string>&,
          const std::vector<std::string>&,
          const std::vector<std::vector<std::string>>&>())
        .def_readwrite("nb_grid_points", &bear::SecondaryEclipseConfig::nb_grid_points)
        .def_readwrite("atmos_boundaries", &bear::SecondaryEclipseConfig::atmos_boundaries)
        .def_readwrite("temperature_profile_model", &bear::SecondaryEclipseConfig::temperature_profile_model)
        .def_readwrite("temperature_profile_parameters", &bear::SecondaryEclipseConfig::temperature_profile_parameters)
        .def_readwrite("radiative_transfer_model", &bear::SecondaryEclipseConfig::radiative_transfer_model)
        .def_readwrite("radiative_transfer_parameters", &bear::SecondaryEclipseConfig::radiative_transfer_parameters)
        .def_readwrite("stellar_spectrum_model", &bear::SecondaryEclipseConfig::radiative_transfer_model)
        .def_readwrite("stellar_model_parameters", &bear::SecondaryEclipseConfig::stellar_model_parameters)
        .def_readwrite("chemistry_model", &bear::SecondaryEclipseConfig::chemistry_model)
        .def_readwrite("chemistry_parameters", &bear::SecondaryEclipseConfig::chemistry_parameters)
        .def_readwrite("cloud_model", &bear::SecondaryEclipseConfig::cloud_model)
        .def_readwrite("cloud_model_parameters", &bear::SecondaryEclipseConfig::cloud_model_parameters)
        .def_readwrite("opacity_species_symbol", &bear::SecondaryEclipseConfig::opacity_species_symbol)
        .def_readwrite("opacity_species_folder", &bear::SecondaryEclipseConfig::opacity_species_folder);

    py::class_<bear::EmissionModelConfig>(m, "EmissionModelConfig", genericConfig)
        .def(py::init<
          const std::string&>())
        .def(py::init<
          const int,
          const double,
          const double,
          const std::string,
          const std::vector<std::string>&,
          const std::string,
          const std::vector<std::string>&,
          const std::vector<std::string>&,
          const std::vector<std::vector<std::string>>&,
          const std::vector<std::string>&,
          const std::vector<std::string>&>())
        .def(py::init<
          const int,
          const double,
          const double,
          const std::string,
          const std::vector<std::string>&,
          const std::string,
          const std::vector<std::string>&,
          const std::vector<std::string>&,
          const std::vector<std::vector<std::string>>&,
          const std::vector<std::string>&,
          const std::vector<std::string>&,
          const std::vector<std::string>&,
          const std::vector<std::vector<std::string>>&>())
        .def_readwrite("nb_grid_points", &bear::EmissionModelConfig::nb_grid_points)
        .def_readwrite("atmos_boundaries", &bear::EmissionModelConfig::atmos_boundaries)
        .def_readwrite("temperature_profile_model", &bear::EmissionModelConfig::temperature_profile_model)
        .def_readwrite("temperature_profile_parameters", &bear::EmissionModelConfig::temperature_profile_parameters)
        .def_readwrite("radiative_transfer_model", &bear::EmissionModelConfig::radiative_transfer_model)
        .def_readwrite("radiative_transfer_parameters", &bear::EmissionModelConfig::radiative_transfer_parameters)
        .def_readwrite("chemistry_model", &bear::EmissionModelConfig::chemistry_model)
        .def_readwrite("chemistry_parameters", &bear::EmissionModelConfig::chemistry_parameters)
        .def_readwrite("cloud_model", &bear::EmissionModelConfig::cloud_model)
        .def_readwrite("cloud_model_parameters", &bear::EmissionModelConfig::cloud_model_parameters)
        .def_readwrite("opacity_species_symbol", &bear::EmissionModelConfig::opacity_species_symbol)
        .def_readwrite("opacity_species_folder", &bear::EmissionModelConfig::opacity_species_folder);

    py::class_<bear::TransmissionPostProcessConfig>(m, "TransmissionPostProcessConfig", genericConfig)
        .def(py::init<
          const std::string&>())
        .def(py::init<
          const bool, 
          const bool, 
          const std::vector<std::string>&>())
        .def_readwrite("save_temperatures", &bear::TransmissionPostProcessConfig::save_temperatures)
        .def_readwrite("save_spectra", &bear::TransmissionPostProcessConfig::save_spectra)
        .def_readwrite("delete_sampler_files", &bear::TransmissionPostProcessConfig::delete_sampler_files);

    py::class_<bear::SecondaryEclipsePostProcessConfig>(m, "SecondaryEclipsePostProcessConfig", genericConfig)
        .def(py::init<
          const std::string&>())
        .def(py::init<
          const bool, 
          const bool, 
          const bool,
          const std::vector<std::string>&>())
        .def_readwrite("save_temperatures", &bear::SecondaryEclipsePostProcessConfig::save_temperatures)
        .def_readwrite("save_spectra", &bear::SecondaryEclipsePostProcessConfig::save_spectra)
        .def_readwrite("save_contribution_functions", &bear::SecondaryEclipsePostProcessConfig::save_contribution_functions)
        .def_readwrite("delete_sampler_files", &bear::SecondaryEclipsePostProcessConfig::delete_sampler_files);

    py::class_<bear::EmissionPostProcessConfig>(m, "EmissionPostProcessConfig", genericConfig)
        .def(py::init<
          const std::string&>())
        .def(py::init<
          const bool, 
          const bool,
          const bool, 
          const bool,
          const std::vector<std::string>&>())
        .def_readwrite("save_temperatures", &bear::EmissionPostProcessConfig::save_temperatures)
        .def_readwrite("save_effective_temperatures", &bear::EmissionPostProcessConfig::save_temperatures)
        .def_readwrite("save_spectra", &bear::EmissionPostProcessConfig::save_spectra)
        .def_readwrite("save_contribution_functions", &bear::EmissionPostProcessConfig::save_contribution_functions)
        .def_readwrite("delete_sampler_files", &bear::EmissionPostProcessConfig::delete_sampler_files);

    py::class_<bear::TransmissionModel>(m, "TransmissionModel")
        .def(py::init<
          bear::GlobalConfig*, 
          bear::SpectralGrid*, 
          const size_t, 
          const std::vector<std::string>&, 
          const std::vector<std::string>&>())
        .def("calcSpectrum", &bear::TransmissionModel::calcSpectrum);

    py::class_<bear::SecondaryEclipseModel>(m, "SecondaryEclipseModel")
        .def(py::init<
          bear::GlobalConfig*, 
          bear::SpectralGrid*, 
          const size_t, 
          const std::vector<double>&, 
          const std::vector<double>&, 
          const std::vector<std::string>&, 
          const std::vector<std::string>&>())
        .def("calcSpectrum", &bear::SecondaryEclipseModel::calcSpectrum);

    py::class_<bear::EmissionModel>(m, "EmissionModel")
        .def(py::init<
          bear::GlobalConfig*, 
          bear::SpectralGrid*, 
          const size_t, 
          const std::string,
          const std::vector<std::string>&,
          const std::vector<std::string>&, 
          const std::vector<std::string>&>())
        .def("calcSpectrum", &bear::EmissionModel::calcSpectrum);
}