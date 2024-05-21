
#ifndef _contribution_function_kernels_h
#define _contribution_function_kernels_h


#include <vector>



namespace helios{


void contributionFunctionGPU(double* contribution_function_dev,
                             double* absorption_coeff_device, double* wavenumber_list_device,
                             std::vector<double>& temperature, std::vector<double>& vertical_grid,
                             const size_t& nb_spectral_points);

}


#endif

