
include make.global_options


model_main_: 
	cd $/model_main && $(MAKE) $(MAKEOPTIONS) model_main_
	
	
spectral_grid_: 
	cd $/spectral_grid && $(MAKE) $(MAKEOPTIONS) spectral_grid_
	
global_config_: 
	cd $/config && $(MAKE) $(MAKEOPTIONS) global_config_
	
observations_:
	cd $/observations && $(MAKE) $(MAKEOPTIONS) observations_
	
retrieval_:
	cd $/retrieval && $(MAKE) $(MAKEOPTIONS) retrieval_
	
transport_coeff_:
	cd $/transport_coeff && $(MAKE) $(MAKEOPTIONS) transport_coeff_
	
forward_model_:
	cd $/forward_model && $(MAKE) $(MAKEOPTIONS) forward_model_
	
cdisort_:
	cd $/cdisort && $(MAKE) $(MAKEOPTIONS) cdisort_
	
radiative_transfer_:
	cd $/radiative_transfer && $(MAKE) $(MAKEOPTIONS) radiative_transfer_
	
aux_:
	cd $/additional && $(MAKE) $(MAKEOPTIONS) aux_
	
fastchem_:
#	cd $/fastchem && $(MAKE) $(MAKEOPTIONS) fastchem_
	
cuda_kernels_:
	cd $/CUDA_kernels && $(MAKE) $(MAKEOPTIONS) cuda_kernels_
	
all-c: model_main_ spectral_grid_ global_config_ observations_ retrieval_ transport_coeff_ forward_model_ cdisort_ aux_ radiative_transfer_ fastchem_ cuda_kernels_

	@echo "compiling done."

all-l:
	cd $/obj && $(MAKE) $(MAKEOPTIONS) all-l
	
	@echo "linking done."
	
all: all-c all-l

	@echo "everything is done and fine. enjoy your day!"

	
#and here we clean the mess
clean: clean-binary clean-model_main clean-spectral_grid clean-global_config clean-observations clean-retrieval clean-transport_coeff clean-cuda_kernels clean-forward_model clean-cdisort clean-radiative_transfer clean-fastchem clean-aux

clean-binary: 
	rm -f $(EXECUTABLE_NAME)
	
clean-model_main: 
	cd $/model_main && $(MAKE) clean
	
clean-global_config: 
	cd $/config && $(MAKE) clean

clean-spectral_grid: 
	cd $/spectral_grid && $(MAKE) clean
	
clean-observations: 
	cd $/observations && $(MAKE) clean
	
clean-retrieval: 
	cd $/retrieval && $(MAKE) clean
	
clean-transport_coeff: 
	cd $/transport_coeff && $(MAKE) clean
	
clean-forward_model: 
	cd $/forward_model && $(MAKE) clean
	
clean-cdisort: 
	cd $/cdisort && $(MAKE) clean
	
clean-radiative_transfer: 
	cd $/radiative_transfer && $(MAKE) clean
	
clean-aux: 
	cd $/additional && $(MAKE) clean
	
clean-fastchem: 
#	cd $/fastchem && $(MAKE) clean
	
clean-cuda_kernels: 
	cd $/CUDA_kernels && $(MAKE) clean
	
	
	@echo "all clean. have a nice day!"
	
