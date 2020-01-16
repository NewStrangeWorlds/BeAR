
# get options which hold for all files of the project
D=..
include $D/make.global_options


# create lists of file names
c-files     = $(shell echo *.cu)
o-files-main = $(addprefix ../$(OBJDIR), $(patsubst %,%,$(c-files:.cu=.$(OBJEXT))))
o-files      = $(sort $(o-files-main))
h-files      = $(sort $(shell echo *.h))


# top level targets
cuda_kernels_: $(o-files)

../$(OBJDIR)%.o : %.cu
	@echo "compiling: CUDA kernels --> $<"
	@$(CUDA) $(CUDAFLAGS) $< $(OBJ) $@
	
clean: 
	rm -f $(o-files) 
