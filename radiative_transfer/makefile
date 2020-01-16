
# get options which hold for all files of the project
D=..
include $D/make.global_options


# create lists of file names
c-files     = $(shell echo *.cpp)
o-files-main = $(addprefix ../$(OBJDIR), $(patsubst %,%,$(c-files:.cpp=.$(OBJEXT))))
o-files      = $(sort $(o-files-main))
h-files      = $(sort $(shell echo *.h))
d-files      = $(o-files:.o=.d)


# top level targets
radiative_transfer_: $(o-files)

../$(OBJDIR)%.o : %.cpp
	@echo "compiling: radiative transfer --> $<"
	@$(CXX) $(CXXFLAGS) $< $(OBJ) $@
	
clean: 
	rm -f $(o-files)  $(d-files)

-include $(d-files)
