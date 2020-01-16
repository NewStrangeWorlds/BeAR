
# get options which hold for all files of the project
D=..
include $D/make.global_options


# create lists of file names
c-files     = $(shell echo *.c)
o-files-main = $(addprefix ../$(OBJDIR), $(patsubst %,%,$(c-files:.c=.$(OBJEXT))))
o-files      = $(sort $(o-files-main))
h-files      = $(sort $(shell echo *.h))
d-files      = $(o-files:.o=.d)


# top level targets
cdisort_: $(o-files)

../$(OBJDIR)%.o : %.c
	@echo "compiling: C-Disort --> $<"
	@$(CC) $(CCFLAGS) $< $(OBJ) $@
	
clean: 
	rm -f $(o-files) $(d-files)

-include $(d-files)
