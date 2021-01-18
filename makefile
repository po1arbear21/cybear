# colors
FC_COL = \e[1;37m
IN_COL = \e[1;32m
OU_COL = \e[1;36m
NO_COL = \e[m

# load default values for BUILD, COMPILER and library USE flags (overwrite from command line)
include options.mk

# load compiler configuration
ifeq ($(COMPILER),intel)
  include intel.mk
endif
ifeq ($(COMPILER),gnu)
  include gnu.mk
endif

# define macros for fortran
include macro.mk

# main target
all: dirs

# directories
BUILD_DIR := build/
TRASH_DIR := build/trash/
dirs: $(BUILD_DIR) $(TRASH_DIR)
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)
$(TRASH_DIR):
	@mkdir -p $(TRASH_DIR)

# find all source files in src directory
SOURCES := $(shell find src/ -name '*.f90')

# C source files and objects
SOURCES_C := $(shell find src/ -name '*.c')
vpath %.c $(sort $(dir $(SOURCES_C)))
OBJECTS_C := $(addprefix $(BUILD_DIR), $(addsuffix .o, $(notdir $(SOURCES_C))))

# include directories
FINCLUDE := -I$(BUILD_DIR)

# libraries
include lib.mk

# generate targets and their dependencies (one target per program found in sources)
depend: $(BUILD_DIR).depend
$(BUILD_DIR).depend: $(BUILD_DIR) $(SOURCES)
	@printf "%b" "$(FC_COL)python3 depend.py$(NO_COL) $(IN_COL)$(SOURCES)$(NO_COL) -b $(BUILD_DIR) > $(OU_COL)$(BUILD_DIR).depend$(NO_COL) || rm -f $(OU_COL)$(BUILD_DIR).depend$(NO_COL)\n\n"
	@python3 depend.py $(SOURCES) -b $(BUILD_DIR) > $(BUILD_DIR).depend || rm -f $(BUILD_DIR).depend
include $(BUILD_DIR).depend

# build all targets
all: $(TARGETS)

# rule for anchor files
$(BUILD_DIR)%.anc:
	@printf "%b" "$(FC_COL)$(FC)$(NO_COL) $(FFLAGS) $(FINCLUDE) $(FMODULE) $(BUILD_DIR) $(FSYNTAXONLY) -c $(IN_COL)$<$(NO_COL)\n\n"
	@$(FC) $(FFLAGS) $(FINCLUDE) $(FMODULE) $(BUILD_DIR) $(FSYNTAXONLY) -c $<
	@mv $(notdir $(<:.f90=.i90)) $(BUILD_DIR) 2>/dev/null || true
	@touch $@

# rule for c object files
$(BUILD_DIR)%.c.o: %.c
	@printf "%b" "$(FC_COL)$(CC)$(NO_COL) $(CFLAGS) -c $(IN_COL)$<$(NO_COL) -o $(OU_COL)$@$(NO_COL)\n\n"
	@$(CC) $(CFLAGS) -c $< -o $@

# rule for object files
$(BUILD_DIR)%.o:
	@printf "%b" "$(FC_COL)$(FC)$(NO_COL) $(FFLAGS) $(FINCLUDE) $(FMODULE) $(TRASH_DIR) -c $(IN_COL)$<$(NO_COL) -o $(OU_COL)$@$(NO_COL)\n\n"
	@$(FC) $(FFLAGS) $(FINCLUDE) $(FMODULE) $(TRASH_DIR) -c $< -o $@

clean:
	rm -f $(TRASH_DIR)*.{s,}mod $(BUILD_DIR)*.{anc,i90,mod,smod,o} $(BUILD_DIR).depend $(TARGETS)

doc: all
	ford -e i90 -d $(BUILD_DIR) README.md

.PHONY: all dirs doc clean


#
# DOCUMENTATION
#
# steps for compilation
#   1. create module interface files ".mod".
#     - reason for separating compilation from interface file generation: makes parallization of compilation possible.
#    2. create anchor ".anc" files.
#     - reasons
#       1. barrier for parallel compilation s.t. all interface files are created first (no race conditions)
#       2. makes fortran "include" statements possible??
#    3. create object files
#
