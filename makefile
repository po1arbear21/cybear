# use mumps ?
ifeq ($(USE_MUMPS),true)
ifndef MUMPSROOT
$(error MUMPSROOT is not set)
endif
ifndef METISROOT
$(error METISROOT is not set)
endif
endif

# use ilupack ?
ifeq ($(USE_ILUPACK),true)
ifndef ILUPACKROOT
$(error ILUPACKROOT is not set)
endif
endif

# compiler configuration (overwrite from command line)
COMPILER := intel
ifeq ($(COMPILER),intel)
include intel.mk
else
ifeq ($(COMPILER),gnu)
include gnu.mk
else
$(error COMPILER must be intel or gnu!)
endif
endif

# build configuration (overwrite from command line)
BUILD := debug
ifeq ($(BUILD),debug)
FFLAGS   := $(FFLAGS) $(FDEBUG)
CFLAGS   := $(CFLAGS) $(CDEBUG)
EXT_LIBS := $(EXT_LIBS_DEBUG)
else
ifeq ($(BUILD),release)
FFLAGS   := $(FFLAGS) $(FRELEASE)
CFLAGS   := $(CFLAGS) $(CRELEASE)
EXT_LIBS := $(EXT_LIBS_RELEASE)
else
ifeq ($(BUILD),profile)
FFLAGS   := $(FFLAGS) $(FPROFILE)
CFLAGS   := $(CFLAGS) $(CPROFILE)
EXT_LIBS := $(EXT_LIBS_PROFILE)
else
$(error BUILD must be debug, release or profile!)
endif
endif
endif

# colors
FC_COL = \e[1;37m
IN_COL = \e[1;32m
OU_COL = \e[1;36m
NO_COL = \e[m

# main target
all: dirs

# directories
BUILD_DIR := build/${COMPILER}/${BUILD}/
TRASH_DIR := build/${COMPILER}/${BUILD}/trash/
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

# additional libraries
include lib/quadpack/quadpack.mk
ifeq ($(COMPILER),gnu)

# needs to be outside blas95,lapack95 makefiles, otherwise makefile warnings.
F95_BUILD_DIR := build/${COMPILER}/F95

$(F95_BUILD_DIR):
	@mkdir -p $(F95_BUILD_DIR)

include lib/blas95/blas95.mk
include lib/lapack95/lapack95.mk
endif

# include directories
FINCLUDE = -I$(BUILD_DIR) $(MUMPS_INC)

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
	@printf "%b" "$(FC_COL)$(FC)$(NO_COL) $(FFLAGS) $(FINT64) $(FREAL64) $(FINCLUDE) $(FMODULE) $(BUILD_DIR) $(FSYNTAXONLY) -c $(IN_COL)$<$(NO_COL)\n\n"
	@$(FC) $(FFLAGS) $(FINT64) $(FREAL64) $(FINCLUDE) $(FMODULE) $(BUILD_DIR) $(FSYNTAXONLY) -c $<
	@mv $(notdir $(<:.f90=.i90)) $(BUILD_DIR) 2>/dev/null || true
	@touch $@

# rule for c object files
$(BUILD_DIR)%.c.o: %.c
	@printf "%b" "$(FC_COL)$(CC)$(NO_COL) $(CFLAGS) -c $(IN_COL)$<$(NO_COL) -o $(OU_COL)$@$(NO_COL)\n\n"
	@$(CC) $(CFLAGS) -c $< -o $@

# rule for object files
$(BUILD_DIR)%.o:
	@printf "%b" "$(FC_COL)$(FC)$(NO_COL) $(FFLAGS) $(FINT64) $(FREAL64) $(FINCLUDE) $(FMODULE) $(TRASH_DIR) -c $(IN_COL)$<$(NO_COL) -o $(OU_COL)$@$(NO_COL)\n\n"
	@$(FC) $(FFLAGS) $(FINT64) $(FREAL64) $(FINCLUDE) $(FMODULE) $(TRASH_DIR) -c $< -o $@

clean:
	rm -f $(TRASH_DIR)*.{s,}mod $(BUILD_DIR)*.{anc,i90,mod,smod,o} $(BUILD_DIR).depend $(TARGETS)

doc: all
	ford -e i90 -d $(BUILD_DIR) README.md

clean_all: clean

.PHONY: all dirs doc clean clean_all


#
# DOCUMENTATION
#
# steps for compilation
# 	1. create module interface files ".mod".
# 		- reason for separating compilation from interface file generation: makes parallization of compilation possible.
#		2. create anchor ".anc" files.
# 		- reasons
# 			1. barrier for parallel compilation s.t. all interface files are created first (no race conditions)
# 			2. makes fortran "include" statements possible??
#		3. create object files
#
