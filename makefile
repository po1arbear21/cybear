# fortran compiler flags
FC       := ifort
FFLAGS   := -march=native -real-size 64 -i8 -fpp -warn all -qopenmp -fp-model precise
FDEBUG   := -O0 -mkl=sequential -g -check all -check noarg_temp_created -fpe1 -traceback -debug extended -D DEBUG -init=huge
FRELEASE := -O3 -mkl -ftz

# c compiler flags
CC       := icc
CFLAGS   := -march=native
CDEBUG   := -O0 -g
CRELEASE := -O3

# colors
FC_COL = \e[1;37m
IN_COL = \e[1;32m
OU_COL = \e[1;36m
NO_COL = \e[m

# build configuration (can be overwritten from command line)
BUILD := debug
ifeq ($(BUILD),debug)
FFLAGS := $(FFLAGS) $(FDEBUG)
CFLAGS := $(CFLAGS) $(CDEBUG)
else
ifeq ($(BUILD),release)
FFLAGS := $(FFLAGS) $(FRELEASE)
CFLAGS := $(CFLAGS) $(CRELEASE)
else
$(error BUILD must be debug or release!)
endif
endif

# main target
all: dirs

# directories
BUILD_DIR := build/${BUILD}/
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

# libraries
LIBS := ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a
include lib/arpack/arpack.mk
include lib/expokit/expokit.mk
#include lib/feast/feast.mk
include lib/quadpack/quadpack.mk

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
	@printf "%b" "$(FC_COL)$(FC)$(NO_COL) $(FFLAGS) -module $(BUILD_DIR) -syntax-only -c $(IN_COL)$<$(NO_COL)\n\n"
	@$(FC) $(FFLAGS) -module $(BUILD_DIR) -syntax-only -c $<
	@mv $(notdir $(<:.f90=.i90)) $(BUILD_DIR) 2>/dev/null || true
	@touch $@

# rule for c object files
$(BUILD_DIR)%.c.o: %.c
	@printf "%b" "$(FC_COL)$(CC)$(NO_COL) $(CFLAGS) -c $(IN_COL)$<$(NO_COL) -o $(OU_COL)$@$(NO_COL)\n\n"
	@$(CC) $(CFLAGS) -c $< -o $@

# rule for object files
$(BUILD_DIR)%.o:
	@printf "%b" "$(FC_COL)$(FC)$(NO_COL) $(FFLAGS) -I$(BUILD_DIR) -module $(TRASH_DIR) -c $(IN_COL)$<$(NO_COL) -o $(OU_COL)$@$(NO_COL)\n\n"
	@$(FC) $(FFLAGS) -I$(BUILD_DIR) -module $(TRASH_DIR) -c $< -o $@

clean:
	rm -f $(TRASH_DIR)*.{s,}mod $(BUILD_DIR)*.{anc,i90,mod,smod,o} $(BUILD_DIR).depend $(TARGETS)

doc: all
	ford -e i90 -d $(BUILD_DIR) README.md

clean_all: clean

.PHONY: all dirs doc clean clean_all
