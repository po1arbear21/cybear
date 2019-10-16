# compiler flags
FC       := ifort
FFLAGS   := -march=native -real-size 64 -i8 -fpp -warn all -qopenmp
FDEBUG   := -O0 -mkl=sequential -g -check all -check noarg_temp_created -fpe1 -traceback -debug extended -D DEBUG
FRELEASE := -O2 -mkl -ftz

# build configuration (can be overwritten from command line)
BUILD := debug
ifeq ($(BUILD),debug)
FFLAGS := $(FFLAGS) $(FDEBUG)
else
ifeq ($(BUILD),release)
FFLAGS := $(FFLAGS) $(FRELEASE)
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
	mkdir -p $(BUILD_DIR)
$(TRASH_DIR):
	mkdir -p $(TRASH_DIR)

# find all source files in src directory
SOURCES := $(shell find src/ -name '*.f90')

# libraries
LIBS := ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a
include lib/arpack/arpack.mk
#include lib/expokit/expokit.mk
#include lib/feast/feast.mk
#include lib/ilupack/ilupack.mk
#include lib/quadpack/quadpack.mk

# generate targets and their dependencies (one target per program found in sources)
depend: $(BUILD_DIR).depend
$(BUILD_DIR).depend: $(BUILD_DIR) $(SOURCES)
	./depend/depend $(SOURCES) -b $(BUILD_DIR) > $(BUILD_DIR).depend || rm -f $(BUILD_DIR).depend
include $(BUILD_DIR).depend

# build all targets
all: $(TARGETS)

# rule for anchor files
$(BUILD_DIR)%.anc:
	$(FC) $(FFLAGS) -module $(BUILD_DIR) -syntax-only -c $<
	@mv $(notdir $(<:.f90=.i90)) $(BUILD_DIR) 2>/dev/null || true
	@touch $@

# rule for object files
$(BUILD_DIR)%.o:
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -module $(TRASH_DIR) -c $< -o $@

clean:
	rm -f $(TRASH_DIR)*.mod $(BUILD_DIR)*.anc $(BUILD_DIR)*.i90 $(BUILD_DIR)*.mod $(BUILD_DIR)*.o $(BUILD_DIR).depend $(TARGETS)

clean_all: clean

.PHONY: all dirs clean clean_all