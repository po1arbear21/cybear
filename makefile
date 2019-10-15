FC        := ifort
FFLAGS    := -march=native -real-size 64 -i8 -fpp -warn all -qopenmp
FDEBUG    := -O0 -mkl=sequential -g -check all -check noarg_temp_created -fpe1 -traceback -debug extended -D DEBUG
FRELEASE  := -O2 -mkl -ftz
LIBS      := ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a
BUILD     := debug
BUILD_DIR := build/${BUILD}/
TRASH_DIR := build/trash/

ifeq ($(BUILD),debug)
	FFLAGS := $(FFLAGS) $(FDEBUG)
endif
ifeq ($(BUILD),release)
	FFLAGS := $(FFLAGS) $(FRELEASE)
endif

SOURCES = $(shell find src/ -name '*.f90')

all:

depend: $(BUILD_DIR).depend

$(BUILD_DIR).depend: $(SOURCES)
	./depend/depend $(SOURCES) -b $(BUILD_DIR) > $(BUILD_DIR).depend || rm -f $(BUILD_DIR).depend

include $(BUILD_DIR).depend

all: $(TARGETS)

$(BUILD_DIR)%.anc:
	$(FC) $(FFLAGS) -module $(BUILD_DIR) -syntax-only -c $<
	@mv $(notdir $(<:.f90=.i90)) $(BUILD_DIR) 2>/dev/null || true
	@touch $@

$(BUILD_DIR)%.o:
	$(FC) $(FFLAGS) -I$(BUILD_DIR) -module $(TRASH_DIR) -c $< -o $@

clean:
	rm -f $(TRASH_DIR)*.mod $(BUILD_DIR)*.anc $(BUILD_DIR)*.i90 $(BUILD_DIR)*.mod $(BUILD_DIR)*.o $(BUILD_DIR).depend $(TARGETS)

.PHONY : clean all

# #-----------------------------------------------------------------------------------------------------------------------
# # Additional LIBRARIES
# #-----------------------------------------------------------------------------------------------------------------------

# # arpack
# LIBS := $(LIBS) lib/arpack/libarpack.a
# INCLUDE := $(INCLUDE) -Ilib/arpack/

# # expokit
# LIBS := $(LIBS) lib/expokit/libexpokit.a
# INCLUDE := $(INCLUDE) -Ilib/expokit/

# # feast

# #-----------------------------------------------------------------------------------------------------------------------