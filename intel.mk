# intel fortran compiler flags, removed "-march=native" from flags
FC     := ifort
FFLAGS := -warn all -qopenmp -fp-model precise -real-size 64 -fpp
ifeq ($(BUILD),debug)
  FFLAGS += -O0 -qmkl=sequential -g -check all -check noarg_temp_created -fpe1 -traceback -debug extended -init=huge
endif
ifeq ($(BUILD),release)
  FFLAGS += -O3 -qmkl -ftz
endif
ifeq ($(BUILD),profile)
  FFLAGS += -O3 -qmkl -ftz -g -shared-intel -debug inline-debug-info -parallel-source-info=2
endif
ifeq ($(INTSIZE),64)
  FFLAGS += -i8
endif

# additional fortran flags
FMODULE     := -module
FSYNTAXONLY := -syntax-only

# intel c compiler flags
CC     := icc
CFLAGS := -march=native
ifeq ($(BUILD),debug)
  CFLAGS += -O0 -g
endif
ifeq ($(BUILD),release)
  CFLAGS += -O3
endif
ifeq ($(BUILD),profile)
  CFLAGS += -O3 -g -shared-intel
endif
