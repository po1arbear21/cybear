# intel fortran compiler flags
FC     := ifort
FFLAGS := -march=native -fpp -warn all -qopenmp -fp-model precise
ifeq ($(BUILD),debug)
FFLAGS += -O0 -mkl=sequential -g -check all -check noarg_temp_created -fpe1 -traceback -debug extended -D DEBUG -init=huge
endif
ifeq ($(BUILD),release)
FFLAGS += -O3 -mkl -ftz
endif
ifeq ($(BUILD),profile)
FFLAGS += -O3 -mkl -ftz -g -shared-intel -debug inline-debug-info -parallel-source-info=2
endif

# additional fortran flags
FINT64           := -i8
FREAL64          := -real-size 64
FMODULE          := -module
FSYNTAXONLY      := -syntax-only
FNOGENINTERFACES := -nogen-interfaces
FWARNNOUNUSED    := -warn nounused

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
