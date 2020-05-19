# F95 root directory
F95ROOT := $(MKLROOT)

# intel fortran compiler flags
FC       := ifort
FFLAGS   := -march=native -real-size 64 -i8 -fpp -warn all -qopenmp -fp-model precise
FDEBUG   := -O0 -mkl=sequential -g -check all -check noarg_temp_created -fpe1 -traceback -debug extended -D DEBUG -init=huge
FRELEASE := -O3 -mkl -ftz

# additional flags
FINT64           := -i8
FREAL64          := -real-size 64
FMODULE          := -module
FSYNTAXONLY      := -syntax-only
FNOGENINTERFACES := -nogen-interfaces
FWARNNOUNUSED    := -warn nounused

# intel c compiler flags
CC       := icc
CFLAGS   := -march=native
CDEBUG   := -O0 -g
CRELEASE := -O3

# external libraries
EXT_LIBS_DEBUG := \
	$(F95ROOT)/lib/intel64/libmkl_blas95_ilp64.a \
  $(F95ROOT)/lib/intel64/libmkl_lapack95_ilp64.a

EXT_LIBS_RELEASE := \
	$(F95ROOT)/lib/intel64/libmkl_blas95_ilp64.a \
	$(F95ROOT)/lib/intel64/libmkl_lapack95_ilp64.a
