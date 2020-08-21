# intel fortran compiler flags
FC       := ifort
FFLAGS   := -march=native -fpp -warn all -qopenmp -fp-model precise
FDEBUG   := -O0 -mkl=sequential -g -check all -check noarg_temp_created -fpe1 -traceback -debug extended -D DEBUG -init=huge
FRELEASE := -O3 -mkl -ftz
FPROFILE := -O3 -mkl -ftz -g -shared-intel -debug inline-debug-info -parallel-source-info=2

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
CPROFILE := -O3 -g -shared-intel

# ILUPACK
ifeq ($(USE_ILUPACK),true)
FFLAGS := $(FFLAGS) -D USE_ILUPACK

ILUPACK_LIBS := \
	-Wl,--start-group \
		$(ILUPACKROOT)/INTEL64_long/libilupack.a \
		$(ILUPACKROOT)/INTEL64_long/libamd.a \
		$(ILUPACKROOT)/INTEL64_long/libblaslike.a \
		$(ILUPACKROOT)/INTEL64_long/libmetis.a \
		$(ILUPACKROOT)/INTEL64_long/libsparspak.a \
		$(ILUPACKROOT)/INTEL64_long/libmumps.a \
	-Wl,--end-group
else
ILUPACK_LIBS :=
endif

# external libraries
EXT_LIBS_DEBUG := \
	$(MKLROOT)/lib/intel64/libmkl_blas95_ilp64.a \
  $(MKLROOT)/lib/intel64/libmkl_lapack95_ilp64.a \
	$(ILUPACK_LIBS)

EXT_LIBS_RELEASE := \
	$(MKLROOT)/lib/intel64/libmkl_blas95_ilp64.a \
	$(MKLROOT)/lib/intel64/libmkl_lapack95_ilp64.a \
	$(ILUPACK_LIBS)

EXT_LIBS_PROFILE := $(EXT_LIBS_RELEASE)
