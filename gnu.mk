# F95 root directory
ifndef GNU_F95ROOT
$(error GNU_F95ROOT is not set)
endif
F95ROOT := $(GNU_F95ROOT)

# gnu fortran compiler flags
FC       := gfortran
FFLAGS   := -cpp -ffree-line-length-none -I./ -march=native -Wall -fopenmp \
            -I$(F95ROOT)/include/intel64/ilp64 -I$(MKLROOT)/include
FDEBUG   := -O0 -g3 -ggdb -fcheck=all -fbacktrace -D DEBUG
FRELEASE := -O3

# additional flags
FINT64           := -fdefault-integer-8
FREAL64          := -fdefault-real-8
FMODULE          := -J
FSYNTAXONLY      := -fsyntax-only
FNOGENINTERFACES :=
FWARNNOUNUSED    := -Wno-unused -Wno-unused-dummy-argument

# gnu c compiler flags
CC       := gcc
CFLAGS   := -march=native
CDEBUG   := -O0 -g
CRELEASE := -O3

# external libraries
EXT_LIBS_DEBUG := \
	$(F95ROOT)/lib/intel64/libmkl_blas95_ilp64.a \
  $(F95ROOT)/lib/intel64/libmkl_lapack95_ilp64.a \
	-Wl,--start-group \
  	$(MKLROOT)/lib/intel64/libmkl_gf_ilp64.a \
  	$(MKLROOT)/lib/intel64/libmkl_sequential.a \
  	$(MKLROOT)/lib/intel64/libmkl_core.a \
	-Wl,--end-group \
	-lgomp -lpthread -lm -ldl

EXT_LIBS_RELEASE := \
	$(F95ROOT)/lib/intel64/libmkl_blas95_ilp64.a \
	$(F95ROOT)/lib/intel64/libmkl_lapack95_ilp64.a \
	-Wl,--start-group \
		$(MKLROOT)/lib/intel64/libmkl_gf_ilp64.a \
		$(MKLROOT)/lib/intel64/libmkl_gnu_thread.a \
		$(MKLROOT)/lib/intel64/libmkl_core.a \
	-Wl,--end-group \
	-lgomp -lpthread -lm -ldl