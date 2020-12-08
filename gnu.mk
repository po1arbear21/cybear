# gnu fortran compiler flags
FC       := gfortran
FFLAGS   := -cpp -ffree-line-length-none -I./ -march=native -I./build/gnu/F95 -Wall -fopenmp -I$(MKLROOT)/include -fuse-ld=bfd
FDEBUG   := -O0 -g3 -ggdb -fcheck=all -fbacktrace -D DEBUG
FRELEASE := -O2
FPROFILE := -O2 -g -shared-libgcc

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
CPROFILE := -O3 -g -shared-libgcc

# ILUPACK
ifeq ($(USE_ILUPACK),true)
FFLAGS := $(FFLAGS) -D USE_ILUPACK

ILUPACK_LIBS := \
	-Wl,--start-group \
		$(ILUPACKROOT)/GNU64_long/libilupack.a \
		$(ILUPACKROOT)/GNU64_long/libamd.a \
		$(ILUPACKROOT)/GNU64_long/libblaslike.a \
		$(ILUPACKROOT)/GNU64_long/libmetis.a \
		$(ILUPACKROOT)/GNU64_long/libsparspak.a \
		$(ILUPACKROOT)/GNU64_long/libmumps.a \
	-Wl,--end-group
else
ILUPACK_LIBS :=
endif

# external libraries
EXT_LIBS_DEBUG := \
	-Wl,--start-group \
  	$(MKLROOT)/lib/intel64/libmkl_gf_ilp64.a \
  	$(MKLROOT)/lib/intel64/libmkl_sequential.a \
  	$(MKLROOT)/lib/intel64/libmkl_core.a \
	-Wl,--end-group \
	$(ILUPACK_LIBS) \
	-Wl,--start-group \
  	$(MKLROOT)/lib/intel64/libmkl_gf_ilp64.a \
  	$(MKLROOT)/lib/intel64/libmkl_sequential.a \
  	$(MKLROOT)/lib/intel64/libmkl_core.a \
	-Wl,--end-group \
	-lgomp -lpthread -lm -ldl \

EXT_LIBS_RELEASE := \
	-Wl,--start-group \
		$(MKLROOT)/lib/intel64/libmkl_gf_ilp64.a \
		$(MKLROOT)/lib/intel64/libmkl_gnu_thread.a \
		$(MKLROOT)/lib/intel64/libmkl_core.a \
	-Wl,--end-group \
	$(ILUPACK_LIBS) \
	-Wl,--start-group \
		$(MKLROOT)/lib/intel64/libmkl_gf_ilp64.a \
		$(MKLROOT)/lib/intel64/libmkl_gnu_thread.a \
		$(MKLROOT)/lib/intel64/libmkl_core.a \
	-Wl,--end-group \
	-lgomp -lpthread -lm -ldl \

EXT_LIBS_PROFILE := $(EXT_LIBS_RELEASE)
