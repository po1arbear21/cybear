# gnu fortran compiler flags
FC     := gfortran
FFLAGS := -cpp -ffree-line-length-none -I./ -march=native -Wall -Wno-maybe-uninitialized -fopenmp -fuse-ld=bfd -fdefault-real-8
ifeq ($(BUILD),debug)
  FFLAGS += -O0 -g3 -ggdb -fcheck=all -fbacktrace -finit-integer=-9999 -finit-real=-inf
endif
ifeq ($(BUILD),release)
  FFLAGS += -O2 -flto
endif
ifeq ($(BUILD),profile)
  FFLAGS += -O2 -g -shared-libgcc -flto
endif
ifeq ($(INTSIZE),64)
  FFLAGS += -fdefault-integer-8
endif

# additional fortran flags
FMODULE     := -J
FSYNTAXONLY := -fsyntax-only

# gnu c compiler flags
CC     := gcc
CFLAGS := -march=native
ifeq ($(BUILD),debug)
  CFLAGS += -O0 -g
endif
ifeq ($(BUILD),release)
  CFLAGS += -O3 -flto
endif
ifeq ($(BUILD),profile)
  CFLAGS += -O3 -g -shared-libgcc -flto
endif
