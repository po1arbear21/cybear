# gnu fortran compiler flags
FC     := gfortran
FFLAGS := -cpp -ffree-line-length-none -I./ -march=native -Wall -fopenmp -fuse-ld=bfd
ifeq ($(BUILD),debug)
FFLAGS += -O0 -g3 -ggdb -fcheck=all -fbacktrace -D DEBUG
endif
ifeq ($(BUILD),release)
FFLAGS += -O2
endif
ifeq ($(BUILD),profile)
FFLAGS += -O2 -g -shared-libgcc
endif

# additional fortran flags
FINT64           := -fdefault-integer-8
FREAL64          := -fdefault-real-8
FMODULE          := -J
FSYNTAXONLY      := -fsyntax-only
FNOGENINTERFACES :=
FWARNNOUNUSED    := -Wno-unused -Wno-unused-dummy-argument

# gnu c compiler flags
CC     := gcc
CFLAGS := -march=native
ifeq ($(BUILD),debug)
CFLAGS += -O0 -g
endif
ifeq ($(BUILD),release)
CFLAGS += -O3
endif
ifeq ($(BUILD),profile)
CFLAGS += -O3 -g -shared-libgcc
endif
