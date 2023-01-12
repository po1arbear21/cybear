LIBS     :=
FINCLUDE :=
CINCLUDE :=

ifeq ($(INTSIZE),32)
  ARCH := lp64
endif
ifeq ($(INTSIZE),64)
  ARCH := ilp64
endif

# FEAST
ifeq ($(USE_FEAST),true)
  ifndef FEASTROOT
    $(error FEASTROOT is not set)
  endif
  LIBS += $(FEASTROOT)/lib/$(COMPILER)_$(ARCH)/libfeast.a
endif

# BLAS95 and LAPACK95
ifeq ($(COMPILER),intel)
  LIBS += $(MKLROOT)/lib/intel64/libmkl_blas95_$(ARCH).a $(MKLROOT)/lib/intel64/libmkl_lapack95_$(ARCH).a
endif
ifeq ($(COMPILER),gnu)
  ifndef BLAS95ROOT
    $(error BLAS95ROOT is not set)
  endif
  ifndef LAPACK95ROOT
    $(error LAPACK95ROOT is not set)
  endif
  LIBS     += $(BLAS95ROOT)/lib/gnu_$(ARCH)/libmkl_blas95_$(ARCH).a $(LAPACK95ROOT)/lib/gnu_$(ARCH)/libmkl_lapack95_$(ARCH).a
  FINCLUDE += -I$(BLAS95ROOT)/include/gnu_$(ARCH) -I$(LAPACK95ROOT)/include/gnu_$(ARCH)
endif

# MKL
ifeq ($(COMPILER),gnu)
  ifeq ($(BUILD),debug)
    MKL_LIBS :=                                    \
      -Wl,--start-group                            \
        $(MKLROOT)/lib/intel64/libmkl_gf_$(ARCH).a \
        $(MKLROOT)/lib/intel64/libmkl_sequential.a \
        $(MKLROOT)/lib/intel64/libmkl_core.a       \
      -Wl,--end-group
  endif
  ifneq (,$(filter $(BUILD),release profile))
    MKL_LIBS :=                                    \
      -Wl,--start-group                            \
        $(MKLROOT)/lib/intel64/libmkl_gf_$(ARCH).a \
        $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a \
        $(MKLROOT)/lib/intel64/libmkl_core.a       \
      -Wl,--end-group
  endif
  LIBS     += $(MKL_LIBS)
  FINCLUDE += -I$(MKLROOT)/include
endif

# ARPACK
ifeq ($(USE_ARPACK),true)
  ifndef ARPACKROOT
    $(error ARPACKROOT is not set)
  endif
  LIBS += $(ARPACKROOT)/lib/$(COMPILER)_$(ARCH)/libarpack.a
endif

# EXPOKIT
ifeq ($(USE_EXPOKIT),true)
  ifndef EXPOKITROOT
    $(error EXPOKITROOT is not set)
  endif
  LIBS += $(EXPOKITROOT)/lib/$(COMPILER)_$(ARCH)/expokit.a
endif

# ILUPACK
ifeq ($(USE_ILUPACK),true)
  ifndef ILUPACKROOT
    $(error ILUPACKROOT is not set)
  endif

  # ILUPACK does not support mixed-precision integers. either 32bit or 64bit.
  ifeq ($(IDXSIZE),64)
    ARCH_ILUPACK := ilp64
  else
    ARCH_ILUPACK := lp64
  endif

  LIBS +=                                                                   \
    -Wl,--start-group                                                       \
      $(ILUPACKROOT)/lib/$(COMPILER)_$(ARCH_ILUPACK)/libilupack_mumps.a     \
      $(ILUPACKROOT)/lib/$(COMPILER)_$(ARCH_ILUPACK)/libamd.a               \
      $(ILUPACKROOT)/lib/$(COMPILER)_$(ARCH_ILUPACK)/libblaslike.a          \
      $(ILUPACKROOT)/lib/$(COMPILER)_$(ARCH_ILUPACK)/libcamd.a              \
      $(ILUPACKROOT)/lib/$(COMPILER)_$(ARCH_ILUPACK)/libmetis.a             \
      $(ILUPACKROOT)/lib/$(COMPILER)_$(ARCH_ILUPACK)/libmetisomp.a          \
      $(ILUPACKROOT)/lib/$(COMPILER)_$(ARCH_ILUPACK)/libsparspak.a          \
      $(ILUPACKROOT)/lib/$(COMPILER)_$(ARCH_ILUPACK)/libsuitesparseconfig.a \
      $(ILUPACKROOT)/lib/$(COMPILER)_$(ARCH_ILUPACK)/libmumps.a             \
    -Wl,--end-group
  ifeq ($(COMPILER),gnu)
    LIBS += $(MKL_LIBS)
  endif

  undefine ARCH_ILUPACK
endif

# MUMPS
ifeq ($(USE_MUMPS),true)
  ifndef MUMPSROOT
    $(error MUMPSROOT is not set)
  endif
  ifndef METISROOT
    $(error METISROOT is not set)
  endif
  ifndef SCOTCHROOT
    $(error SCOTCHROOT is not set)
  endif
  LIBS +=                                                    \
    -Wl,--start-group                                        \
      $(MUMPSROOT)/lib/$(COMPILER)_$(ARCH)/libdmumps.a       \
      $(MUMPSROOT)/lib/$(COMPILER)_$(ARCH)/libzmumps.a       \
      $(MUMPSROOT)/lib/$(COMPILER)_$(ARCH)/libmumps_common.a \
      $(MUMPSROOT)/lib/$(COMPILER)_$(ARCH)/libpord.a         \
      $(MUMPSROOT)/lib/$(COMPILER)_$(ARCH)/libmpiseq.a       \
      $(METISROOT)/lib/$(COMPILER)_ilp64/libmetis.a          \
      $(SCOTCHROOT)/lib/$(COMPILER)_ilp64/libesmumps.a       \
      $(SCOTCHROOT)/lib/$(COMPILER)_ilp64/libscotch.a        \
      $(SCOTCHROOT)/lib/$(COMPILER)_ilp64/libscotcherr.a     \
    -Wl,--end-group
  ifeq ($(COMPILER),gnu)
    LIBS += $(MKL_LIBS)
  endif
  FINCLUDE += -I$(MUMPSROOT)/include/$(COMPILER)_$(ARCH)
endif

# QUADPACK
ifeq ($(USE_QUADPACK),true)
  ifndef QUADPACKROOT
    $(error QUADPACKROOT is not set)
  endif
  LIBS += $(QUADPACKROOT)/lib/$(COMPILER)_$(ARCH)/quadpack.a
endif

# MPFR
ifeq ($(USE_MPFR),true)
  ifndef GMPROOT
    $(error GMPROOT is not set)
  endif
  ifndef MPFRROOT
    $(error MPFRROOT is not set)
  endif
  LIBS += $(MPFRROOT)/lib/$(COMPILER)_ilp64/libmpfr.a $(GMPROOT)/lib/$(COMPILER)_ilp64/libgmp.a
  CINCLUDE += -I$(MPFRROOT)/include/$(COMPILER)_ilp64 -I$(GMPROOT)/include/$(COMPILER)_ilp64
endif

ifeq ($(USE_TRIANGLE),true)
  LIBS += $(TRIANGLEROOT)/lib/$(COMPILER)_lp64/triangle.a
endif

# additional libraries
ifeq ($(COMPILER), gnu)
  LIBS += -lgomp -lpthread -lm -ldl
endif
