LIBS :=

# BLAS95 and LAPACK95
ifeq ($(COMPILER),intel)
LIBS += $(MKLROOT)/lib/intel64/libmkl_blas95_ilp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_ilp64.a
endif
ifeq ($(COMPILER),gnu)
ifndef BLAS95ROOT
$(error BLAS95ROOT is not set)
endif
ifndef LAPACK95ROOT
$(error LAPACK95ROOT is not set)
endif
LIBS     += $(BLAS95ROOT)/lib/gnu/libmkl_blas95_ilp64.a $(LAPACK95ROOT)/lib/gnu/libmkl_lapack95_ilp64.a
FINCLUDE += -I$(BLAS95ROOT)/include/gnu -I$(LAPACK95ROOT)/include/gnu
endif

# MKL
ifeq ($(COMPILER),gnu)
ifeq ($(BUILD),debug)
MKL_LIBS := \
	-Wl,--start-group \
  	$(MKLROOT)/lib/intel64/libmkl_gf_ilp64.a \
  	$(MKLROOT)/lib/intel64/libmkl_sequential.a \
  	$(MKLROOT)/lib/intel64/libmkl_core.a \
	-Wl,--end-group
endif
ifneq (,$(filter $(BUILD),release profile))
MKL_LIBS := \
	-Wl,--start-group \
		$(MKLROOT)/lib/intel64/libmkl_gf_ilp64.a \
		$(MKLROOT)/lib/intel64/libmkl_gnu_thread.a \
		$(MKLROOT)/lib/intel64/libmkl_core.a \
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
FFLAGS += -D USE_ARPACK
LIBS   += $(ARPACKROOT)/lib/$(COMPILER)/libarpack.a
endif

# EXPOKIT
ifeq ($(USE_EXPOKIT),true)
ifndef EXPOKITROOT
$(error EXPOKITROOT is not set)
endif
FFLAGS += -D USE_EXPOKIT
LIBS   += $(EXPOKITROOT)/lib/$(COMPILER)/expokit.a
endif

# FEAST
ifeq ($(USE_FEAST),true)
ifndef FEASTROOT
$(error FEASTROOT is not set)
endif
FFLAGS += -D USE_FEAST
LIBS   += $(FEASTROOT)/lib/$(COMPILER)/libfeast.a
endif

# ILUPACK
ifeq ($(USE_ILUPACK),true)
ifndef ILUPACKROOT
$(error ILUPACKROOT is not set)
endif
FFLAGS += -D USE_ILUPACK
LIBS   += \
	-Wl,--start-group \
		$(ILUPACKROOT)/lib/$(COMPILER)/libilupack.a \
		$(ILUPACKROOT)/lib/$(COMPILER)/libamd.a \
		$(ILUPACKROOT)/lib/$(COMPILER)/libblaslike.a \
		$(ILUPACKROOT)/lib/$(COMPILER)/libmetis.a \
		$(ILUPACKROOT)/lib/$(COMPILER)/libsparspak.a \
		$(ILUPACKROOT)/lib/$(COMPILER)/libmumps.a \
	-Wl,--end-group
ifeq ($(COMPILER),gnu)
LIBS += $(MKL_LIBS)
endif
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
FFLAGS += -D USE_MUMPS
LIBS   += \
	-Wl,--start-group \
		$(MUMPSROOT)/lib/$(COMPILER)/libdmumps.a \
		$(MUMPSROOT)/lib/$(COMPILER)/libzmumps.a \
		$(MUMPSROOT)/lib/$(COMPILER)/libmumps_common.a \
		$(MUMPSROOT)/lib/$(COMPILER)/libpord.a \
		$(MUMPSROOT)/lib/$(COMPILER)/libmpiseq.a \
		$(METISROOT)/lib/$(COMPILER)/libmetis.a \
		$(SCOTCHROOT)/lib/$(COMPILER)/libesmumps.a \
		$(SCOTCHROOT)/lib/$(COMPILER)/libscotch.a \
		$(SCOTCHROOT)/lib/$(COMPILER)/libscotcherr.a \
	-Wl,--end-group
FINCLUDE += -I$(MUMPSROOT)/include/$(COMPILER)
endif

# QUADPACK
ifeq ($(USE_QUADPACK),true)
ifndef QUADPACKROOT
$(error QUADPACKROOT is not set)
endif
FFLAGS += -D USE_QUADPACK
LIBS   += $(QUADPACKROOT)/lib/$(COMPILER)/quadpack.a
endif

# additional libraries
ifeq ($(COMPILER), gnu)
LIBS += -lgomp -lpthread -lm -ldl
endif

