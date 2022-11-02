M4 := m4

# force m4_ prefix for all builtin macros
M4FLAGS := -P

# generate synclines
M4FLAGS += -s

# DEBUG
ifeq ($(BUILD),debug)
  M4FLAGS += -Dm4_debug
endif

# INTSIZE
ifneq ($(INTSIZE),32)
  ifneq ($(INTSIZE),64)
    $(error INTSIZE must be 32 or 64)
  endif
endif
M4FLAGS += -Dm4_intsize=$(INTSIZE)

# IDXSIZE
ifneq ($(IDXSIZE),32)
  ifneq ($(IDXSIZE),64)
    $(error IDXSIZE must be 32 or 64)
  endif
endif
ifeq ($(IDXSIZE),32)
  ifeq ($(INTSIZE),64)
    $(error IDXSIZE must be >= INTSIZE)
  endif
endif
M4FLAGS += -Dm4_idxsize=$(IDXSIZE)

# compiler flags
ifeq ($(COMPILER),intel)
  M4FLAGS += -Dm4_intel
endif
ifeq ($(COMPILER),gnu)
  M4FLAGS += -Dm4_gnu
endif

# library flags
ifeq ($(USE_FEAST),true)
  M4FLAGS += -Dm4_feast
endif
ifeq ($(USE_ARPACK),true)
  M4FLAGS += -Dm4_arpack
endif
ifeq ($(USE_EXPOKIT),true)
  M4FLAGS += -Dm4_expokit
endif
ifeq ($(USE_ILUPACK),true)
  M4FLAGS += -Dm4_ilupack
endif
ifeq ($(USE_MUMPS),true)
  M4FLAGS += -Dm4_mumps
endif
ifeq ($(USE_QUADPACK),true)
  M4FLAGS += -Dm4_quadpack
endif
ifeq ($(USE_MPFR),true)
  M4FLAGS += -Dm4_mpfr
endif
