# DEBUG
ifeq ($(BUILD),debug)
FFLAGS += -DDEBUG
endif

# INTSIZE
ifeq ($(INTSIZE),32)
FFLAGS += -DINTSIZE32
endif
ifeq ($(INTSIZE),64)
FFLAGS += -DINTSIZE64
endif

# IDXSIZE
ifeq ($(IDXSIZE),32)
ifeq ($(INTSIZE),64)
$(error IDXSIZE must be >= INTSIZE)
endif
FFLAGS += -DIDXSIZE32
endif
ifeq ($(IDXSIZE),64)
FFLAGS += -DIDXSIZE64
endif