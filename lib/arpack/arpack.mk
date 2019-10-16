ARPACK_SRC_DIR   := $(dir $(lastword $(MAKEFILE_LIST)))
ARPACK_BUILD_DIR := $(BUILD_DIR)arpack/
ARPACK_SOURCES   := $(shell find $(ARPACK_SRC_DIR) -name '*.f')
ARPACK_OBJECTS   := $(addprefix $(ARPACK_BUILD_DIR), $(notdir $(ARPACK_SOURCES:.f=.o)))
ARPACK_LIB       := $(ARPACK_BUILD_DIR)libarpack.a
ARPACK_FFLAGS    := $(FFLAGS) -nogen-interfaces -warn nounused
LIBS             := $(LIBS) $(ARPACK_LIB)

dirs: $(ARPACK_BUILD_DIR)
$(ARPACK_BUILD_DIR):
	mkdir -p $(ARPACK_BUILD_DIR)

$(ARPACK_LIB): $(ARPACK_OBJECTS)
	$(AR) rcs $(ARPACK_LIB) $(ARPACK_OBJECTS)

$(ARPACK_BUILD_DIR)%.o: $(ARPACK_SRC_DIR)%.f
	$(FC) $(ARPACK_FFLAGS) -c $< -o $@

arpack: $(ARPACK_LIB)

clean_arpack:
	rm -f $(ARPACK_BUILD_DIR)*.o $(ARPACK_BUILD_DIR)*.a

clean_all: clean_arpack

.PHONY: arpack clean_arpack