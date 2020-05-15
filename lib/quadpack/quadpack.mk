QUADPACK_SRC_DIR   := $(dir $(lastword $(MAKEFILE_LIST)))
QUADPACK_BUILD_DIR := $(BUILD_DIR)quadpack/
QUADPACK_SOURCES   := $(shell find $(QUADPACK_SRC_DIR) -name '*.f')
QUADPACK_OBJECTS   := $(addprefix $(QUADPACK_BUILD_DIR), $(notdir $(QUADPACK_SOURCES:.f=.o)))
QUADPACK_LIB       := $(QUADPACK_BUILD_DIR)libquadpack.a
QUADPACK_FFLAGS    := $(FFLAGS) -nogen-interfaces -warn nounused
LIBS               := $(LIBS) $(QUADPACK_LIB)

dirs: $(QUADPACK_BUILD_DIR)
$(QUADPACK_BUILD_DIR):
	@mkdir -p $(QUADPACK_BUILD_DIR)

$(QUADPACK_LIB): $(QUADPACK_OBJECTS)
	@printf "%b" "$(FC_COL)$(AR)$(NO_COL) rcs $(OU_COL)$(QUADPACK_LIB)$(NO_COL) $(IN_COL)$(QUADPACK_OBJECTS)$(NO_COL)\n\n"
	@$(AR) rcs $(QUADPACK_LIB) $(QUADPACK_OBJECTS)

$(QUADPACK_BUILD_DIR)%.o: $(QUADPACK_SRC_DIR)%.f
	@printf "%b" "$(FC_COL)$(FC)$(NO_COL) $(QUADPACK_FFLAGS) -c $(IN_COL)$<$(NO_COL) -o $(OU_COL)$@$(NO_COL)\n\n"
	@$(FC) $(QUADPACK_FFLAGS) -c $< -o $@

quadpack: $(QUADPACK_LIB)

clean_quadpack:
	rm -f $(QUADPACK_BUILD_DIR)*.o $(QUADPACK_BUILD_DIR)*.a

clean_all: clean_quadpack

.PHONY: quadpack clean_quadpack