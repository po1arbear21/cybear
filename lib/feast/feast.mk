FEAST_SRC_DIR   := $(dir $(lastword $(MAKEFILE_LIST)))
FEAST_BUILD_DIR := $(BUILD_DIR)feast/
FEAST_SOURCES   := $(shell find $(FEAST_SRC_DIR) -name '*.f90')
FEAST_OBJECTS   := $(addprefix $(FEAST_BUILD_DIR), $(notdir $(FEAST_SOURCES:.f90=.o)))
FEAST_LIB       := $(FEAST_BUILD_DIR)libfeast.a
FEAST_FFLAGS    := $(FFLAGS) -nogen-interfaces -warn nounused
LIBS            := $(LIBS) $(FEAST_LIB)

dirs: $(FEAST_BUILD_DIR)
$(FEAST_BUILD_DIR):
	@mkdir -p $(FEAST_BUILD_DIR)

$(FEAST_LIB): $(FEAST_OBJECTS)
	@printf "%b" "$(FC_COL)$(AR)$(NO_COL) rcs $(OU_COL)$(FEAST_LIB)$(NO_COL) $(IN_COL)$(FEAST_OBJECTS)$(NO_COL)\n\n"
	@$(AR) rcs $(FEAST_LIB) $(FEAST_OBJECTS)

$(FEAST_BUILD_DIR)%.o: $(FEAST_SRC_DIR)%.f90
	@printf "%b" "$(FC_COL)$(FC)$(NO_COL) $(FEAST_FFLAGS) -c $(IN_COL)$<$(NO_COL) -o $(OU_COL)$@$(NO_COL)\n\n"
	@$(FC) $(FEAST_FFLAGS) -c $< -o $@

feast: $(FEAST_LIB)

clean_feast:
	rm -f $(FEAST_BUILD_DIR)*.o $(FEAST_BUILD_DIR)*.a

clean_all: clean_feast

.PHONY: feast clean_feast