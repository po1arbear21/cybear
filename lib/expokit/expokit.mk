EXPOKIT_SRC_DIR   := $(dir $(lastword $(MAKEFILE_LIST)))
EXPOKIT_BUILD_DIR := $(BUILD_DIR)expokit/
EXPOKIT_SOURCES   := $(shell find $(EXPOKIT_SRC_DIR) -name '*.f')
EXPOKIT_OBJECTS   := $(addprefix $(EXPOKIT_BUILD_DIR), $(notdir $(EXPOKIT_SOURCES:.f=.o)))
EXPOKIT_LIB       := $(EXPOKIT_BUILD_DIR)libexpokit.a
EXPOKIT_FFLAGS    := $(FFLAGS) -nogen-interfaces -warn nounused
LIBS             := $(LIBS) $(EXPOKIT_LIB)

dirs: $(EXPOKIT_BUILD_DIR)
$(EXPOKIT_BUILD_DIR):
	@mkdir -p $(EXPOKIT_BUILD_DIR)

$(EXPOKIT_LIB): $(EXPOKIT_OBJECTS)
	@printf "%b" "$(FC_COL)$(AR)$(NO_COL) rcs $(OU_COL)$(EXPOKIT_LIB)$(NO_COL) $(IN_COL)$(EXPOKIT_OBJECTS)$(NO_COL)\n\n"
	@$(AR) rcs $(EXPOKIT_LIB) $(EXPOKIT_OBJECTS)

$(EXPOKIT_BUILD_DIR)%.o: $(EXPOKIT_SRC_DIR)%.f
	@printf "%b" "$(FC_COL)$(FC)$(NO_COL) $(EXPOKIT_FFLAGS) -c $(IN_COL)$<$(NO_COL) -o $(OU_COL)$@$(NO_COL)\n\n"
	@$(FC) $(EXPOKIT_FFLAGS) -c $< -o $@

expokit: $(EXPOKIT_LIB)

clean_expokit:
	rm -f $(EXPOKIT_BUILD_DIR)*.o $(EXPOKIT_BUILD_DIR)*.a

clean_all: clean_expokit

.PHONY: expokit clean_expokit
