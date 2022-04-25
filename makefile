# set default values for BUILD, COMPILER and library USE flags (overwritten by command line)
include options.mk

# set flags for M4 macro processor
include m4.mk

# set compiler flags
ifeq ($(COMPILER),intel)
	include intel.mk
endif
ifeq ($(COMPILER),gnu)
	include gnu.mk
endif

# libraries
include lib.mk

# colors
include color.mk

# main target
all: dirs depend

# directories
BUILD_DIR := build
TMP_DIR := $(BUILD_DIR)/tmp
dirs: $(BUILD_DIR) $(TMP_DIR)
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)
$(TMP_DIR):
	@mkdir -p $(TMP_DIR)

# find all fortran and c source files in src directory
F90_SOURCES := $(shell find src -name '*.f90')
C_SOURCES := $(shell find src -name '*.c')

# create list of c objects
vpath %.c $(sort $(dir $(C_SOURCES)))
C_OBJECTS := $(addprefix $(BUILD_DIR)/, $(addsuffix .o, $(notdir $(C_SOURCES))))

# generate targets + dependencies (one target per program found in sources)
depend: $(BUILD_DIR)/depend.mk depend.py
$(BUILD_DIR)/depend.mk: $(BUILD_DIR) $(F90_SOURCES)
	@printf "%b" "$(PYTHON_P) depend.py $(INP_COLOR)$(F90_SOURCES)$(OFF_COLOR) > $(OUT_COLOR)$(BUILD_DIR)/depend.mk$(OFF_COLOR) || $(RM_P) -f $(OUT_COLOR)$(BUILD_DIR)/depend.mk$(OFF_COLOR)\n\n"
	@$(PYTHON) depend.py $(F90_SOURCES) > $(BUILD_DIR)/depend.mk || rm -f $(BUILD_DIR)/depend.mk
include $(BUILD_DIR)/depend.mk

# build all targets
all: $(TARGETS)

# rule for m4 files
$(BUILD_DIR)/%.f90:
	@printf "%b" "$(M4_P) $(M4FLAGS) $(INP_COLOR)$<$(OFF_COLOR) > $(OUT_COLOR)$@$(OFF_COLOR)\n\n"
	@$(M4) $(M4FLAGS) -I$(dir $<) $< > $@

# rule for anchor files
$(BUILD_DIR)/%.anc:
	@printf "%b" "$(FC_P) $(FFLAGS) $(FINCLUDE) $(FMODULE) $(BUILD_DIR) $(FSYNTAXONLY) -c $(INP_COLOR)$<$(OFF_COLOR)\n\n"
	@$(FC) $(FFLAGS) $(FINCLUDE) $(FMODULE) $(BUILD_DIR) $(FSYNTAXONLY) -c $<
	@mv $(notdir $(<:.f90=.i90)) $(BUILD_DIR) 2>/dev/null || true
	@touch $@

# rule for c object files
$(BUILD_DIR)/%.c.o: %.c
	@printf "%b" "$(CC_P) $(CFLAGS) -c $(INP_COLOR)$<$(OFF_COLOR) -o $(OUT_COLOR)$@$(OFF_COLOR)\n\n"
	@$(CC) $(CFLAGS) -c $< -o $@

# rule for object files
$(BUILD_DIR)/%.o:
	@printf "%b" "$(FC_P) $(FFLAGS) $(FINCLUDE) $(FMODULE) $(TMP_DIR) -c $(INP_COLOR)$<$(OFF_COLOR) -o $(OUT_COLOR)$@$(OFF_COLOR)\n\n"
	@$(FC) $(FFLAGS) $(FINCLUDE) $(FMODULE) $(TMP_DIR) -c $< -o $@

# delete files in build folder
clean:
	@printf "%b" "$(RM_P) -f $(INP_COLOR)$(TMP_DIR)/*.{s,}mod $(BUILD_DIR)/*.{f90,i90,anc,mod,smod,o,a} $(BUILD_DIR)/depend.mk $(TARGETS)$(OFF_COLOR)\n"
	@rm -f $(TMP_DIR)/*.{s,}mod $(BUILD_DIR)/*.{f90,i90,anc,mod,smod,o,a} $(BUILD_DIR)/depend.mk $(TARGETS)

# search and delete accidentally created mod files in source folder
clean_mod:
	@printf "%b" "$(FIND_P) $(INP_COLOR)src$(OFF_COLOR) -type f -name '*.mod' -delete\n"
	@find src -type f -name '*.mod' -delete

doc: all
	@printf "%b" "$(FORD_P) -d $(BUILD_DIR) README.md $(FINCLUDE) -I$(MKLROOT)/include\n"
	@ford -d $(BUILD_DIR) README.md $(FINCLUDE) -I$(MKLROOT)/include

.PHONY: all dirs doc clean clean_mod
