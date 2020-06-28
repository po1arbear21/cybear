ILUPACK_BUILD_DIR 		:= $(BUILD_DIR)ilupack/
ILUPACK_LIB_NAME  		:= libilupack.a
ILUPACK_LIB       		:= $(ILUPACK_BUILD_DIR)$(ILUPACK_LIB_NAME)
ILUPACK_CODE_DIR_REL  := lib/ilupack/
ILUPACK_CODE_DIR_ABS  := $(shell pwd)/lib/ilupack/
ILUPACK_ARCHIVES 			:= Dilupack Zilupack amd blaslike metis sparspak
ILUPACK_ARCHIVES_REL  := $(addprefix $(ILUPACK_CODE_DIR_REL), $(foreach lname, $(ILUPACK_ARCHIVES), lib$(lname).a))
ILUPACK_ARCHIVES_ABS  := $(addprefix $(ILUPACK_CODE_DIR_ABS), $(foreach lname, $(ILUPACK_ARCHIVES), lib$(lname).a))
LIBS              		:= $(LIBS) $(ILUPACK_LIB)

# reason for cd: ar only works in current dir...
$(ILUPACK_LIB): $(ILUPACK_ARCHIVES_REL) | $(ILUPACK_BUILD_DIR)
	cd $(ILUPACK_BUILD_DIR) 																	   \
	$(foreach fname, $(ILUPACK_ARCHIVES_ABS), && ar -x $(fname) && ar qc $(ILUPACK_LIB_NAME) *.o && rm -f *.o)

$(ILUPACK_BUILD_DIR):
	mkdir -p $@

ilupack: $(ILUPACK_LIB)

clean_ilupack:
	rm -f $(ILUPACK_BUILD_DIR)*.{a,o}

clean_all: clean_ilupack

.PHONY: ilupack clean_ilupack


#
#  DOCUMENTATION
#
# ILUPACK_ARCHIVES:
#	 	- ordering is important: ilupacks before other libs!
#
# why all the "_REL", "_ABS"?
# 	- "ar" can only extract in the current directory. thus, we need to work with absolute file paths.
# 	- relative path is needed as a make target. check if libs have been copied into repo folder.
#
# what is the "|" in the prerequisite?
#		it is an "order-only prereq". advantage:
# 		if file is created in directory
# 		-> directory time stamp is updated
# 		-> usually this will force rebuilding but not for "order-only prereq".
