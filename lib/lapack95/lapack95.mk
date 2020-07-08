LAPACK95_SRC_DIR   := $(dir $(lastword $(MAKEFILE_LIST)))
LAPACK95_BUILD_DIR := build/${COMPILER}/F95
LAPACK95_LIB       := $(LAPACK95_BUILD_DIR)/lapack95.a
LIBS               := $(LIBS) $(LAPACK95_LIB)

$(LAPACK95_BUILD_DIR):
	@mkdir -p $(LAPACK95_BUILD_DIR)

$(LAPACK95_LIB):
	cd $(LAPACK95_SRC_DIR) && make libintel64 INSTALL_DIR=../../$(LAPACK95_BUILD_DIR) interface=ilp64 FC=gfortran
	mv $(LAPACK95_BUILD_DIR)/lib/intel64/libmkl_lapack95_ilp64.a $(LAPACK95_LIB) && rm -r $(LAPACK95_BUILD_DIR)/lib
	mv $(LAPACK95_BUILD_DIR)/include/intel64/ilp64/*.mod $(LAPACK95_BUILD_DIR)/ && rm -r $(LAPACK95_BUILD_DIR)/include

lapack95: $(LAPACK95_LIB)

clean_lapack95:
	rm -f $(LAPACK95_BUILD_DIR)/*.a $(LAPACK95_BUILD_DIR)/*.mod

clean_all: clean_lapack95

.PHONY: lapack95 clean_lapack95