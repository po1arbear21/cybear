BLAS95_SRC_DIR   := $(dir $(lastword $(MAKEFILE_LIST)))
BLAS95_BUILD_DIR := $(F95_BUILD_DIR)
BLAS95_LIB       := $(BLAS95_BUILD_DIR)/blas95.a
LIBS             := $(LIBS) $(BLAS95_LIB)

$(BLAS95_LIB):
	cd $(BLAS95_SRC_DIR) && make libintel64 INSTALL_DIR=../../$(BLAS95_BUILD_DIR) interface=ilp64 FC=gfortran
	mv $(BLAS95_BUILD_DIR)/lib/intel64/libmkl_blas95_ilp64.a $(BLAS95_LIB) && rm -r $(BLAS95_BUILD_DIR)/lib
	mv $(BLAS95_BUILD_DIR)/include/intel64/ilp64/*.mod $(BLAS95_BUILD_DIR)/ && rm -r $(BLAS95_BUILD_DIR)/include

blas95: $(BLAS95_LIB)

clean_blas95:
	rm -f $(BLAS95_BUILD_DIR)/*.a $(BLAS95_BUILD_DIR)/*.mod

clean_all: clean_blas95

.PHONY: blas95 clean_blas95
