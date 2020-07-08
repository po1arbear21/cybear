#===============================================================================
# Copyright 2005-2020 Intel Corporation.
#
# This software and the related documents are Intel copyrighted  materials,  and
# your use of  them is  governed by the  express license  under which  they were
# provided to you (License).  Unless the License provides otherwise, you may not
# use, modify, copy, publish, distribute,  disclose or transmit this software or
# the related documents without Intel's prior written permission.
#
# This software and the related documents  are provided as  is,  with no express
# or implied  warranties,  other  than those  that are  expressly stated  in the
# License.
#===============================================================================

#  Content:
#      Intel(R) Math Kernel Library LAPACK95 library creation
#******************************************************************************

help:
	@echo "Usage: make {libia32|libintel64} INSTALL_DIR=<path>"
	@echo "            [interface=interface_name] [FC=compiler_name]"
	@echo
	@echo "Intel(R) Fortran Compiler as default"
	@echo
	@echo "INSTALL_DIR    - library and the .mod files will be built and installed"
	@echo "                 in subdirectories of <path> the same way as if <path>"
	@echo "                 were the Intel(R) MKL installation directory."
	@echo "interface_name - can be lp64 or ilp64 for intel64. Default value is lp64."
	@echo "FC             - can be ifort, gfortran or pgf95. Default value is ifort."

##------------------------------------------------------------------------------
## examples of using:
##
## make libia32 INSTALL_DIR=lib95    - build  by Intel(R) Fortran Compiler  (as default)
##                                   for 32-bit applications, static linking
##
## make libintel64 INSTALL_DIR=lib95 - build  by Intel(R) Fortran Compiler  (as default)
##                                   for Intel(R) 64 applications, static linking
##------------------------------------------------------------------------------

include lapack95.lst

iface_95 = lapack.f90
mod95 = lapack95.mod
iface_77 = lapack_interfaces.f90

ifndef FC
  FC=ifort
else
  ifeq ($(FC),f77)
    FC=ifort
  endif
endif

opts0    =
optsc    = -c $(opts0)
sufobj   = o
suflib   = a
sufmod   = mod
src_path = source
LIBR     = ar -rs

ifndef MKLROOT
MKLROOT = ../..
endif

obj_path = $(INSTALL_DIR)/lib/$(IA)/obj_lapack95_$(IA)_$(interface)
obj77_path = $(obj_path)/obj77

inc_path = $(MKLROOT)/include

MODULE_DIR = $(INSTALL_DIR)/include

objs    = $(src_lapack95:.f90=.$(sufobj))
objects = $(addprefix $(obj_path)/,$(objs))
mod77 = $(obj77_path)/$(iface_77:.f90=.obj)

FOPTS =

ifeq ($(FC),pgf95)
    auto=-mp
    # space after module
    opt_mod = -module #
else ifeq ($(FC),gfortran)
    auto=-fautomatic
    opt_mod = -J
else
    auto=-auto
    # space after module
    opt_mod = -module #
endif

FOPTS += $(auto)

ifeq ($(PLAT),lnx32)
    IA = ia32
    mod_path = $(MODULE_DIR)/ia32
    mkl_lapack95_lib      = libmkl_lapack95.$(suflib)
    FOPTS += $(if $(filter gfortran, $(FC)), -m32)
else
    IA = intel64
    ifeq ($(interface),ilp64)
        FOPTS += $(if $(filter gfortran, $(FC)), -fdefault-integer-8, -i8)
        mod_path = $(MODULE_DIR)/intel64/ilp64
        mkl_lapack95_lib      = libmkl_lapack95_ilp64.$(suflib)
    else
        mod_path = $(MODULE_DIR)/intel64/lp64
        mkl_lapack95_lib      = libmkl_lapack95_lp64.$(suflib)
    endif
endif

check:
ifndef PLAT
	$(error Undefined platform: PLAT=$(PLAT). Set one of: lnx32 lnx32e)
endif

vpath %.f90       $(src_path)
vpath %.f90       $(inc_path)
vpath %.$(sufobj) $(obj_path)

$(obj_path)/%.$(sufobj): %.f90
	$(FC) $(SPEC_OPT) $(FOPTS) $(opt_mod)$(obj77_path) $(optsc) -o $(obj_path)/$(@F) $<

$(mod_path)/$(mod95): $(inc_path)/$(iface_95)
	mkdir -p $(obj77_path)
	mkdir -p $(mod_path)
	$(FC) $(SPEC_OPT) $(FOPTS) $(opt_mod)$(mod_path) $(optsc) -o $(obj77_path)/mkl_lapack.$(sufobj) $<

.PHONY: build lib mod clean check checkdir lib32 libia32 libem64t libintel64 objects clean_obj cleania32 cleanintel64

libia32 lib32: checkdir
	make PLAT=lnx32 build MKLROOT=$(MKLROOT) INSTALL_DIR=$(INSTALL_DIR) FC=$(FC)
libintel64 libem64t: checkdir
	make PLAT=lnx32e build interface=$(interface) MKLROOT=$(MKLROOT) INSTALL_DIR=$(INSTALL_DIR) FC=$(FC)

build:	check lib mod clean_obj
lib:	$(INSTALL_DIR)/lib/$(IA)/$(mkl_lapack95_lib)
mod:	$(mod_path)/$(mod95)

$(INSTALL_DIR)/lib/$(IA)/$(mkl_lapack95_lib):
	mkdir -p $(obj77_path)
	$(FC) $(SPEC_OPT) $(FOPTS) $(opt_mod)$(obj77_path) $(optsc) -o $(mod77) $(src_path)/$(iface_77)
	make objects
	cd $(obj_path) ; $(LIBR) $(mkl_lapack95_lib) *.$(sufobj)
	mv $(obj_path)/$(mkl_lapack95_lib) $(INSTALL_DIR)/lib/$(IA)

objects:$(objects)

checkdir:
ifndef INSTALL_DIR
	$(error Not given INSTALL_DIR)
endif

clean: cleania32 cleanintel64
cleania32: checkdir
	-rm -rf $(INSTALL_DIR)/lib/ia32/obj_lapack95_32_*
	-rm -f $(INSTALL_DIR)/lib/ia32/libmkl_lapack95.$(suflib)
cleanintel64: checkdir
	-rm -rf $(INSTALL_DIR)/lib/intel64/obj_lapack95_intel64_*
	-rm -f $(INSTALL_DIR)/lib/intel64/libmkl_lapack95_*lp64.$(suflib)

clean_obj: check
	-rm -rf $(obj_path)
