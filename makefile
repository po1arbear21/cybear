FC             = ifort
FFLAGS0        = -march=native -real-size 64 -i8 -fpp -module build/ -warn all -qopenmp
FFLAGS_DEBUG   = -O0 -mkl=sequential -g -check all -check noarg_temp_created -fpe1 -traceback -debug extended -D DEBUG
FFLAGS_RELEASE = -O2 -mkl
LIBS           = ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a

TARGET  = main
DIRS    = $(shell find src/ -type d)
SOURCES = $(shell find src/ -name '*.f90')
OBJECTS = $(addprefix build/, $(notdir $(SOURCES:.f90=.o)))

TEST_TARGET  = test_main
TEST_DIRS    = $(DIRS) $(shell find test/ -type d)
TEST_SOURCES = $(SOURCES) $(shell find test/ -name '*.f90')
TEST_OBJECTS = $(filter-out build/$(TARGET).o, $(addprefix build/, $(notdir $(TEST_SOURCES:.f90=.o))))

# convert spaces in $(DIRS) to colons for makedepf90
TEST_DIRS_COLON = $(subst $() $(),:,$(TEST_DIRS))

# set default test to all
USE_TEST := $(if $(USE_TEST),$(USE_TEST),all)

.PHONY: debug release test clean depend

debug: FFLAGS = $(FFLAGS0) $(FFLAGS_DEBUG)
debug: $(TARGET)

release: FFLAGS = $(FFLAGS0) $(FFLAGS_RELEASE)
release: $(TARGET)

test: FFLAGS = $(FFLAGS0) $(FFLAGS_DEBUG) -D USE_TEST=$(USE_TEST)
test: $(TEST_TARGET)

include .depend

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

$(TEST_TARGET): $(TEST_OBJECTS)
	$(FC) $(FFLAGS) -o $(TEST_TARGET) $(TEST_OBJECTS) $(LIBS)

build/%.o:
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f build/*.o build/*.mod $(TARGET) .depend

depend .depend:
	./makedepf90 $(TEST_SOURCES) -b build/ -I $(TEST_DIRS_COLON) > .depend
