FC     = ifort
FFLAGS = -march=native -O2 -fpp -module build/

.PHONY: clean

depend: src/depend.f90 build/src_file.o
	$(FC) $(FFLAGS) src/depend.f90 -o depend build/src_file.o

build/src_file.o: src/src_file.f90 src/vector_def.f90.inc src/vector_imp.f90.inc
	$(FC) $(FFLAGS) -c src/src_file.f90 -o build/src_file.o

clean:
	rm -f build/*.o build/*.mod depend
