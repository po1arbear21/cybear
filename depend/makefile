FC            = ifort
DEPEND_FFLAGS = -real-size 64 -i8 -march=native -O2 -fpp -module depend/build/

.PHONY: clean_depend

depend/depend: depend/src/depend.f90 depend/build/src_file.o
	$(FC) $(DEPEND_FFLAGS) depend/src/depend.f90 -o depend/depend depend/build/src_file.o

depend/build/src_file.o: depend/src/src_file.f90 depend/src/vector_def.f90.inc depend/src/vector_imp.f90.inc
	$(FC) $(DEPEND_FFLAGS) -c depend/src/src_file.f90 -o depend/build/src_file.o

clean_depend:
	rm -f depend/build/*.{o,mod} depend/depend

clean_all: clean_depend
