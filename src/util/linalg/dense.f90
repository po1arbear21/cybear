#include "../macro.f90.inc"

submodule (matrix_m) dense_m

contains

#define T real
#define TT real
#include "dense_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "dense_imp.f90.inc"

end submodule
