#include "../macro.f90.inc"

submodule (matrix_m) triang_m

contains

#define T real
#define TT real
#include "triang_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "triang_imp.f90.inc"

end submodule
