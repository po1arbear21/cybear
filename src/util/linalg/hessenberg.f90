#include "../macro.f90.inc"

submodule (matrix_m) hessenberg_m

contains

#define T real
#define TT real
#include "hessenberg_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "hessenberg_imp.f90.inc"

end submodule
