#include "../macro.f90.inc"

submodule (matrix_m) sparse_m

contains

#define T real
#define TT real
#include "sparse_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "sparse_imp.f90.inc"

end submodule
