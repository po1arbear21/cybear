#include "../macro.f90.inc"

submodule (matrix_m) block_m

contains

#define T real
#define TT real
#include "block_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "block_imp.f90.inc"

end submodule
