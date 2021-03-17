#include "../macro.f90.inc"

submodule (matrix_m) band_m

contains

#define T real
#define TT real
#include "band_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "band_imp.f90.inc"

end submodule
