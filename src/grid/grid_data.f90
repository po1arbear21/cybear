#include "../util/macro.f90.inc"

submodule (grid_m) grid_data_m

contains

#define T int
#define TT integer
#include "grid_data_imp.f90.inc"

#define T log
#define TT logical
#define TLOG
#include "grid_data_imp.f90.inc"

#define T real
#define TT real
#include "grid_data_imp.f90.inc"

#define T cmplx
#define TT complex
#include "grid_data_imp.f90.inc"

end submodule
