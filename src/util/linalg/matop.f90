#include "../macro.f90.inc"

module matop_m
  use matrix_m
  implicit none

#define T real
#define TT real
#include "matop_def.f90.inc"

#define T cmplx
#define TT complex
#include "matop_def.f90.inc"

contains

#define T real
#define TT real
#include "matop_imp.f90.inc"

#define T cmplx
#define TT complex
#include "matop_imp.f90.inc"

end module