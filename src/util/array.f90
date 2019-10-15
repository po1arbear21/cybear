module array_m
  implicit none

#define T int
#define TT integer
#include "array_def.f90.inc"

#define T real
#define TT real
#include "array_def.f90.inc"

#define T cmplx
#define TT complex
#include "array_def.f90.inc"

end module
