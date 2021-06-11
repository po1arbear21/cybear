module array_m

  use string_m, only: string

  implicit none

  private
  ! public (see array_def.f90.inc)

#define T int
#define TT integer
#include "array_def.f90.inc"

#define T log
#define TT logical
#include "array_def.f90.inc"

#define T string
#define TT type(string)
#include "array_def.f90.inc"

#define T real
#define TT real
#include "array_def.f90.inc"

#define T cmplx
#define TT complex
#include "array_def.f90.inc"

end module
