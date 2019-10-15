module vector_m
  use string_m
  implicit none

#define T int
#define TT integer
#include "vector_def.f90.inc"

#define T log
#define TT logical
#include "vector_def.f90.inc"

#define T string
#define TT type(string)
#include "vector_def.f90.inc"

#define T real
#define TT real
#include "vector_def.f90.inc"

#define T cmplx
#define TT complex
#include "vector_def.f90.inc"

contains

#define T int
#define TT integer
#include "vector_imp.f90.inc"

#define T log
#define TT logical
#include "vector_imp.f90.inc"

#define T string
#define TT type(string)
#include "vector_imp.f90.inc"

#define T real
#define TT real
#include "vector_imp.f90.inc"

#define T cmplx
#define TT complex
#include "vector_imp.f90.inc"

end module
