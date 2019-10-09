module vector_m
  implicit none

#define T int
#define TT integer
#include "src/util/vector_def.f90.inc"

#define T log
#define TT logical
#include "src/util/vector_def.f90.inc"

#define T char
#define TT character(len=64)
#include "src/util/vector_def.f90.inc"

#define T real
#define TT real
#include "src/util/vector_def.f90.inc"

#define T cmplx
#define TT complex
#include "src/util/vector_def.f90.inc"

contains

#define T int
#define TT integer
#include "src/util/vector_imp.f90.inc"

#define T log
#define TT logical
#include "src/util/vector_imp.f90.inc"

#define T char
#define TT character(len=64)
#include "src/util/vector_imp.f90.inc"

#define T real
#define TT real
#include "src/util/vector_imp.f90.inc"

#define T cmplx
#define TT complex
#include "src/util/vector_imp.f90.inc"

end module
