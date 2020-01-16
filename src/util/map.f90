module map_m
  use vector_m
  implicit none

#define T int
#define TT integer
#define U int
#define UU integer
#include "map_def.f90.inc"

#define T int
#define TT integer
#define U logical
#define UU logical
#include "map_def.f90.inc"

#define T int
#define TT integer
#define U double
#define UU doubleprecision
#include "map_def.f90.inc"

#define T int
#define TT integer
#define U complex
#define UU doublecomplex
#include "map_def.f90.inc"

#define T char
#define TT character(len=64)
#define U int
#define UU integer
#include "map_def.f90.inc"

#define T char
#define TT character(len=64)
#define U logical
#define UU logical
#include "map_def.f90.inc"

#define T char
#define TT character(len=64)
#define U double
#define UU doubleprecision
#include "map_def.f90.inc"

#define T char
#define TT character(len=64)
#define U complex
#define UU doublecomplex
#include "map_def.f90.inc"

#define T int
#define TT integer
#define U vector_int
#define UU type(vector_int)
#include "map_def.f90.inc"

contains

#define T int
#define TT integer
#define U int
#define UU integer
#include "map_imp.f90.inc"

#define T int
#define TT integer
#define U logical
#define UU logical
#include "map_imp.f90.inc"

#define T int
#define TT integer
#define U double
#define UU doubleprecision
#include "map_imp.f90.inc"

#define T int
#define TT integer
#define U complex
#define UU doublecomplex
#include "map_imp.f90.inc"

#define T char
#define TT character(len=*)
#define U int
#define UU integer
#include "map_imp.f90.inc"

#define T char
#define TT character(len=*)
#define U logical
#define UU logical
#include "map_imp.f90.inc"

#define T char
#define TT character(len=*)
#define U double
#define UU doubleprecision
#include "map_imp.f90.inc"

#define T char
#define TT character(len=*)
#define U complex
#define UU doublecomplex
#include "map_imp.f90.inc"

#define T int
#define TT integer
#define U vector_int
#define UU type(vector_int)
#include "map_imp.f90.inc"

end module
