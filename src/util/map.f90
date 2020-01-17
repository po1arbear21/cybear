module map_m
  use vector_m
  implicit none

#define T string
#define TT type(string)
#define U int
#define UU integer
#include "map_def.f90.inc"

#define T string
#define TT type(string)
#define U string
#define UU type(string)
#include "map_def.f90.inc"

#define T string
#define TT type(string)
#define U real
#define UU real
#include "map_def.f90.inc"

#define T string
#define TT type(string)
#define U cmplx
#define UU complex
#include "map_def.f90.inc"

contains

#define T string
#define TT type(string)
#define U int
#define UU integer
#include "map_imp.f90.inc"

#define T string
#define TT type(string)
#define U string
#define UU type(string)
#include "map_imp.f90.inc"

#define T string
#define TT type(string)
#define U real
#define UU real
#include "map_imp.f90.inc"

#define T string
#define TT type(string)
#define U cmplx
#define UU complex
#include "map_imp.f90.inc"

end module
