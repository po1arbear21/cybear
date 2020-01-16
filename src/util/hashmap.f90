module hashmap_m
  use util_m
  use vector_m
  implicit none

#define T int
#define TT integer
#include "hashmap_def.f90.inc"

contains

#define T int
#define TT integer
#include "hashmap_imp.f90.inc"

end module
