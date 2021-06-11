#include "macro.f90.inc"

module qsort_m

  use error_m,  only: assert_failed
  use string_m, only: string, operator(<), operator(<=), operator(>)

  implicit none

  private
  public qsort

  interface qsort
    module procedure :: qsort_int, qsort_string, qsort_real
  end interface

contains

#define T int
#define TT integer
#include "qsort_imp.f90.inc"

#define T string
#define TT type(string)
#include "qsort_imp.f90.inc"

#define T real
#define TT real
#include "qsort_imp.f90.inc"

end module
