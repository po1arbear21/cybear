#include "macro.f90.inc"

module bin_search_m

  use error_m,         only: program_error, assert_failed
  use iso_fortran_env, only: int32, int64

  implicit none

  private
  public bin_search
  public BS_NEAR, BS_LESS, BS_GREAT

  interface bin_search
    module procedure :: bin_search_int32, bin_search_int64, bin_search_real
  end interface

  ! binary search modes
  integer, parameter :: BS_NEAR  = 1
  integer, parameter :: BS_LESS  = 2
    !! returned index yields smaller/equal value wrt query point, i.e. x(i0) <= xq
    !! (as long as xq is within range)
  integer, parameter :: BS_GREAT = 3
    !! returned index yields larger/equal value wrt query point, i.e. x(i0) >= xq
    !! (as long as xq is within range)

contains

#define T int32
#define TT integer(int32)
#include "bin_search_imp.f90.inc"

#define T int64
#define TT integer(int64)
#include "bin_search_imp.f90.inc"

#define T real
#define TT real
#include "bin_search_imp.f90.inc"

end module
