module bin_search_m
  use error_m
  use sparse_idx_m
  implicit none

  private
  public :: bin_search
  public :: BS_NEAR, BS_LESS, BS_GREAT

  interface bin_search
    module procedure :: bin_search_int, bin_search_idx, bin_search_real
  end interface

  ! binary search modes
  integer, parameter :: BS_NEAR  = 1
  integer, parameter :: BS_LESS  = 2
  integer, parameter :: BS_GREAT = 3

contains

#define T int
#define TT integer
#include "bin_search_imp.f90.inc"

#define T idx
#define TT integer(kind=SPARSE_IDX)
#include "bin_search_imp.f90.inc"

#define T real
#define TT real
#include "bin_search_imp.f90.inc"

end module