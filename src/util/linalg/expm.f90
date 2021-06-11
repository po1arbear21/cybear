#include "../macro.f90.inc"

module expm_m

  use error_m,  only: assert_failed
  use lapack95
  use matrix_m, only: dense_real, matrix_add

  implicit none

  private
  public :: expm

  interface expm
    module procedure :: expm_dense_real
  end interface

contains

#define T  real
#define TT real
#define M  dense
#include "expm_imp.f90.inc"

end module
