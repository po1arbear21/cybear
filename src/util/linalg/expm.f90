#include "../macro.f90.inc"

module expm_m
  use matrix_m
  implicit none

  interface expm
    module procedure :: expm_dense_real
  end interface

contains

#define T  real
#define TT real
#define M  dense
#include "expm_imp.f90.inc"

end module
