#include "../macro.f90.inc"

module sqrtm_m
  use schur_m
  implicit none

  interface sqrtm
    module procedure :: sqrtm_dense_real
    module procedure :: sqrtm_hessenberg_real
  end interface

contains

#define T  real
#define TT real
#define M  dense
#include "sqrtm_imp.f90.inc"

#define T  real
#define TT real
#define M  hessenberg
#include "sqrtm_imp.f90.inc"

end module
