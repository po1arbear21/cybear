#include "../macro.f90.inc"

module sqrtm_m

  use error_m,  only: assert_failed
  use matrix_m, only: dense_real, hessenberg_real, matrix_convert, triang_real
  use schur_m,  only: schur

  implicit none

  private
  public sqrtm

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
