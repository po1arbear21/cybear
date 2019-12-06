#include "../assert.f90.inc"

module sqrtm_m
  use matrix_m
  implicit none

  interface sqrtm
  ! module procedure :: sqrtm_dense_cmplx
  ! module procedure :: sqrtm_dense_real
  module procedure :: sqrtm_hessenberg_cmplx
  module procedure :: sqrtm_hessenberg_real
  end interface

contains

#define T real
#define TT real
#include "sqrtm_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "sqrtm_imp.f90.inc"

end module
