#include "../macro.f90.inc"

module schur_m
  use matrix_m
  implicit none

  interface schur
    module procedure :: schur_dense_real
    module procedure :: schur_dense_cmplx
    module procedure :: schur_hessenberg_real
    module procedure :: schur_hessenberg_cmplx
  end interface

contains

#define T real
#define TT real
#include "schur_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "schur_imp.f90.inc"

end module
