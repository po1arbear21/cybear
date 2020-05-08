#include "../macro.f90.inc"

module matop_m
  use matrix_m
  implicit none

  private
  public :: matop_real
  public :: matop_cmplx
  public :: single_matop_real
  public :: single_matop_cmplx
  public :: chain_matop_real
  public :: chain_matop_cmplx

#define T real
#define TT real
#include "matop_def.f90.inc"

#define T cmplx
#define TT complex
#include "matop_def.f90.inc"

contains

#define T real
#define TT real
#include "matop_imp.f90.inc"

#define T cmplx
#define TT complex
#include "matop_imp.f90.inc"

end module