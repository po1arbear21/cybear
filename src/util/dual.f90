module dual_m
  implicit none

  private
  public dual_1, dual_2, dual_3, dual_4, dual_5, dual_6, dual_7, dual_8
  public operator(+), operator(-), operator(*), operator(/), operator(**), operator(.dot.)
  public abs, sqrt, exp, log, sin, cos, tan

#define N 1
#include "dual_def.f90.inc"
#define N 2
#include "dual_def.f90.inc"
#define N 3
#include "dual_def.f90.inc"
#define N 4
#include "dual_def.f90.inc"
#define N 5
#include "dual_def.f90.inc"
#define N 6
#include "dual_def.f90.inc"
#define N 7
#include "dual_def.f90.inc"
#define N 8
#include "dual_def.f90.inc"

contains

#define N 1
#include "dual_imp.f90.inc"
#define N 2
#include "dual_imp.f90.inc"
#define N 3
#include "dual_imp.f90.inc"
#define N 4
#include "dual_imp.f90.inc"
#define N 5
#include "dual_imp.f90.inc"
#define N 6
#include "dual_imp.f90.inc"
#define N 7
#include "dual_imp.f90.inc"
#define N 8
#include "dual_imp.f90.inc"

end module
