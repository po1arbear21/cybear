#include "../util/macro.f90.inc"

module ilupack_m
  !! wrapper around ILUPACK. uses interfaces for ILUPACK defined in ilupack_interf_m.

  use error_m
  use util_m

  implicit none

  private
  public :: ilupack_real
  public :: ilupack_cmplx

#define T0 D
#define T  real
#define TT real
#include "ilupack_def.f90.inc"

#define T0 Z
#define T  cmplx
#define TT complex
#include "ilupack_def.f90.inc"

  character(*), parameter :: ILUPACK_FACTOR_ERROR(-7:-1) = [  &
    & "buffers are too small                               ", &
    & "zero column encountered                             ", &
    & "zero row encountered                                ", &
    & "Illegal value for lfil                              ", &
    & "matrix U overflow, increase elbow and retry         ", &
    & "matrix L overflow, increase elbow and retry         ", &
    & "Error. input matrix may be wrong.                   "  ]

  character(*), parameter :: ILUPACK_SOLVER_ERROR(-3:-1) = [  &
    & "algorithm breaks down                               ", &
    & "not enough work space                               ", &
    & "too many iterations                                 "  ]


contains

#define T0 D
#define T  real
#define TT real
#include "ilupack_imp.f90.inc"

#define T0 Z
#define T  cmplx
#define TT complex
#include "ilupack_imp.f90.inc"

end module
