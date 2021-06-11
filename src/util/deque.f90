module deque_m

  use string_m, only: string

  implicit none

  private
  public deque_int
  public deque_log
  public deque_string
  public deque_real
  public deque_cmplx

#define T int
#define TT integer
#include "deque_def.f90.inc"

#define T log
#define TT logical
#include "deque_def.f90.inc"

#define T string
#define TT type(string)
#include "deque_def.f90.inc"

#define T real
#define TT real
#include "deque_def.f90.inc"

#define T cmplx
#define TT complex
#include "deque_def.f90.inc"

contains

#define T int
#define TT integer
#include "deque_imp.f90.inc"

#define T log
#define TT logical
#include "deque_imp.f90.inc"

#define T string
#define TT type(string)
#include "deque_imp.f90.inc"

#define T real
#define TT real
#include "deque_imp.f90.inc"

#define T cmplx
#define TT complex
#include "deque_imp.f90.inc"

end module
