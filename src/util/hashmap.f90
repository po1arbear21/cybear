#include "macro.f90.inc"

module hashmap_m

  use error_m
  use string_m, only: string
  use util_m,   only: hash
  use vector_m, only: vector_int, vector_log, vector_string, vector_real, vector_cmplx

  implicit none

  private
  public hashmap_cmplx
  public hashmap_int
  public hashmap_log
  public hashmap_real
  public hashmap_string

#define T int
#define TT integer
#include "hashmap_def.f90.inc"

#define T log
#define TT logical
#include "hashmap_def.f90.inc"

#define T string
#define TT type(string)
#include "hashmap_def.f90.inc"

#define T real
#define TT real
#include "hashmap_def.f90.inc"

#define T cmplx
#define TT complex
#include "hashmap_def.f90.inc"

contains

#define T int
#define TT integer
#include "hashmap_imp.f90.inc"

#define T log
#define TT logical
#include "hashmap_imp.f90.inc"

#define T string
#define TT type(string)
#include "hashmap_imp.f90.inc"

#define T real
#define TT real
#include "hashmap_imp.f90.inc"

#define T cmplx
#define TT complex
#include "hashmap_imp.f90.inc"

end module
