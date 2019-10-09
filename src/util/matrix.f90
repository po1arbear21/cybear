#include "src/util/assert.f90.inc"
#include "src/util/macro.f90.inc"

module matrix_m
  use blas95
  use lapack95
  use mkl_spblas
  use pardiso_m
  use qsort_m
  use vector_m
  implicit none

#define T real
#define TT real
#include "src/util/matrix_def.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "src/util/matrix_def.f90.inc"

contains

#define T real
#define TT real
#include "src/util/matrix_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "src/util/matrix_imp.f90.inc"

end module
