#include "../macro.f90.inc"

module matrix_m
  use array_m
  use blas95
  use error_m
  use high_precision_m
  use lapack95
  use omp_lib
  use pardiso_m
  use qsort_m
  use vector_m
  implicit none

#define T real
#define TT real
#include "matrix_def.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "matrix_def.f90.inc"

contains

#define T real
#define TT real
#include "matrix_imp.f90.inc"

#define T real
#define TT real
#include "dense_imp.f90.inc"

#define T real
#define TT real
#include "sparse_imp.f90.inc"

#define T real
#define TT real
#include "band_imp.f90.inc"

#define T real
#define TT real
#include "hessenberg_imp.f90.inc"

#define T real
#define TT real
#include "triang_imp.f90.inc"

#define T real
#define TT real
#include "block_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "matrix_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "dense_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "sparse_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "band_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "hessenberg_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "triang_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "block_imp.f90.inc"

end module
