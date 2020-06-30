#include "../macro.f90.inc"

module matrix_m
  use array_m
  use bin_search_m
  use blas95
  use error_m
  use high_precision_m
  use lapack95
  use omp_lib
  use pardiso_m
  use qsort_m
  use vector_m
  implicit none

  private
  public :: matrix_real
  public :: matrix_cmplx
  public :: matrix_alc_real
  public :: matrix_alc_cmplx
  public :: matrix_ptr_real
  public :: matrix_ptr_cmplx
  public :: dense_real
  public :: dense_cmplx
  public :: dense_eye_real
  public :: dense_eye_cmplx
  public :: sparse_real
  public :: sparse_cmplx
  public :: sparse_eye_real
  public :: sparse_eye_cmplx
  public :: sparse_zero_real
  public :: sparse_zero_cmplx
  public :: spbuild_real
  public :: spbuild_cmplx
  public :: band_real
  public :: band_cmplx
  public :: band_eye_real
  public :: band_eye_cmplx
  public :: hessenberg_real
  public :: hessenberg_cmplx
  public :: triang_real
  public :: triang_cmplx
  public :: block_real
  public :: block_cmplx

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
