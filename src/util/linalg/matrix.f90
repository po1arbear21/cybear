#include "../macro.f90.inc"

module matrix_m

  use bin_search_m,     only: bin_search
  use blas95
  use error_m,          only: assert_failed, program_error
  use high_precision_m, only: hp_dot
#ifdef USE_ILUPACK
  use ilupack_m
#endif
  use iso_fortran_env,  only: int32, int64
  use lapack95
#ifdef USE_MUMPS
  use mumps_m,          only: create_mumps_handle_c, create_mumps_handle_r, destruct_mumps_handle_c, &
    &                         destruct_mumps_handle_r, mumps_factorize, mumps_solve
#endif
  use omp_lib,          only: omp_get_thread_num, omp_get_num_threads
  use pardiso_m,        only: create_pardiso_handle, destruct_pardiso_handle, pardiso_factorize, pardiso_solve
  use qsort_m,          only: qsort
  use sparse_idx_m,     only: sparse_idx
  use util_m,           only: int2str
  use vector_m,         only: vector_cmplx, vector_int, vector_log, vector_real

  implicit none

  private
  public SOLVER_PARDISO
#ifdef USE_MUMPS
  public SOLVER_MUMPS
#endif
#ifdef USE_ILUPACK
  public SOLVER_ILUPACK
#endif
  public default_solver
  public matrix_real
  public matrix_cmplx
  public matrix_ptr_real
  public matrix_ptr_cmplx
  public dense_real
  public dense_cmplx
  public dense_ptr_real
  public dense_ptr_cmplx
  public dense_eye_real
  public dense_eye_cmplx
  public sparse_real
  public sparse_cmplx
  public sparse_ptr_real
  public sparse_ptr_cmplx
  public sparse_eye_real
  public sparse_eye_cmplx
  public sparse_zero_real
  public sparse_zero_cmplx
  public spbuild_real
  public spbuild_cmplx
  public band_real
  public band_cmplx
  public band_ptr_real
  public band_ptr_cmplx
  public band_eye_real
  public band_eye_cmplx
  public hessenberg_real
  public hessenberg_cmplx
  public hessenberg_ptr_real
  public hessenberg_ptr_cmplx
  public triang_real
  public triang_cmplx
  public triang_ptr_real
  public triang_ptr_cmplx
  public block_real
  public block_cmplx
  public block_ptr_real
  public block_ptr_cmplx

  public matrix_add
  public matrix_convert
  public matrix_diag
  public matrix_approx

  ! sparse solvers
  integer, parameter :: SOLVER_PARDISO = 1
#ifdef USE_MUMPS
  integer, parameter :: SOLVER_MUMPS   = 2
#endif
#ifdef USE_ILUPACK
  integer, parameter :: SOLVER_ILUPACK = 3
#endif
  integer            :: default_solver = SOLVER_PARDISO

  ! matrix types
#define T real
#define TT real
#include "matrix_def.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "matrix_def.f90.inc"

#define T real
#define TT real
#include "band_def.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "band_def.f90.inc"

#define T real
#define TT real
#include "dense_def.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "dense_def.f90.inc"

#define T real
#define TT real
#include "hessenberg_def.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "hessenberg_def.f90.inc"

#define T real
#define TT real
#include "sparse_def.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "sparse_def.f90.inc"

#define T real
#define TT real
#include "triang_def.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "triang_def.f90.inc"

#define T real
#define TT real
#include "block_def.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "block_def.f90.inc"

  ! matrix arithmetics
#define T real
#define TT real
#include "matrix_arith_def.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "matrix_arith_def.f90.inc"

  ! matrix conversions
#define T real
#define TT real
#include "matrix_conv_def.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "matrix_conv_def.f90.inc"

contains

#define T real
#define TT real
#include "matrix_imp.f90.inc"

#define T cmplx
#define TT complex
#define TCMPLX
#include "matrix_imp.f90.inc"

end module
