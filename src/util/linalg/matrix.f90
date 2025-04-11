m4_include(../macro.f90.inc)

module matrix_m

  use bin_search_m,     only: bin_search
  use blas95
  use error_m,          only: assert_failed, program_error
  use high_precision_m, only: hp_dot
  m4_ifdef({m4_ilupack},{
  use ilupack_m
  })
  use iso_fortran_env,  only: int32, int64
  m4_ifdef({m4_klu2},{
  use klu2_m,           only: create_klu2_handle_c, create_klu2_handle_r, destruct_klu2_handle_c, &
    &                         destruct_klu2_handle_r, klu2_factorize, klu2_solve
  })
  use lapack95
  m4_ifdef({m4_mumps},{
  use mumps_m,          only: create_mumps_handle_c, create_mumps_handle_r, destruct_mumps_handle_c, &
    &                         destruct_mumps_handle_r, mumps_factorize, mumps_solve
  })
  use lsqr_m,           only: lsqr_solver_ez
  use omp_lib,          only: omp_get_thread_num, omp_get_num_threads
  use pardiso_m,        only: create_pardiso_handle, destruct_pardiso_handle, pardiso_factorize, pardiso_solve
  use qsort_m,          only: qsort
  use sparse_idx_m,     only: sparse_idx
  m4_ifdef({m4_spike},{
  use spike_m, only: spike_params, spike_factorize, spike_solve
  })
  use util_m,           only: int2str
  use vector_m,         only: vector_cmplx, vector_int, vector_log, vector_real

  implicit none

  private
  public SPSOLVER_PARDISO
  m4_ifdef({m4_klu2},{
  public SPSOLVER_KLU2
  })
  m4_ifdef({m4_mumps},{
  public SPSOLVER_MUMPS
  })
  m4_ifdef({m4_ilupack},{
  public SPSOLVER_ILUPACK
  })
  public BSOLVER_LAPACK
  m4_ifdef({m4_spike},{
  public BSOLVER_SPIKE
  })
  public DSOLVER_LAPACK
  public SOLVER_GMRES
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
  public block_real
  public block_cmplx
  public block_ptr_real
  public block_ptr_cmplx

  public matrix_add
  public matrix_mul
  public matrix_convert
  public matrix_diag
  public matrix_approx

  ! sparse solvers
  integer, parameter :: SPSOLVER_PARDISO = 1
  m4_ifdef({m4_klu2},{
  integer, parameter :: SPSOLVER_KLU2    = 2
  })
  m4_ifdef({m4_mumps},{
  integer, parameter :: SPSOLVER_MUMPS   = 3
  })
  m4_ifdef({m4_ilupack},{
  integer, parameter :: SPSOLVER_ILUPACK = 4
  })

  ! band solvers
  integer, parameter :: BSOLVER_LAPACK = 5
  m4_ifdef({m4_spike},{
  integer, parameter :: BSOLVER_SPIKE = 6
  })

  ! dense solvers
  integer, parameter :: DSOLVER_LAPACK = 7

  ! iterative general solvers
  integer, parameter :: SOLVER_GMRES = 8

  ! type list
  m4_define({m4_list},{
    m4_X(real)
    m4_X(cmplx)
  })

  ! include matrix type definitions
  m4_define({m4_X},{
    m4_define({T},$1)
    m4_include(matrix_def.f90.inc)
    m4_include(band_def.f90.inc)
    m4_include(dense_def.f90.inc)
    m4_include(sparse_def.f90.inc)
    m4_include(block_def.f90.inc)
    m4_undefine({T})
  })
  m4_list

  ! matrix arithmetics and conversions
  m4_define({m4_X},{
    m4_define({T},$1)
    m4_include(matrix_arith_def.f90.inc)
    m4_include(matrix_conv_def.f90.inc)
    m4_undefine({T})
  })
  m4_list

contains

  ! include matrix type procedures (subtypes in submodules)
  m4_define({m4_X},{
    m4_define({T},$1)
    m4_include(matrix_imp.f90.inc)
    m4_undefine({T})
  })
  m4_list

end module
