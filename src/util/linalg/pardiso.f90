m4_include(../macro.f90.inc)

module pardiso_m

  use deque_m,      only: deque_int
  use error_m,      only: assert_failed, program_error
  use sparse_idx_m, only: sparse_idx

  implicit none

  ! Privat vor Staat (Olaf Scholz does not approve of this message)

  private
  public create_pardiso_handle
  public destruct_pardiso_handle
  public pardiso_factorize
  public pardiso_solve

  ! include pardiso or pardiso_64 interface
  include "mkl_pardiso.fi"
  m4_define({PARDISO_ROUTINE},{m4_ifelse(m4_intsize,m4_idxsize,{pardiso},{pardiso_64})})

  type pardiso_handle
    type(MKL_PARDISO_HANDLE) :: pt(64)
      !! internal PARDISO handle
    integer(SPARSE_IDX)      :: maxfct
      !! number of matrices with identical sparse structure
    integer(SPARSE_IDX)      :: mnum
      !! which matrix should be used in solution phase
    integer(SPARSE_IDX)      :: mtype
      !! matrix type
    integer(SPARSE_IDX)      :: nrows
      !! number of equations
    integer(SPARSE_IDX)      :: msglvl
      !! message level
    integer(SPARSE_IDX)      :: iparm(64)
      !! PARDISO parameters
    logical                  :: factorized = .false.
      !! factorization complete flag
  end type

  integer, parameter   :: PARDISO_NUM_HANDLES = 128
  type(pardiso_handle) :: pardiso_handles(PARDISO_NUM_HANDLES)
  type(deque_int)      :: pardiso_free_handles

  interface pardiso_factorize
    module procedure :: pardiso_factorize_r
    module procedure :: pardiso_factorize_c
  end interface

  interface pardiso_solve
    module procedure :: pardiso_solve_r
    module procedure :: pardiso_solve_c
  end interface

  character(*), parameter :: PARDISO_ERROR(-12:0) = [                      &
    "pardiso_64 called from 32-bit library                              ", &
    "read/write error with OOC files                                    ", &
    "error opening OOC files                                            ", &
    "not enough memory for OOC                                          ", &
    "32-bit integer overflow problem                                    ", &
    "diagonal matrix is singular                                        ", &
    "reordering failed                                                  ", &
    "unclassified (internal) error                                      ", &
    "zero pivot, numerical factorization or iterative refinement problem", &
    "reordering problem                                                 ", &
    "not enough memory                                                  ", &
    "input inconsistent                                                 ", &
    "no error                                                           "  &
  ]

contains

  function create_pardiso_handle(nrows, cmplx) result(h)
    integer, intent(in) :: nrows
      !! number of rows
    logical, intent(in) :: cmplx
      !! complex/real flag
    integer             :: h
      !! return pardiso handle (index)

    integer :: i, mtype, iparm(64)

    ! get free handle
    !$omp critical (omp_pardiso_handles)
    if (.not. allocated(pardiso_free_handles%d)) then
      call pardiso_free_handles%init(PARDISO_NUM_HANDLES, x = [(i, i=1, PARDISO_NUM_HANDLES)])
    end if
    if (pardiso_free_handles%n < 1) call program_error("No free pardiso handles!")
    h = pardiso_free_handles%front()
    call pardiso_free_handles%pop_front()
    !$omp end critical (omp_pardiso_handles)

    associate (p => pardiso_handles(h))
      ! maximum number of factors with idential sparse structure
      p%maxfct = 1

      ! which matrix should be used in solution phase (1 <= mnum <= maxfct)
      p%mnum = 1

      ! matrix type: real unsymmetric or complex unsymmetric
      if (cmplx) then
        p%mtype = 13
      else
        p%mtype = 11
      end if

      ! number of rows
      p%nrows = nrows

      ! do not print statistical information
      p%msglvl = 0

      ! init pt, iparm
      mtype = int(p%mtype, kind = kind(mtype))
      call pardisoinit(p%pt, mtype, iparm)
      p%iparm = int(iparm, kind = SPARSE_IDX)

      ! ======================
      !     ipar
      ! ======================

      ! p%iparm( 1) = 1  ! Use default values
      ! p%iparm( 2) = 2  ! Fill-in reducing ordering for the input matrix (default: METIS)
      ! p%iparm( 4) = 0  ! Preconditioned CGS/CG
      ! p%iparm( 5) = 0  ! User permutation
      ! p%iparm( 6) = 0  ! Write solution on x (default: right-hand side b unchanged)
      ! p%iparm( 7)      ! output: Actual number of iterative refinement steps performed
        p%iparm( 8) = 9  ! Maximal number of iterative refinement steps
      ! p%iparm(10)      ! Pivoting permutation (default: 1 or 3 depending on cmplx)
      ! p%iparm(11) = 1  ! Scaling (default: on for unsymmetric)
      ! p%iparm(12) = 0  ! Solve with transposed or conjugate transposed matrix (default: off)
      ! p%iparm(13) = 1  ! Improved accuracy using non-symmetric weighted matching (default: on for unsymmetric)
      ! p%iparm(14)      ! output: Number of perturbed pivots
      ! p%iparm(15)      ! output: Peak memory on symbolic factorization
      ! p%iparm(16)      ! output: Permanent memory on symbolic factorization
      ! p%iparm(17)      ! output: Size of factors/Peak memory on numerical factorization and solution
      ! p%iparm(18) = -1 ! Report number of non-zero elements in the factors (default: on)
      ! p%iparm(19) = 0  ! Report number of floating point operations for factorization (default: off)
      ! p%iparm(20)      ! output: Report CG/CGS diagnostics
      ! p%iparm(21) = 1  ! Pivoting for symmetric indefinite matrices (default: Bunch-Kaufman)
      ! p%iparm(22)      ! output: Inertia, number of positive eigenvalues for symmetric indefinite matrices
      ! p%iparm(23)      ! output: Inertia, number of negative eigenvalues for symmetric indefinite matrices
        p%iparm(24) = 1  ! Parallel factorization control (default: classic)
        p%iparm(25) = 2  ! Parallel forward/backward solve control (default: sequential)
      ! p%iparm(27) = 0  ! Matrix checker (default: off)
      ! p%iparm(28) = 0  ! Single or double precision (default: double precision)
      ! ....

      p%factorized = .false.
    end associate
  end function

  subroutine destruct_pardiso_handle(h)
    integer, intent(inout) :: h
      !! pardiso handle (index)

    ! local variables
    integer(SPARSE_IDX) :: phase, error, ia(0), ja(0), perm(0), nrhs
    real                :: a(0), b(0), x(0)
    complex             :: ca(0), cb(0), cx(0)

    associate (p => pardiso_handles(h))
      if (any(p%pt%dummy /= 0)) then
        ! termination and release of memory
        nrhs = int(0, kind = SPARSE_IDX)
        phase = -1 ! release internal memory
        if (p%mtype == 11) then
          call PARDISO_ROUTINE{}(p%pt, p%maxfct, p%mnum, p%mtype, phase, p%nrows, a, ia, ja, perm, nrhs, &
            &                    p%iparm, p%msglvl, b, x, error)
        else if (p%mtype == 13) then
          call PARDISO_ROUTINE{}(p%pt, p%maxfct, p%mnum, p%mtype, phase, p%nrows, ca, ia, ja, perm, nrhs, &
            &                    p%iparm, p%msglvl, cb, cx, error)
        end if

        if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))
      end if
      p%factorized = .false.
    end associate

    !$omp critical (omp_pardiso_handles)
    call pardiso_free_handles%push_back(h)
    !$omp end critical (omp_pardiso_handles)

    h = 0
  end subroutine

  subroutine pardiso_factorize_r(h, ia, ja, a)
    !! Factorize real matrix with PARDISO.

    integer,             intent(in) :: h
      !! pardiso handle (index)
    integer(SPARSE_IDX), intent(in) :: ia(:)
    integer, target,     intent(in) :: ja(:)
    real,                intent(in) :: a(:)

    ! local variables
    integer(SPARSE_IDX)          :: error, phase, perm(0), nrhs
    integer(SPARSE_IDX), pointer :: japtr(:)
    real                         :: dum(0)

    associate (p => pardiso_handles(h))
      ! make sure the matrix type is real
      m4_assert(p%mtype == 11)

      ! convert ja if necessary
      m4_ifelse(m4_intsize,m4_idxsize,{japtr => ja},{
        allocate (japtr(size(ja, kind = SPARSE_IDX)))
        japtr = int(ja, kind = SPARSE_IDX)
      })

      nrhs = int(0, kind = SPARSE_IDX)

      ! reordering and symbolic factorization
      phase = 11 ! analysis
      call PARDISO_ROUTINE{}(p%pt, p%maxfct, p%mnum, p%mtype, phase, p%nrows, &
        &                    a, ia, japtr, perm, nrhs, p%iparm, p%msglvl, dum, dum, error)
      if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))

      ! factorization
      phase = 22 ! numerical factorization
      call PARDISO_ROUTINE{}(p%pt, p%maxfct, p%mnum, p%mtype, phase, p%nrows, &
        &                    a, ia, japtr, perm, nrhs, p%iparm, p%msglvl, dum, dum, error)
      if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))
      p%factorized = .true.

      ! free memory for ja copy
      m4_ifelse(m4_intsize,m4_idxsize,,{deallocate(japtr)})
    end associate
  end subroutine

  subroutine pardiso_factorize_c(h, ia, ja, a)
    !! Factorize complex matrix with PARDISO.

    integer,             intent(in) :: h
      !! pardiso handle (index)
    integer(SPARSE_IDX), intent(in) :: ia(:)
    integer, target,     intent(in) :: ja(:)
    complex,             intent(in) :: a(:)

    ! local variables
    integer(SPARSE_IDX)          :: error, phase, perm(0), nrhs
    integer(SPARSE_IDX), pointer :: japtr(:)
    complex                      :: dum(0)

    associate (p => pardiso_handles(h))
      ! make sure matrix type is complex
      m4_assert(p%mtype == 13)

      ! convert ja if necessary
      m4_ifelse(m4_intsize,m4_idxsize,{japtr => ja},{
        allocate (japtr(size(ja, kind = SPARSE_IDX)))
        japtr = int(ja, kind = SPARSE_IDX)
      })

      nrhs = int(0, kind = SPARSE_IDX)

      ! reordering and symbolic factorization
      phase = 11 ! analysis
      call PARDISO_ROUTINE{}(p%pt, p%maxfct, p%mnum, p%mtype, phase, p%nrows, &
        &                    a, ia, japtr, perm, nrhs, p%iparm, p%msglvl, dum, dum, error)
      if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))

      ! factorization
      phase = 22 ! numerical factorization
      call PARDISO_ROUTINE{}(p%pt, p%maxfct, p%mnum, p%mtype, phase, p%nrows, &
        &                    a, ia, japtr, perm, nrhs, p%iparm, p%msglvl, dum, dum, error)
      if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))
      p%factorized = .true.

      ! free memory for ja copy
      m4_ifelse(m4_intsize,m4_idxsize,,{deallocate(japtr)})
    end associate
  end subroutine

  subroutine pardiso_solve_r(h, ia, ja, a, b, x)
    !! Solve real system with PARDISO.

    integer,             intent(in)  :: h
      !! pardiso handle (index)
    integer(SPARSE_IDX), intent(in)  :: ia(:)
    integer, target,     intent(in)  :: ja(:)
    real,                intent(in)  :: a(:)
    real,                intent(in)  :: b(:)
    real,                intent(out) :: x(:)

    ! local variables
    integer(SPARSE_IDX)          :: error, phase, perm(0)
    integer(SPARSE_IDX), pointer :: japtr(:)
    real                         :: b_(size(b))

    associate (p => pardiso_handles(h))
      ! make sure matrix type is real and matrix is factorized
      m4_assert(p%mtype == 11)
      m4_assert(p%factorized)

      ! convert ja if necessary
      m4_ifelse(m4_intsize,m4_idxsize,{japtr => ja},{
        allocate (japtr(size(ja, kind = SPARSE_IDX)))
        japtr = int(ja, kind = SPARSE_IDX)
      })

      ! copy data
      b_ = b

      ! back substitution and iterative refinement
      phase = 33 ! solve, iterative refinement
      call PARDISO_ROUTINE{}(p%pt, p%maxfct, p%mnum, p%mtype, phase, p%nrows, &
        &                    a, ia, japtr, perm, size(b, kind = SPARSE_IDX)/p%nrows, p%iparm, p%msglvl, b_, x, error)
      if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))

      ! free memory for ja copy
      m4_ifelse(m4_intsize,m4_idxsize,,{deallocate(japtr)})
    end associate
  end subroutine

  subroutine pardiso_solve_c(h, ia, ja, a, b, x)
    !! Solve complex system with PARDISO.

    integer,             intent(in)  :: h
      !! pardiso handle (index)
    integer(SPARSE_IDX), intent(in)  :: ia(:)
    integer, target,     intent(in)  :: ja(:)
    complex,             intent(in)  :: a(:)
    complex,             intent(in)  :: b(:)
    complex,             intent(out) :: x(:)

    ! local variables
    integer(SPARSE_IDX)          :: error, phase, perm(0)
    integer(SPARSE_IDX), pointer :: japtr(:)
    complex                      :: b_(size(b))

    associate (p => pardiso_handles(h))
      ! make sure matrix type is complex and matrix is factorized
      m4_assert(p%mtype == 13)
      m4_assert(p%factorized)

      ! convert ja if necessary
      m4_ifelse(m4_intsize,m4_idxsize,{japtr => ja},{
        allocate (japtr(size(ja, kind = SPARSE_IDX)))
        japtr = int(ja, kind = SPARSE_IDX)
      })

      ! copy data
      b_ = b

      ! back substitution and iterative refinement
      phase = 33 ! solve, iterative refinement
      call PARDISO_ROUTINE{}(p%pt, p%maxfct, p%mnum, p%mtype, phase, p%nrows, &
        &                    a, ia, japtr, perm, size(b, kind = SPARSE_IDX)/p%nrows, p%iparm, p%msglvl, b_, x, error)
      if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))

      ! free memory for ja copy
      m4_ifelse(m4_intsize,m4_idxsize,,{deallocate(japtr)})
    end associate
  end subroutine

end module
