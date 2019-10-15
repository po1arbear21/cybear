#include "assert.f90.inc"

module pardiso_m
  use error_m
  implicit none

  ! include pardiso interface
  include 'mkl_pardiso.fi'

  character(len=*), parameter :: PARDISO_ERROR(-12:0) = (/                   &
    & "pardiso_64 called from 32-bit library                              ", &
    & "read/write error with OOC files                                    ", &
    & "error opening OOC files                                            ", &
    & "not enough memory for OOC                                          ", &
    & "32-bit integer overflow problem                                    ", &
    & "diagonal matrix is singular                                        ", &
    & "reordering failed                                                  ", &
    & "unclassified (internal) error                                      ", &
    & "zero pivot, numerical factorization or iterative refinement problem", &
    & "reordering problem                                                 ", &
    & "not enough memory                                                  ", &
    & "input inconsistent                                                 ", &
    & "no error                                                           "  &
    &/)

  type pardiso_handle
    type(MKL_PARDISO_HANDLE) :: pt(64)
      !! internal PARDISO handle
    integer                  :: maxfct
      !! number of matrices with identical sparse structure
    integer                  :: mnum
      !! which matrix should be used in solution phase
    integer                  :: mtype
      !! matrix type
    integer                  :: nrows
      !! number of equations
    integer                  :: msglvl
      !! message level
    integer                  :: iparm(64)
      !! PARDISO parameters
    logical                  :: factorized = .false.
      !! factorization complete flag

    ! 3 array csr
    integer, allocatable :: ia(:)
    integer, allocatable :: ja(:)
    real,    allocatable :: a(:)
    complex, allocatable :: ac(:)
  contains
    procedure :: init          => pardiso_init
    procedure :: destruct      => pardiso_destruct
    procedure :: pardiso_factorize_r
    procedure :: pardiso_factorize_c
    generic   :: factorize     => pardiso_factorize_r, pardiso_factorize_c
    procedure :: pardiso_solve_r
    procedure :: pardiso_solve_c
    generic   :: solve         => pardiso_solve_r, pardiso_solve_c
  end type

contains

  subroutine pardiso_init(this, nrows, cmplx)
    !! Inits the pardiso_handle object.
    !!
    !! Assumes unsymmetric matrices.

    class(pardiso_handle), target, intent(out) :: this
    integer,                       intent(in)  :: nrows
      !! number of rows
    logical,                       intent(in)  :: cmplx
      !! complex/real flag

    ! maximum number of factors with idential sparse structure
    this%maxfct = 1

    ! which matrix should be used in solution phase (1 <= mnum <= maxfct)
    this%mnum = 1

    ! matrix type: real unsymmetric or complex unsymmetric
    if (cmplx) then
      this%mtype = 13
    else
      this%mtype = 11
    end if

    ! number of rows
    this%nrows = nrows

    ! do not print statistical information
    this%msglvl = 0

    ! init pt, iparm
    call pardisoinit(this%pt, this%mtype, this%iparm)

    ! ======================
    !     ipar
    ! ======================

    ! this%iparm( 1) = 1  ! Use default values
    ! this%iparm( 2) = 2  ! Fill-in reducing ordering for the input matrix (default: METIS)
    ! this%iparm( 4) = 0  ! Preconditioned CGS/CG
    ! this%iparm( 5) = 0  ! User permutation
    ! this%iparm( 6) = 0  ! Write solution on x (default: right-hand side b unchanged)
    ! this%iparm( 7)      ! output: Actual number of iterative refinement steps performed
      this%iparm( 8) = 9  ! Maximal number of iterative refinement steps
    ! this%iparm(10)      ! Pivoting permutation (default: 1 or 3 depending on cmplx)
    ! this%iparm(11) = 1  ! Scaling (default: on for unsymmetric)
    ! this%iparm(12) = 0  ! Solve with transposed or conjugate transposed matrix (default: off)
    ! this%iparm(13) = 1  ! Improved accuracy using non-symmetric weighted matching (default: on for unsymmetric)
    ! this%iparm(14)      ! output: Number of perturbed pivots
    ! this%iparm(15)      ! output: Peak memory on symbolic factorization
    ! this%iparm(16)      ! output: Permanent memory on symbolic factorization
    ! this%iparm(17)      ! output: Size of factors/Peak memory on numerical factorization and solution
    ! this%iparm(18) = -1 ! Report number of non-zero elements in the factors (default: on)
    ! this%iparm(19) = 0  ! Report number of floating point operations for factorization (default: off)
    ! this%iparm(20)      ! output: Report CG/CGS diagnostics
    ! this%iparm(21) = 1  ! Pivoting for symmetric indefinite matrices (default: Bunch-Kaufman)
    ! this%iparm(22)      ! output: Inertia, number of positive eigenvalues for symmetric indefinite matrices
    ! this%iparm(23)      ! output: Inertia, number of negative eigenvalues for symmetric indefinite matrices
    ! this%iparm(24) = 0  ! Parallel factorization control (default: classic)
    ! this%iparm(25) = 0  ! Parallel forward/backward solve control (default: sequential)
    ! this%iparm(27) = 0  ! Matrix checker (default: off)
    ! ....

    this%factorized = .false.
  end subroutine

  subroutine pardiso_destruct(this)
    !! Release internal memory.

    class(pardiso_handle), intent(inout) :: this

    ! local variables
    integer :: phase, error, ia(0), ja(0), perm(0)
    real    :: a(0), b(0), x(0)
    complex :: ca(0), cb(0), cx(0)

    if (any(this%pt%dummy /= 0)) then
      ! termination and release of memory
      phase = -1 ! release internal memory
      if (this%mtype == 11) then
        call pardiso(this%pt, this%maxfct, this%mnum, this%mtype, phase, this%nrows, a, ia, ja, perm, 0, &
                     this%iparm, this%msglvl, b, x, error)
      else if (this%mtype == 13) then
        call pardiso(this%pt, this%maxfct, this%mnum, this%mtype, phase, this%nrows, ca, ia, ja, perm, 0, &
                     this%iparm, this%msglvl, cb, cx, error)
      end if

      if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))
    end if
    this%factorized = .false.
  end subroutine

  subroutine pardiso_factorize_r(this, ia0, ia1, ja, a)
    !! Factorize real matrix with PARDISO.

    class(pardiso_handle), intent(inout) :: this
    integer,               intent(in)    :: ia0(:)
    integer,               intent(in)    :: ia1(:)
    integer,               intent(in)    :: ja(:)
    real,                  intent(in)    :: a(:)

    ! local variables
    integer :: i, j, k, n, error, phase, perm(0)
    real    :: dum(0)

    ! make sure the matrix type is real
    ASSERT(this%mtype == 11)

    ! total number of elements
    n = 0
    do i = 1, this%nrows
      n = n + ia1(i) - ia0(i)
    end do

    ! allocate 3 array csr
    allocate (this%ia(this%nrows+1))
    allocate (this%ja(n))
    allocate (this%a(n))

    ! convert 4 array csr to 3 array csr
    k = 0
    do i = 1, this%nrows
      this%ia(i) = k + 1
      do j = ia0(i), ia1(i)-1
        k = k + 1
        this%ja(k) = ja(k)
        this%a( k) = a( k)
      end do
    end do
    this%ia(this%nrows+1) = k + 1

    ! reordering and symbolic factorization
    phase = 11 ! analysis
    call pardiso(this%pt, this%maxfct, this%mnum, this%mtype, phase, this%nrows, &
                 this%a, this%ia, this%ja, perm, 0, this%iparm, this%msglvl, dum, dum, error)
    if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))

    ! factorization
    phase = 22 ! numerical factorization
    call pardiso(this%pt, this%maxfct, this%mnum, this%mtype, phase, this%nrows, &
                 this%a, this%ia, this%ja, perm, 0, this%iparm, this%msglvl, dum, dum, error)
    if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))
    this%factorized = .true.
  end subroutine

  subroutine pardiso_factorize_c(this, ia0, ia1, ja, a)
    !! Factorize complex matrix with PARDISO.

    class(pardiso_handle), intent(inout) :: this
    integer,               intent(in)    :: ia0(:)
    integer,               intent(in)    :: ia1(:)
    integer,               intent(in)    :: ja(:)
    complex,               intent(in)    :: a(:)

    ! local variables
    integer :: i, j, k, n, error, phase, perm(0)
    complex :: dum(0)

    ! make sure matrix type is complex
    ASSERT(this%mtype == 13)

    ! total number of elements
    n = 0
    do i = 1, this%nrows
      n = n + ia1(i) - ia0(i)
    end do

    ! allocate 3 array csr
    allocate (this%ia(this%nrows+1))
    allocate (this%ja(n))
    allocate (this%ac(n))

    ! convert 4 array csr to 3 array csr
    k = 0
    do i = 1, this%nrows
      this%ia(i) = k + 1
      do j = ia0(i), ia1(i)-1
        k = k + 1
        this%ja(k) = ja(k)
        this%ac(k) = a( k)
      end do
    end do
    this%ia(this%nrows+1) = k + 1

    ! reordering and symbolic factorization
    phase = 11 ! analysis
    call pardiso(this%pt, this%maxfct, this%mnum, this%mtype, phase, this%nrows, &
                 this%ac, this%ia, this%ja, perm, 0, this%iparm, this%msglvl, dum, dum, error)
    if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))

    ! factorization
    phase = 22 ! numerical factorization
    call pardiso(this%pt, this%maxfct, this%mnum, this%mtype, phase, this%nrows, &
                 this%ac, this%ia, this%ja, perm, 0, this%iparm, this%msglvl, dum, dum, error)
    if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))
    this%factorized = .true.
  end subroutine

  subroutine pardiso_solve_r(this, b, x)
    !! Solve real system with PARDISO.

    class(pardiso_handle), intent(in)  :: this
    real,                  intent(in)  :: b(:)
    real,                  intent(out) :: x(:)

    ! local variables
    integer                  :: error, phase, perm(0), iparm_(64)
    real                     :: b_(size(b))
    type(MKL_PARDISO_HANDLE) :: pt_(64)

    ! make sure matrix type is real and matrix is factorized
    ASSERT(this%mtype == 11)
    ASSERT(this%factorized)

    ! copy data
    pt_    = this%pt
    iparm_ = this%iparm
    b_     = b

    ! back substitution and iterative refinement
    phase = 33 ! solve, iterative refinement
    call pardiso(pt_, this%maxfct, this%mnum, this%mtype, phase, this%nrows, &
                 this%a, this%ia, this%ja, perm, size(b)/this%nrows, iparm_, this%msglvl, b_, x, error)
    if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))

    ! make sure pardiso does not change the handle
    ASSERT(all(this%pt%dummy == pt_%dummy))
  end subroutine

  subroutine pardiso_solve_c(this, b, x)
    !! Solve complex system with PARDISO.

    class(pardiso_handle), intent(in)  :: this
    complex,               intent(in)  :: b(:)
    complex,               intent(out) :: x(:)

    ! local variables
    integer                  :: error, phase, perm(0), iparm_(64)
    complex                  :: b_(size(b))
    type(MKL_PARDISO_HANDLE) :: pt_(64)

    ! make sure matrix type is complex and matrix is factorized
    ASSERT(this%mtype == 13)
    ASSERT(this%factorized)

    ! copy data
    pt_    = this%pt
    iparm_ = this%iparm
    b_     = b

    ! back substitution and iterative refinement
    phase = 33 ! solve, iterative refinement
    call pardiso(pt_, this%maxfct, this%mnum, this%mtype, phase, this%nrows, &
                 this%ac, this%ia, this%ja, perm, size(b)/this%nrows, iparm_, this%msglvl, b_, x, error)
    if (error /= 0) call program_error("Error in pardiso: "//trim(PARDISO_ERROR(error)))

    ! make sure pardiso does not change the handle
    ASSERT(all(this%pt%dummy == pt_%dummy))
  end subroutine

end module