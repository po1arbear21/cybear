#include "../macro.f90.inc"

module preconditioner_m

  use error_m,         only: assert_failed, program_error
  use iso_fortran_env, only: real64
  use matop_m,         only: matop_real
  use matrix_m,        only: sparse_real
  use sparse_idx_m,    only: SPARSE_IDX
  use util_m,          only: int2str

  implicit none

  private
  public preconditioner
  public mkl_ilu0

  type, extends(matop_real), abstract :: preconditioner
    real, allocatable :: tmp(:)
      !! temp array needed for solving
  contains
    procedure(prec_apply    ), deferred :: apply
    procedure(prec_factorize), deferred :: factorize
  end type

  type mkl_ilu0_options
    logical :: abort_zero_diag     = .true.                                                   ! ipar 31
      !! abort when zero diagonal element occurs? or replace 0 value with diag_threshold?
    real    :: zero_diag_threshold = 1e-16                                                    ! dpar 31
      !! (only used when abort_zero_diag is false)
      !! diagonal element is interpreted as 0 when below this threshold
    real    :: zero_diag_replace   = 1e-10                                                    ! dpar 32
      !! (only used when abort_zero_diag is false)
      !! replace diagonal element with this value when below threshold
  end type

  type, extends(preconditioner) :: mkl_ilu0
    !! ilu0 preconditioner
    !!
    !! exact matrix:   A = (ia,ja,a)
    !! preconditioner: B = (ia,ja,b)
    !!      (ia,ja) can be saved via pointer to A or by copying data

    type(sparse_real), pointer :: A_ptr => null()
      !! exact matrix A (either pointer or copy)
    type(sparse_real)          :: A_copy
      !! exact matrix A (either pointer or copy)
    real, allocatable          :: b(:)
      !! ilu0 matrix B: values b

    type(mkl_ilu0_options) :: opt
      !! options
  contains
    generic   :: init      => mkl_ilu0_init
    procedure :: apply     => mkl_ilu0_apply
    procedure :: factorize => mkl_ilu0_factorize
    procedure :: exec1     => mkl_ilu0_exec1
    procedure :: exec2     => mkl_ilu0_exec2

    procedure, private :: mkl_ilu0_init
  end type

  ! preconditioner
  abstract interface
    subroutine prec_apply(this, ipar, dpar)
      import :: preconditioner
      class(preconditioner), intent(in)    :: this
      integer,               intent(inout) :: ipar(128)
      real,                  intent(inout) :: dpar(128)
    end subroutine

    subroutine prec_factorize(this, ipar, dpar)
      import :: preconditioner
      class(preconditioner), intent(inout) :: this
      integer,               intent(in)    :: ipar(128)
      real,                  intent(in)    :: dpar(128)
    end subroutine
  end interface

  ! mkl ilu routines
  interface
    ! - copied from /opt/exports/intel20/compilers_and_libraries_2020.1.217/linux/mkl/include/mkl_rci.fi
    ! - added kind information and changed formatting
    ! - including interfaces by '#include "mkl_rci.fi"' doesnt work for gfortran b.c. real data type is in interfaces
    !   defined as double precision which is real(kind=16) when using "-fdefault-real-8" compiler flag
    subroutine dcsrilu0(n, a, ia, ja, alu, ipar, dpar, ierr)
      import :: real64, SPARSE_IDX
      integer             :: n
      real(real64)        :: a(*)
      integer(SPARSE_IDX) :: ia(*)
      integer             :: ja(*)
      real(real64)        :: alu(*)
      integer             :: ipar(*)
      real(real64)        :: dpar(*)
      integer             :: ierr
    end subroutine

    subroutine dcsrilut(n, a, ia, ja, alut, ialut, jalut, tol, maxfil, ipar, dpar, ierr)
      import :: real64, SPARSE_IDX
      integer             :: n
      real(real64)        :: a(*)
      integer(SPARSE_IDX) :: ia(*)
      integer             :: ja(*)
      real(real64)        :: alut(*)
      integer(SPARSE_IDX) :: ialut(*)
      integer             :: jalut(*)
      real(real64)        :: tol
      integer             :: maxfil
      integer             :: ipar(*)
      real(real64)        :: dpar(*)
      integer             :: ierr
    end subroutine
  end interface

  ! mkl spblas routines
  interface
    ! - copied from /opt/exports/intel20/compilers_and_libraries_2020.1.217/linux/mkl/include/mkl_spblas.fi
    ! - added kind information and changed formatting
    ! - including interfaces by '#include "mkl_spblas.fi"' doesnt work for gfortran b.c. real data type is in interfaces
    !   defined as double precision which is real(kind=16) when using "-fdefault-real-8" compiler flag

    subroutine mkl_dcsrtrsv(uplo, transa, diag, m, a, ia, ja, x, y)
      import :: real64, SPARSE_IDX
      character           :: uplo
      character           :: transa
      character           :: diag
      integer             :: m
      real(real64)        :: a(*)
      integer(SPARSE_IDX) :: ia(*)
      integer             :: ja(*)
      real(real64)        :: x(*)
      real(real64)        :: y(*)
    end
  end interface

  character(*), parameter :: MKL_ILU0_ERROR(-106:-101) = [                              &
    & "ja not in ascending order                                                     ", &
    & "system size: n <= 0                                                           ", &
    & "insufficient memory                                                           ", &
    & "small diagonal element: overflow or bad approximation possible                ", &
    & "zero diagonal encountered                                                     ", &
    & "error. at least one diagonal element is omitted from the matrix in CSR3 format"  ]

contains

  subroutine mkl_ilu0_init(this, A, copy)
    !! init ilu0 preconditioner
    !!
    !! exact matrix:   A = (ia,ja,a)
    !! preconditioner: B = (ia,ja,b)
    !!
    !! Options: Either save A by pointer or copy its data.
    !! Ask yourself if A will still be available when we solve with this preconditioner.

    class(mkl_ilu0),           intent(out) :: this
      !! ilu matrix B
    type(sparse_real), target, intent(in)  :: A
      !! exact matrix A
    logical, optional,         intent(in)  :: copy
      !! save copy of A (vs only pointer to A). (default: false)

    logical :: copy_

    ! optional arg
    copy_ = .false.
    if (present(copy)) copy_ = copy

    ! init base
    call this%init("", A%nrows)

    ! init data
    if (copy_) then
      call this%A_copy%init(A%nrows)
      this%A_copy%a  = A%a
      this%A_copy%ia = A%ia
      this%A_copy%ja = A%ja
    else
      this%A_ptr => A
    end if

    allocate (this%b(size(A%a)), this%tmp(A%nrows))
  end subroutine

  subroutine mkl_ilu0_apply(this, ipar, dpar)
    !! change gmres ipar/dpar according to options
    class(mkl_ilu0), intent(in)    :: this
    integer,         intent(inout) :: ipar(128)
      !! integer parameters. see mkl doc
    real,            intent(inout) :: dpar(128)
      !! real parameters. see mkl doc

    if (this%opt%abort_zero_diag) then
      ipar(31) = 0

    else
      ipar(31) = 1
      dpar(31) = this%opt%zero_diag_threshold
      dpar(32) = this%opt%zero_diag_replace
    end if
  end subroutine

  subroutine mkl_ilu0_factorize(this, ipar, dpar)
    !! compute ilu0 factorization
    class(mkl_ilu0), intent(inout) :: this
    integer,         intent(in)    :: ipar(128)
      !! integer parameters. see mkl doc
    real,            intent(in)    :: dpar(128)
      !! real parameters. see mkl doc

    integer :: ierr

    if (associated(this%A_ptr)) then
      call dcsrilu0(this%nrows, this%A_ptr%a,  this%A_ptr%ia,  this%A_ptr%ja,  this%b, ipar, dpar, ierr)
    else
      call dcsrilu0(this%nrows, this%A_copy%a, this%A_copy%ia, this%A_copy%ja, this%b, ipar, dpar, ierr)
    end if

    if (ierr /= 0) call program_error(MKL_ILU0_ERROR(ierr))
  end subroutine

  subroutine mkl_ilu0_exec1(this, x, y)
    !! solve for given vector
    !!
    !!    this*y = x   =>   y = this^-1 * x
    class(mkl_ilu0), intent(in)  :: this
    real, target,    intent(in)  :: x(:)
      !! input vector
    real, target,    intent(out) :: y(:)
      !! output vector

    if (associated(this%A_ptr)) then
      call mkl_dcsrtrsv('L', 'N', 'U', this%nrows, this%b, this%A_ptr%ia,  this%A_ptr%ja,  x,        this%tmp)
      call mkl_dcsrtrsv('U', 'N', 'N', this%nrows, this%b, this%A_ptr%ia,  this%A_ptr%ja,  this%tmp, y       )
    else
      call mkl_dcsrtrsv('L', 'N', 'U', this%nrows, this%b, this%A_copy%ia, this%A_copy%ja, x,        this%tmp)
      call mkl_dcsrtrsv('U', 'N', 'N', this%nrows, this%b, this%A_copy%ia, this%A_copy%ja, this%tmp, y       )
    end if
  end subroutine

  subroutine mkl_ilu0_exec2(this, x, y)
    !! solve using this preconditioner given multiple rhs
    !!
    !!    this*y = x   =>   y = this^-1 * x
    class(mkl_ilu0), intent(in)  :: this
    real, target,    intent(in)  :: x(:,:)
      !! input matrix
    real, target,    intent(out) :: y(:,:)
      !! output matrix

    integer :: i

    ASSERT(size(x, dim=1) == this%nrows)
    ASSERT(all(shape(x) == shape(y))   )

    do i = 1, size(x, dim=2)
      call this%exec1(x(:,i), y(:,i))
    end do
  end subroutine

end module
