#ifdef USE_ILUPACK
#include "../util/macro.f90.inc"

#ifdef IDXSIZE64
#ifdef INTSIZE32
#error "ILUPACK does not support different integer sizes!"
#endif
#endif

module ilupack_m

  use error_m
  use iso_fortran_env, only: int64
  use util_m

  implicit none

  private
  public ilupack_handle
  public create_ilupack_handle
  public destruct_ilupack_handle
  public ilupack_factorize
  public ilupack_solve
  public get_ilupack_handle_ptr

  type ilupack_handle
    ! ilupack options (can be set manually)

    character(20) :: ordering
      !! string indicating which reordering to apply
      !! multilevel orderings
      !! 'amd'           Approximate Minimum Degree
      !! 'mmd'           Minimum Degree
      !! 'rcm'           Reverse Cuthill-McKee
      !! 'metisn'        Metis multilevel nested dissection by nodes
      !! 'metise'        Metis multilevel nested dissection by edges
      !! 'pq'            ddPQ strategy by Saad
      !! default: 'amd
    real          :: elbow
      !! elbow space factor for the fill computed during the ILU
      !! default: 10
    integer       :: lfil
      !! number of fill elements per column in L resp. rows in U
      !! default: n+1
    integer       :: lfilS
      !! maximum number of nonzeros per row in the approximate Schur complement
      !! default: n+1
    integer       :: matching
      !! maximum weight matching
      !! default: 1 == is different from zero == matching turned on
    integer       :: maxit
      !! maximum number of iterative steps
      !! default: (500??)
    integer       :: nrestart
      !! restart length for GMRES
      !! default: 30
    real          :: droptol
      !! threshold for ILU
      !! default: 1e-2
    real          :: droptolS
      !! threshold for the approximate Schur complements
      !! info: should be a smaller than droptol (~10%?)
      !! default: 1e-3
    real          :: condest
      !! norm bound for the inverse factors L^{-1}, U^{-1}
      !! default: 1e2
    real          :: restol
      !! relative error for the backward error (SPD case: relative energy norm) used during the iterative solver
      !! default: sqrt(eps)
    integer       :: mixedprecision
      !! use single precision preconditioner if == 1
      !! default: 0

    ! parameters (should not be set manually)
    integer :: n
      !! size of the system. matrix is of (nxn) size.
    integer :: it
      !! number of iterations at solver step

    ! ilupack workspace variables
    integer, allocatable :: ind(:)
      !! array of size n
      !! default: 0
    integer(int64)       :: param
      !! parameter C-pointer casted to integer
    integer(int64)       :: prec
      !! preconditioner C-pointer casted to integer

    integer, allocatable :: ia(:)
    integer, allocatable :: ja(:)
    real,    allocatable :: ar(:)
    complex, allocatable :: ac(:)
    logical              :: cmplx
      !! real or complex

  contains
    procedure :: print_opts => ilupack_print_opts
    procedure :: info       => ilupack_info
    procedure :: nnz        => ilupack_nnz
  end type

#define T ilupack_handle
#define TT type(ilupack_handle)
#include "../util/vector_def.f90.inc"

  type(vector_ilupack_handle), target :: ilupack_handles
  type(vector_int)                    :: free_ilupack_handles

! external interfaces
#define T0 D
#define TT real
#include "ilupack_def.f90.inc"

#define T0 Z
#define TT complex
#include "ilupack_def.f90.inc"

  interface create_ilupack_handle
    module procedure :: create_ilupack_handle_r
    module procedure :: create_ilupack_handle_c
  end interface

  interface ilupack_solve
    module procedure :: ilupack_solve_r
    module procedure :: ilupack_solve_c
  end interface

  character(*), parameter :: ILUPACK_FACTOR_ERROR(-7:-1) = [  &
    & "buffers are too small                               ", &
    & "zero column encountered                             ", &
    & "zero row encountered                                ", &
    & "Illegal value for lfil                              ", &
    & "matrix U overflow, increase elbow and retry         ", &
    & "matrix L overflow, increase elbow and retry         ", &
    & "Error. input matrix may be wrong.                   "  ]

  character(*), parameter :: ILUPACK_SOLVER_ERROR(-3:-1) = [  &
    & "algorithm breaks down                               ", &
    & "not enough work space                               ", &
    & "too many iterations                                 "  ]

contains

#define T ilupack_handle
#define TT type(ilupack_handle)
#include "../util/vector_imp.f90.inc"

  function create_ilupack_handle_rc() result(h)
    !! get free ilupack integer handle

    integer :: h

    if (.not. allocated(ilupack_handles%d)) call ilupack_handles%init(0, c = 4)

    if (free_ilupack_handles%n > 0) then
      h = free_ilupack_handles%d(free_ilupack_handles%n)
      call free_ilupack_handles%resize(free_ilupack_handles%n-1)

    else
      block
        type(ilupack_handle) :: p

        call ilupack_handles%push(p)
      end block
      h = ilupack_handles%n
    end if
  end function

  function create_ilupack_handle_r(ia, ja, a) result(h)
    !! create real ilupack handle

    integer, intent(in) :: ia(:)
    integer, intent(in) :: ja(:)
    real,    intent(in) :: a(:)
    integer             :: h
      !! return ilupack integer handle

    ASSERT(size(a) == size(ja))

    h = create_ilupack_handle_rc()

    associate (p => ilupack_handles%d(h))
      p%cmplx = .false.

      p%n = size(ia) - 1
      allocate (p%ind(p%n), source = 0)

      ! copy matrix
      p%ia = ia
      p%ja = ja
      p%ar = a

      call dgnlamginit(p%n,              p%ia,      p%ja,       p%ar,      p%matching, &
        &              p%ordering,       p%droptol, p%droptolS, p%condest, p%restol,   &
        &              p%maxit,          p%elbow,   p%lfil,     p%lfilS,   p%nrestart, &
        &              p%mixedprecision, p%ind                                         )
    end associate
  end function

  function create_ilupack_handle_c(ia, ja, a) result(h)
    !! create complex ilupack handle

    integer, intent(in) :: ia(:)
    integer, intent(in) :: ja(:)
    complex, intent(in) :: a(:)
    integer             :: h
      !! return ilupack integer handle

    ASSERT(size(a) == size(ja))

    h = create_ilupack_handle_rc()

    associate (p => ilupack_handles%d(h))
      p%cmplx = .true.

      p%n = size(ia) - 1
      allocate (p%ind(p%n), source = 0)

      ! copy matrix
      p%ia = ia
      p%ja = ja
      p%ac = a

      call zgnlamginit(p%n,              p%ia,      p%ja,       p%ac,      p%matching, &
        &              p%ordering,       p%droptol, p%droptolS, p%condest, p%restol,   &
        &              p%maxit,          p%elbow,   p%lfil,     p%lfilS,   p%nrestart, &
        &              p%mixedprecision, p%ind                                         )
    end associate
  end function

  subroutine destruct_ilupack_handle(h)
    integer, intent(inout) :: h
      !! ilupack integer handle

    associate (p => ilupack_handles%d(h))
      if (p%cmplx) then
        call zgnlamgdelete(p%param, p%prec)
      else
        call dgnlamgdelete(p%param, p%prec)
      end if

      if (allocated(p%ind)) deallocate (p%ind)
      if (allocated(p%ia )) deallocate (p%ia)
      if (allocated(p%ja )) deallocate (p%ja)
      if (allocated(p%ar )) deallocate (p%ar)
      if (allocated(p%ac )) deallocate (p%ac)
    end associate

    if (.not. allocated(free_ilupack_handles%d)) call free_ilupack_handles%init(0, c = 4)
    call free_ilupack_handles%push(h)
    h = 0
  end subroutine

  subroutine ilupack_factorize(h)
    integer, intent(in) :: h
      !! ilupack integer handle

    integer :: ierr

    associate (p => ilupack_handles%d(h))
      if (p%cmplx) then
        ierr = zgnlamgfactor(p%param,          p%prec,                                       &
          &                  p%n,              p%ia,      p%ja,       p%ac,      p%matching, &
          &                  p%ordering,       p%droptol, p%droptolS, p%condest, p%restol,   &
          &                  p%maxit,          p%elbow,   p%lfil,     p%lfilS,   p%nrestart, &
          &                  p%mixedprecision, p%ind                                         )

      else
        ierr = dgnlamgfactor(p%param,          p%prec,                                       &
          &                  p%n,              p%ia,      p%ja,       p%ar,      p%matching, &
          &                  p%ordering,       p%droptol, p%droptolS, p%condest, p%restol,   &
          &                  p%maxit,          p%elbow,   p%lfil,     p%lfilS,   p%nrestart, &
          &                  p%mixedprecision, p%ind                                         )
      end if
    end associate

    if      ((ierr > -8) .and. (ierr < 0)) then
      call program_error("Error in ilupack: " // trim(ILUPACK_FACTOR_ERROR(ierr)))
    else if ( ierr /= 0                  ) then
      call program_error("Error in ilupack: zero pivot encountered at step number " // int2str(ierr))
    end if
  end subroutine

  subroutine ilupack_solve_r(h, b, x)
    !! Solve real system with ILUPACK.

    integer, intent(in)  :: h
      !! ilupack handle (index)
    real,    intent(in)  :: b(:)
      !! one or multiple right-hand sides as flat array
    real,    intent(out) :: x(:)
      !! output solution as flat array

    integer :: ierr, nrhs, i, i0, i1

    associate (p => ilupack_handles%d(h))
      ASSERT(size(x) == size(b))
      ASSERT(mod(size(b), p%n) == 0)

      nrhs = size(b) / p%n

      i1 = 0
      do i = 1, nrhs
        i0 = i1 + 1
        i1 = i1 + p%n

        ! initial guess
        call dgnlamgsol(p%param, p%prec, b(i0:i1), x(i0:i1), p%n)

        p%it = p%maxit
        ierr = dgnlamgsolver(p%param,          p%prec,    b(i0:i1),   x(i0:i1),              &
          &                  p%n,              p%ia,      p%ja,       p%ar,      p%matching, &
          &                  p%ordering,       p%droptol, p%droptolS, p%condest, p%restol,   &
          &                  p%it,             p%elbow,   p%lfil,     p%lfilS,   p%nrestart, &
          &                  p%mixedprecision, p%ind                                         )

        if      ((ierr >= -3) .and. (ierr < 0)) then
          call program_error("Error in ilupack: " // trim(ILUPACK_SOLVER_ERROR(ierr)))
        else if ( ierr /=  0                  ) then
          call program_error("Error in ilupack: solver exited with error code " // int2str(ierr))
        end if
      end do
    end associate
  end subroutine

  subroutine ilupack_solve_c(h, b, x)
    !! Solve complex system with ILUPACK.

    integer, intent(in)  :: h
      !! ilupack handle (index)
    complex, intent(in)  :: b(:)
      !! one or multiple right-hand sides as flat array
    complex, intent(out) :: x(:)
      !! output solution as flat array

    integer :: ierr, nrhs, i, i0, i1

    associate (p => ilupack_handles%d(h))
      ASSERT(size(x) == size(b))
      ASSERT(mod(size(b), p%n) == 0)

      nrhs = size(b) / p%n

      i1 = 0
      do i = 1, nrhs
        i0 = i1 + 1
        i1 = i1 + p%n

        ! initial guess
        call zgnlamgsol(p%param, p%prec, b(i0:i1), x(i0:i1), p%n)

        p%it = p%maxit
        ierr = zgnlamgsolver(p%param,          p%prec,    b(i0:i1),   x(i0:i1),              &
          &                  p%n,              p%ia,      p%ja,       p%ac,      p%matching, &
          &                  p%ordering,       p%droptol, p%droptolS, p%condest, p%restol,   &
          &                  p%it,             p%elbow,   p%lfil,     p%lfilS,   p%nrestart, &
          &                  p%mixedprecision, p%ind                                         )

        if      ((ierr >= -3) .and. (ierr < 0)) then
          call program_error("Error in ilupack: " // trim(ILUPACK_SOLVER_ERROR(ierr)))
        else if ( ierr /=  0                  ) then
          call program_error("Error in ilupack: solver exited with error code " // int2str(ierr))
        end if
      end do
    end associate
  end subroutine

  subroutine get_ilupack_handle_ptr(h, ilu)
    !! get pointer to ilupack object handle (change options, print info etc.)

    integer,                       intent(in)  :: h
      !! integer ilupack handle
    type(ilupack_handle), pointer, intent(out) :: ilu
      !! output ilupack handle object

    ilu => ilupack_handles%d(h)
  end subroutine

  subroutine ilupack_print_opts(this)
    !! Print ilupack options.

    class(ilupack_handle), intent(in) :: this

    print '(A)',         'ordering:      "' // this%ordering // '"'
    print '(A,ES24.16)', 'elbow:          ',   this%elbow
    print '(A,I24)',     'lfil:           ',   this%lfil
    print '(A,I24)',     'lfilS:          ',   this%lfilS
    print '(A,I24)',     'matching:       ',   this%matching
    print '(A,I24)',     'maxit:          ',   this%maxit
    print '(A,I24)',     'nrestart:       ',   this%nrestart
    print '(A,ES24.16)', 'droptol:        ',   this%droptol
    print '(A,ES24.16)', 'droptolS:       ',   this%droptolS
    print '(A,ES24.16)', 'condest:        ',   this%condest
    print '(A,ES24.16)', 'restol:         ',   this%restol
    print '(A,I24)',     'mixedprecision: ',   this%mixedprecision
  end subroutine

  subroutine ilupack_info(this)
    !! Print ILUPACK information.

    class(ilupack_handle), intent(in) :: this

    if (this%cmplx) then
      call zgnlamginfo(this%param, this%prec, this%n, this%ia, this%ja, this%ac)
    else
      call dgnlamginfo(this%param, this%prec, this%n, this%ia, this%ja, this%ar)
    end if
  end subroutine

  function ilupack_nnz(this) result(nnz)
    !! Number of non-zero entries.

    class(ilupack_handle), intent(in) :: this
    integer                           :: nnz

    if (this%cmplx) then
      nnz = zgnlamgnnz(this%param, this%prec)
    else
      nnz = dgnlamgnnz(this%param, this%prec)
    end if
  end function

end module
#endif
