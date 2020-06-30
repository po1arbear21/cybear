#include "../util/macro.f90.inc"

module ilupack_m
  !! wrapper around ILUPACK. uses interfaces for ILUPACK defined in ilupack_interf_m.

  use ilupack_interfaces_m
  use error_m
  use util_m

  implicit none

  private
  public :: ilupack
  public :: ilupack_opt

  type ilupack_opt
    !! ilupack options.
    !! usage:
    !!    1. call init: sets default parameters
    !!    [2. overwrite values manually]

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
    integer       :: elbow
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

  contains
    procedure :: init  => ilupack_opt_init
    procedure :: print => ilupack_opt_print
  end type

  type ilupack
    !! ilupack wrapper.
    !! usage:
    !!    1. call init: sets default parameters
    !!    [1b. overwrite opt values manually]
    !!    2. call factor
    !!    [2a. call info: writes info about multilevel structure]
    !!    [2b. get fill-in-factor]
    !!    3. call solve
    !!    [3b. call solver for multiple rhs]
    !!    4. call delete

    type(ilupack_opt) :: opt

    integer           :: n
      !! size of the system. matrix is of (nxn) size.
    integer           :: it
      !! number of iterations at solver step

    integer           :: param
      !! parameter C-pointer casted to integer
    integer           :: prec
      !! preconditioner C-pointer casted to integer

  contains
    procedure :: init   => ilupack_init
    procedure :: factor => ilupack_factor
    procedure :: info   => ilupack_info
    procedure :: nnz    => ilupack_nnz
    procedure :: solve  => ilupack_solve
    procedure :: delete => ilupack_delete
  end type

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

  function ilupack_opt_init(this, a, ia, ja) result(n)
    !! init ilupack options.

    class(ilupack_opt), intent(out) :: this
    real,               intent(in)  ::  a(:)
    integer,            intent(in)  :: ia(:)
    integer,            intent(in)  :: ja(:)
    integer                         :: n
      !! system size

    n = size(ia)-1

    call dgnlamginit(n,             ia,           ja,            a,            this%matching, &
      &              this%ordering, this%droptol, this%droptolS, this%condest, this%restol,   &
      &              this%maxit,    this%elbow,   this%lfil,     this%lfilS,   this%nrestart  )

  end function

  subroutine ilupack_opt_print(this)
    !! print ilupack options.

    class(ilupack_opt), intent(in) :: this

    print '(A)',        'ordering: "' // this%ordering // '"'
    print '(A,I)',      'elbow:     ',   this%elbow
    print '(A,I)',      'lfil:      ',   this%lfil
    print '(A,I)',      'lfilS:     ',   this%lfilS
    print '(A,I)',      'matching:  ',   this%matching
    print '(A,I)',      'maxit:     ',   this%maxit
    print '(A,I)',      'nrestart:  ',   this%nrestart
    print '(A,E20.10)', 'droptol:   ',   this%droptol
    print '(A,E20.10)', 'droptolS:  ',   this%droptolS
    print '(A,E20.10)', 'condest:   ',   this%condest
    print '(A,E20.10)', 'restol:    ',   this%restol
  end subroutine

  subroutine ilupack_init(this, a, ia, ja)
    !! init ilupack.

    class(ilupack), intent(out) :: this
    real,           intent(in)  ::  a(:)
    integer,        intent(in)  :: ia(:)
    integer,        intent(in)  :: ja(:)

    ASSERT(size(a) == size(ja))

    this%n = this%opt%init(a, ia, ja)
  end subroutine

  subroutine ilupack_factor(this, a, ia, ja)
    class(ilupack), intent(inout) :: this
    real,           intent(inout) ::  a(:)
    integer,        intent(inout) :: ia(:)
    integer,        intent(inout) :: ja(:)

    integer :: ierr

    ASSERT(size( a) == size(ja))
    ASSERT(size(ia) == this%n+1)

    ierr = dgnlamgfactor(this%param,        this%prec,                                                                &
      &                  this%n,            ia,               ja,                a,                this%opt%matching, &
      &                  this%opt%ordering, this%opt%droptol, this%opt%droptolS, this%opt%condest, this%opt%restol,   &
      &                  this%opt%maxit,    this%opt%elbow,   this%opt%lfil,     this%opt%lfilS,   this%opt%nrestart  )

    if      ((ierr > -8) .and. (ierr < 0)) then
      call program_error("Error in ilupack: " // trim(ILUPACK_FACTOR_ERROR(ierr)))
    else if ( ierr /= 0                  ) then
      call program_error("Error in ilupack: zero pivot encountered at step number " // int2str(ierr))
    end if
  end subroutine

  subroutine ilupack_info(this, a, ia, ja)
    class(ilupack), intent(in) :: this
    real,           intent(in) ::  a(:)
    integer,        intent(in) :: ia(:)
    integer,        intent(in) :: ja(:)

    ASSERT(size( a) == size(ja))
    ASSERT(size(ia) == this%n+1)

    call dgnlamginfo(this%param, this%prec, this%n, ia, ja, a)
  end subroutine

  integer function ilupack_nnz(this)
    class(ilupack), intent(in) :: this

    ilupack_nnz = dgnlamgnnz(this%param, this%prec)
  end function

  subroutine ilupack_solve(this, a, ia, ja, rhs, sol)
    class(ilupack), intent(inout) :: this
    real,           intent(in)    ::   a(:)
    integer,        intent(in)    ::  ia(:)
    integer,        intent(in)    ::  ja(:)
    real,           intent(in)    :: rhs(:)
      !! right hand side  of length n
    real,           intent(inout) :: sol(:)
      !! solution vector of length n, must be initialized on input with an initial guess

    integer :: ierr

    ASSERT(size(  a) == size(ja))
    ASSERT(size( ia) == this%n+1)
    ASSERT(size(rhs) == this%n  )
    ASSERT(size(sol) == this%n  )

    this%it = this%opt%maxit

    ierr = dgnlamgsolver(this%param,        this%prec,        rhs,               sol,                                 &
      &                  this%n,            ia,               ja,                a,                this%opt%matching, &
      &                  this%opt%ordering, this%opt%droptol, this%opt%droptolS, this%opt%condest, this%opt%restol,   &
      &                  this%it,           this%opt%elbow,   this%opt%lfil,     this%opt%lfilS,   this%opt%nrestart  )

    if      ((ierr >= -3) .and. (ierr < 0)) then
      call program_error("Error in ilupack: " // trim(ILUPACK_SOLVER_ERROR(ierr)))
    else if ( ierr <  -3                  ) then
      call program_error("Error in ilupack: solver exited with error code " // int2str(ierr))
    end if
  end subroutine

  subroutine ilupack_delete(this)
    class(ilupack), intent(inout) :: this

    call dgnlamgdelete(this%param, this%prec)
  end subroutine

end module
