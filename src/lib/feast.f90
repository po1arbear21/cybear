#include "../util/macro.f90.inc"

module feast_m
  use error_m,         only: assert_failed
  use iso_fortran_env, only: real64, real32
  use util_m,          only: int2str

  implicit none

  private
  public feast
  public feast_option

  type feast_option
    integer :: runtime_comments  = 0                                                   ! fpm(1)
      !! print runtime comments
      !!         0: off
      !!         1: on-screen
      !!      -n<0: Write/Append comments in the file: feast<n>.log
    real    :: atol              = 1e-12                                               ! fpm(3)
      !! stopping convergence criteria
    integer :: max_refinement    = 20                                                  ! fpm(4)
      !! maximum numer of FEAST refinement loops
    logical :: trace_convergence = .false.                                             ! fpm(6)
      !! Convergence criteria (for solutions in the search contour)
      !!    true:  Using relative error on the trace epsout, i.e. epsout < epsilon
      !!    false: Using relative residual res, i.e. max_i res(i) < epsilon
    integer :: n_contour         = 16                                                  ! fpm(8)
      !! #contour points for non-Herm./poly. FEAST (full-contour)
      !! if gauss_integration: values permitted (2 to 40, 48, 64, 80, 96, 112)
      !! else:                 all values permitted (>2)
    logical :: gauss_integration = .false.                                             ! fpm(16)
      !! Gauss or trapezoidal integration
    complex :: mid               = 0
      !! coordinate center
    real    :: rad(2)            = [1, 1]
      !! radius in horizontal/vertical directions (before rotation is applied)
    real    :: ang               = 0
      !! rotation angle in degree
    integer :: N
      !! system size
    integer :: M0
      !! overestimate for number of evals in contour
  contains
    procedure :: init  => feast_option_init
    procedure :: apply => feast_option_apply
    procedure :: get   => feast_option_get
  end type

  ! FEAST internal interfaces
  interface
    subroutine feastinit(fpm)
      integer, dimension(*) :: fpm
    end subroutine

    subroutine zfeast_grci(ijob, N, Ze, work1, work2, Aq, Bq, fpm, epsout, loop, Emid, r, M0, E, X, M, res, info)
      import real64
      integer         :: ijob
      integer         :: N
      complex(real64) :: Ze
      complex(real64) :: work1(N,*)
      complex(real64) :: work2(N,*)
      complex(real64) :: Aq(M0,*)
      complex(real64) :: Bq(M0,*)
      integer         :: fpm(*)
      real(real64)    :: epsout
      integer         :: loop
      complex(real64) :: Emid
      real(real64)    :: r
      integer         :: M0
      complex(real64) :: E(*)
      complex(real64) :: X(N,*)
      integer         :: M
      real(real64)    :: res(*)
      integer         :: info
    end subroutine
  end interface

  ! user-defined procedures
  interface
    subroutine feast_factorize(Ze, ctrans)
      !! Factorize Ze*B-A
      import :: real64
      complex(real64),   intent(in) :: Ze
      logical, optional, intent(in) :: ctrans
        !! factorize (Ze*B-A)^H. default: .false.
    end subroutine

    subroutine feast_solve(xb, ctrans)
      !! Solve (Ze*B-A)*x=b
      import :: real64
      complex(real64),   intent(inout) :: xb(:,:)
        !! input:  b
        !! output: x
      logical, optional, intent(in)    :: ctrans
        !! solve (Ze*B-A)^H*x=b. default: .false.
    end subroutine

    subroutine feast_mulvec(x, y, ctrans)
      !! Matrix-vector multiplication: y <- M*x
      import :: real64
      complex(real64),   intent(in)  :: x(:,:)
        !! input vector x
      complex(real64),   intent(out) :: y(:,:)
        !! output vector y
      logical, optional, intent(in)  :: ctrans
        !! y <- M^H*x. default: .false.
    end subroutine
  end interface

contains

  subroutine feast(opt, fact_Az, solve_Az, mulvec_A, eval, evecR, evecL, mulvec_B)
    !! compute eigenvalue decomposition in given contour
    type(feast_option),                     intent(in)  :: opt
      !! feast options
    procedure(feast_factorize)                          :: fact_Az
      !! factorize Az=(Ze*B-A)
    procedure(feast_solve)                              :: solve_Az
      !! solve Az*x=b
    procedure(feast_mulvec)                             :: mulvec_A
      !! multiply y<-A*x
    complex(real64), allocatable,           intent(out) :: eval(:)
      !! eigenvalues
    complex(real64), allocatable, optional, intent(out) :: evecR(:,:)
      !! right eigenvectors
    complex(real64), allocatable, optional, intent(out) :: evecL(:,:)
      !! left eigenvectors
    procedure(feast_mulvec),      optional              :: mulvec_B
      !! multiply y<-B*x. default: B=I

    complex(real64)              :: Ze, Emid
    complex(real64), allocatable :: work1(:,:), work2(:,:), Aq(:,:), Bq(:,:), E(:), X(:,:)
    integer                      :: fpm(64), info, M, M0, loop, i, j, ijob
    real(real64)                 :: epsout, r
    real(real64),    allocatable :: res(:)

    ! set options
    call feastinit(fpm)
    call opt%apply(fpm)
    fpm(15) = merge(0, 1, present(evecL))

    ! initialization
    ijob = -1
    M0   = opt%M0
    allocate (work1(opt%N,M0), work2(opt%N,M0), Aq(M0,M0), Bq(M0,M0), E(M0), res(M0))
    allocate (X(opt%N,M0 * merge(2, 1, present(evecL))))

    ! get parameters
    call opt%get(Emid, r)

    do while (ijob /= 0)
      call zfeast_grci(ijob, opt%N, Ze, work1, work2, Aq, Bq, fpm, epsout, loop, Emid, r, M0, E, X, M, res, info)
      if (.not. check_info(info)) return

      select case (ijob)
        case (10)
          ! Factorize the complex matrix Az <= (ZeB-A)  --  or factorize a preconditioner of ZeB-A
          ! REMARK: Az can be formed and factorized using single precision arithmetic
          call fact_Az(Ze)

        case (11)
          ! Solve the linear system with fpm(23) rhs; Az * Qz=work2(1:N,1:fpm(23))
          ! Result (in place) in work2 <= Qz(1:N,1:fpm(23))
          ! REMAKRS:  - Solve can be perfomed in single precision
          !           - Low accuracy iterative solvers are ok
          call solve_Az(work2(:,:fpm(23)))

        case (20)
          ! [Optional: *only if* needed by case (21)]
          ! Factorize the complex matrix Az^H
          ! REMARKS:  -The matrix Az from case(10) cannot be overwritten
          !           -case(20) becomes obsolete if the solve in case(21) can be perfomed
          !            by reusing the factorization in case(10)
          call fact_Az(Ze, ctrans=.true.)

        case (21)
          ! Solve the linear system with fpm(23) rhs; Az^H * Qz=work2(1:N,1:fpm(23))
          ! Result (in place) in work2 <= Qz(1:N,1:fpm(23))
          call solve_Az(work2(:,:fpm(23)), ctrans=.true.)

        case (30)
          ! Perform multiplication A   ∗ X(1:N,i:j) result in work1(1:N,i:j) where i=fpm(24) and j=fpm(24)+fpm(25)-1
          i = fpm(24)
          j = fpm(24)+fpm(25)-1
          call mulvec_A(X(:,i:j), work1(:,i:j))

        case (31)
          ! Perform multiplication A^H ∗ X(1:N,i:j) result in work1(1:N,i:j) where i=fpm(34) and j=fpm(34)+fpm(35)-1
          i = fpm(34)
          j = fpm(34)+fpm(35)-1
          call mulvec_A(X(:,i:j), work1(:,i:j), ctrans=.true.)

        case (40)
          ! Perform multiplication B   ∗ X(1:N,i:j) result in work1(1:N,i:j) where i=fpm(24) and j=fpm(24)+fpm(25)-1
          ! REMARK: user must set work1(1:N,i:j)=X(1:N,i:j) if B=I
          i = fpm(24)
          j = fpm(24)+fpm(25)-1
          if (present(mulvec_B)) then
            call mulvec_B(X(:,i:j), work1(:,i:j))
          else
            work1(:,i:j) = X(:,i:j)
          end if

        case (41)
          ! Perform multiplication B^H ∗ X(1:N,i:j) result in work1(1:N,i:j) where i=fpm(34) and j=fpm(34)+fpm(35)-1
          ! REMARK: user must set work1(1:N,i:j)=X(1:N,i:j) if B=I
          i = fpm(34)
          j = fpm(34)+fpm(35)-1
          if (present(mulvec_B)) then
            call mulvec_B(X(:,i:j), work1(:,i:j), ctrans=.true.)
          else
            work1(:,i:j) = X(:,i:j)
          end if
      end select
    end do

    ! collect results
    eval = E(:M)
    if (present(evecR)) evecR = X(:,    :M)
    if (present(evecL)) evecL = X(:,M0+1:M0+M)
  end subroutine

  subroutine feast_option_init(this, N, M0)
    !! init feast options
    class(feast_option), intent(out) :: this
    integer,             intent(in)  :: N
      !! system size
    integer,             intent(in)  :: M0
      !! overestimate for number of evals in contour

    ASSERT(N  >= M0)
    ASSERT(M0 >   0)

    this%N  = N
    this%M0 = M0
  end subroutine

  subroutine feast_option_apply(this, fpm)
    !! applies values from options object to feast integer parameter array
    class(feast_option), intent(in)  :: this
    integer,             intent(out) :: fpm(64)
      !! feast parameter array

    fpm( 1) = this%runtime_comments
    fpm( 3) = -floor(log10(this%atol))              ! fpm3: eps = 10^-fpm3
    fpm( 4) = this%max_refinement
    fpm( 6) = merge(0, 1, this%trace_convergence)
    fpm( 8) = this%n_contour
    fpm(16) = merge(0, 1, this%gauss_integration)
    fpm(18) = nint(this%rad(2)/this%rad(1)*100)     ! fpm18 = vertical axis / horizontal axis * 100
    fpm(19) = nint(this%ang)                        ! fpm19 = rotation in degree
  end subroutine

  subroutine feast_option_get(this, Emid, r)
    !! return parameters needed for FEAST rci call
    class(feast_option), intent(in)  :: this
    complex(real64),     intent(out) :: Emid
      !! coordinate center
    real(real64),        intent(out) :: r
      !! horizontal radius

    Emid = this%mid
    r    = this%rad(1)
  end subroutine

  function check_info(info)
    !! FEAST returns status via info variable. its code is checked here. in case of errors/warnings messages are printed.
    integer, intent(in) :: info
      !! status code returned from FEAST
    logical             :: check_info
      !! does info show succesfull run?

    character(*), parameter :: ERROR_20X(200:202) = [ &
      "Problem with Emin, Emax or Emid, r",           &
      "Problem with size of subspace M0  ",           &
      "Problem with size of the system N "            &
    ]
    character(*), parameter :: WARNING_X(1:7) = [                                          &
      "No Eigenvalue found in the search interval                                       ", &
      "No Convergence (#iteration loops>fpm(4))                                         ", &
      "Size of the subspace M0 is too small (M0<=M)                                     ", &
      "Only the subspace has been returned using fpm(14)=1                              ", &
      "Only stochastic estimation of #eigenvalues returned fpm(14)=2                    ", &
      "FEAST converges but subspace is not bi-orthonormal                               ", &
      "The search for extreme eigenvalues has failed, search contour must be set by user"  &
    ]
    character(*), parameter :: ERROR_NEG_X(-3:-1) = [                         &
      "Internal error of the reduced eigenvalue solver                     ", &
      "Internal error of the inner system solver in FEAST Driver interfaces", &
      "Internal error conversion single/double                             "  &
    ]

    check_info = .false.
    if      ( info ==    0                    ) then
      check_info = .true.
    else if ( info >=  200                    ) then
      print *, '[FEAST ERROR] '//ERROR_20X(info)
    else if ( info >   100                    ) then
      print *, '[FEAST ERROR] Problem with ith value of the input FEAST parameter (i.e fpm(i)). i: '//int2str(info-100)
    else if ((info >=    1) .and. (info <=  7)) then
      print *, '[FEAST WARNING] '//WARNING_X(info)
    else if ((info >=   -3) .and. (info <= -1)) then
      print *, '[FEAST ERROR] '//ERROR_NEG_X(info)
    else if ( info <  -100                    ) then
      print *, '[FEAST ERROR] Problem with the ith argument of the FEAST interface. i: '//int2str(-info-100)
    else
      print *, '[FEAST] info code not catched.'
    end if
  end function

end module
