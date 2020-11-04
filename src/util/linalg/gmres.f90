#include "../macro.f90.inc"

module gmres_m
  ! fixme gmres for complex arguemnts works by scaling the complex arrays up to real arrays of double the length.
  !       we should find a routine which works directly on the complex arrays. should be cleaner/faster.

  use error_m
  use matop_m, only: matop_real, matop_cmplx
  use util_m,  only: int2str

  implicit none

  private
  public gmres
  public gmres_options

  type gmres_options
    integer :: ipar(128)
      !! mkl's gmres ipar
      !!
      !! ipar( 4): current iteration number
      !! ipar( 5): maximal iteration number
      !! ipar( 8): stopping test wrt iteration number
      !!             == 0: doesnt perform  stopping test
      !!             /= 0:        performs stopping test: ipar(4) <= ipar(5)                  *default*
      !! ipar( 9): stopping test wrt residual
      !!             == 0: doesnt perform  stopping test                                      *default*
      !!             /= 0:        performs stopping test: dpar(5) <= dpar(4)
      !! ipar(10): user defined stopping test
      !!             == 0: user-defined stopping test is not requested
      !!             /= 0: user-defined stopping test is     requested                        *default*
      !! ipar(11): preconditioning
      !!             == 0: preconditioning is not used                                        *default*
      !!             /= 0: preconditioning is     used
      !! ipar(12): zero-norm test of vector
      !!             == 0: user-defined zero-norm test is requested                           *default*
      !!             /= 0: automatic    zero-norm test is performed: dpar(7) <= dpar(8)
      !! ipar(13): solution vector
      !!             <  0: x and b not updated
      !!             == 0: writes solution into argument x                                    *default*
      !!             >  0: writes solution into argument b
      !! ipar(15): non-restarted gmres iterations

    real    :: dpar(128)
      !! mkl's gmres dpar
      !!
      !! dpar(1): relative residual tolerance. default: 1e-6
      !! dpar(2): absolute residual tolerance. default: 0.0
      !! dpar(3): work variable = 2-norm of initial residual.
      !! dpar(4): work variable = dpar(1)*dpar(3)+dpar(2)
      !! dpar(5): work variable = 2-norm of current residual
      !! dpar(7): work variable = 2-norm of currently generated vector
      !! dpar(8): tolerance for zero-norm test

  contains
    generic   :: init   => gmres_options_init_c, gmres_options_init_r
    generic   :: check  => gmres_options_check_c, gmres_options_check_r
    procedure :: print  => gmres_options_print
    procedure :: pprint => gmres_options_pprint
    procedure, private :: gmres_options_init_c, gmres_options_init_r
    procedure, private :: gmres_options_check_c, gmres_options_check_r
  end type

  interface
    ! - copied from /opt/exports/intel20/compilers_and_libraries_2020.1.217/linux/mkl/include/mkl_rci.fi
    ! - added kind information and changed formatting
    ! - including interfaces by '#include "mkl_rci.fi"' doesnt work for gfortran b.c. real data type is in interfaces
    !   defined as double precision which is real(kind=16) when using "-fdefault-real-8" compiler flag

    subroutine dfgmres(n, x, b, rci_request, ipar, dpar, tmp)
      integer(kind=8) :: n, rci_request, ipar(*)
      real(   kind=8) :: x(*), b(*), tmp(*), dpar(*)
    end subroutine

    subroutine dfgmres_init(n, x, b, rci_request, ipar, dpar, tmp)
      integer(kind=8) :: n, rci_request, ipar(*)
      real(   kind=8) :: x(*), b(*), tmp(*), dpar(*)
    end subroutine

    subroutine dfgmres_check(n, x, b, rci_request, ipar, dpar, tmp)
      integer(kind=8) :: n, rci_request, ipar(*)
      real(   kind=8) :: x(*), b(*), tmp(*), dpar(*)
    end subroutine

    subroutine dfgmres_get(n, x, b, rci_request, ipar, dpar, tmp, itercount)
      integer(kind=8) :: n, rci_request, itercount, ipar(*)
      real(   kind=8) :: x(*), b(*), tmp(*), dpar(*)
    end subroutine

  end interface

  interface gmres
    module procedure :: gmres_r
    module procedure :: gmres_c
  end interface

  ! fixme should we make this type more easily available? inside of matop.f90?
  ! fixme because this matop has worktime variables which are altered by exec1 the intents of exec1 and even the matops are alle intent(inout).
  !       e.g. see intent in arnoldi.f90 matop A which could be intent(in).
  !       e.g. all exec1, exec2: matops could be: class(matop), intent(in!!!) :: this
  type, extends(matop_real) :: matop_c2r
    !! matop_cmplx should use gmres for real(!) residuals of double length
    !!
    !! idea:  - init
    !!          - save complex matop
    !!          - alloc worktime variables
    !!        - mat_vec operation (exec1)
    !!          - gets real of double length as input
    !!          - convert input: 1st half is real, 2nd half the imaginary part of the actual vector
    !!          - compute mat_vec with saved complex matop
    !!          - convert result: complex vector to [real, imag] as real array of double length
    !!          - output: double length real array

    private
    class(matop_cmplx), pointer     :: mop_c => null()
      !! complex matop
    complex,            allocatable :: work_x(:), work_y(:)
      !! worktime variables
  contains
    procedure :: matop_c2r_init
    generic   :: init  => matop_c2r_init
    procedure :: exec1 => matop_c2r_exec1
    procedure :: exec2 => matop_c2r_exec2
  end type

contains

  subroutine matop_c2r_init(this, mop_c)
    class(matop_c2r),   intent(out)        :: this
    class(matop_cmplx), intent(in), target :: mop_c

    ! set pointer to complex matop
    this%mop_c => mop_c

    ! allocate worktime variables
    allocate (this%work_x(this%mop_c%ncols), this%work_y(this%mop_c%nrows))
  end subroutine

  subroutine matop_c2r_exec1(this, x, y)
    !! compute mul_vec of input and saved mat_op.
    !!
    !! 1) given double length real array (realpart, imagpart)
    !! 2) convert to complex array
    !! 3) mul_vec with saved mat_op in complex space
    !! 4) convert result to double length real array (realpart, imagpart)

    class(matop_c2r), intent(inout) :: this
    real,             intent(in)    :: x(:)
      !! input vector ([real(x_c), aimag(x_c)])
    real,             intent(out)   :: y(:)
      !! output vector ([real(y_c), aimag(y_c)])

    this%work_x = cmplx(x(:size(x)/2), x(size(x)/2+1:))
    call this%mop_c%exec(this%work_x, this%work_y)
    y = [real(this%work_y), aimag(this%work_y)]
  end subroutine

  subroutine matop_c2r_exec2(this, x, y)
    !! routine is not used in gmres which is why we dont implement it here.
    !! need to create this placeholder anyways as it is marked as deferred in parent type.

    class(matop_c2r), intent(inout) :: this
    real,             intent(in)    :: x(:,:)
    real,             intent(out)   :: y(:,:)

    IGNORE(this)
    IGNORE(x)
    y = 0 ! to remove warnings for "intent(out) variable not assigned value"
    call program_error("not implemented.")
  end subroutine

  subroutine gmres_options_init_c(this, x, b)
    !! inits parameters ipar, dpar using mkl's dfgmres_init.
    !!
    !! wrapper around init routine for real arguments.
    !! converts complex array to real array of double length: [realpart, imagpart].

    class(gmres_options), intent(out) :: this
    complex,              intent(in)  :: x(:)
      !! initial guess
    complex,              intent(in)  :: b(:)
      !! rhs

    call this%init([real(x), aimag(x)], [real(b), aimag(b)])
  end subroutine

  subroutine gmres_options_init_r(this, x, b)
    !! inits parameters ipar, dpar using mkl's dfgmres_init.

    class(gmres_options), intent(out) :: this
    real,                 intent(in)  :: x(:)
      !! initial guess
    real,                 intent(in)  :: b(:)
      !! rhs

    integer           :: n, ipar15, rci_request
    real, allocatable :: tmp_dummy(:)
      !! dummy variable. init sets first couple elements to 0. will be done manually in gmres call...

    n = size(x)
    ASSERT(n == size(b))

    ! FIXME why do we need to set this beforehand? maybe whole calling structure should be differnt?
    !   reason: - dfgmres_init needs tmp to be of right size.
    !           - correct size can be computed by ipar(15) which is not set before calling init???
    !           - thus, set it here manually
    ipar15 = min(150, n)      ! default value according to mkl's documentation

    allocate (tmp_dummy(tmp_size(n, ipar15)))

    ! gmres init call
    call dfgmres_init(n, x, b, rci_request, this%ipar, this%dpar, tmp_dummy)
    if (rci_request /= 0) call program_error('gmres init failed. rci_request: ' // int2str(rci_request))
  end subroutine

  subroutine gmres_options_check_c(this, x, b)
    !! checks parameters ipar, dpar using mkl's dfgmres_check.
    !!
    !! wrapper around check routine for real arguments.
    !! converts complex array to real array of double length: [realpart, imagpart].

    class(gmres_options), intent(inout) :: this
    complex,              intent(in)    :: x(:)
      !! initial guess
    complex,              intent(in)    :: b(:)
      !! rhs

    call this%check([real(x), aimag(x)], [real(b), aimag(b)])
  end subroutine

  subroutine gmres_options_check_r(this, x, b)
    !! checks parameters ipar, dpar using mkl's dfgmres_check.

    class(gmres_options), intent(inout) :: this
    real,                 intent(in)    :: x(:)
      !! initial guess
    real,                 intent(in)    :: b(:)
      !! rhs

    integer           :: n, rci_request
    real, allocatable :: tmp_dummy(:)
      !! dummy variable. init sets first couple elements to 0. will be done manually in gmres call...

    n = size(x)
    ASSERT(n == size(b))

    ! init vars: dfgmres_init sets first elements to 0. we do that here manually.
    allocate (tmp_dummy(tmp_size(n, this%ipar(15))), source=0.0)

    ! dfgmres_init sets it to 0
    rci_request = 0

    ! check params
    call dfgmres_check(n, x, b, rci_request, this%ipar, this%dpar, tmp_dummy)
    if (rci_request /= 0) call program_error('gmres check params failed. rci_request: ' // int2str(rci_request))
  end subroutine

  subroutine gmres_options_print(this)
    !! prints parameters but only those that correspond to user-changeable data, e.g. ipar(16:128) are work variables.

    class(gmres_options), intent(in) :: this

    integer :: i

    do i = 1, 15
      print '(A, I5, I5     )', ' i, ipar(i)', i, this%ipar(i)
    end do

    do i = 1, 8
      print '(A, I5, ES24.16)', ' i, dpar(i)', i, this%dpar(i)
    end do
  end subroutine

  subroutine gmres_options_pprint(this)
    !! **pretty** prints parameters but only those that correspond to user-changeable data, e.g. ipar(16:128) are
    !! work variables.

    class(gmres_options), intent(in) :: this

    print *, 'ipar( 1) - system size:                                  ',  this%ipar( 1)
    print *, 'ipar( 3) - rci_request:                                  ',  this%ipar( 3)
    print *, 'ipar( 4) - current iteration number:                     ',  this%ipar( 4)
    print *, 'ipar( 5) - maximal iteration number:                     ',  this%ipar( 5)
    print *, 'ipar( 8) - performing stopping test wrt iteration number?', (this%ipar( 8) /= 0)
    print *, 'ipar( 9) - performing stopping test wrt residual?        ', (this%ipar( 9) /= 0)
    print *, 'ipar(10) - requesting user-defined stopping test?        ', (this%ipar(10) /= 0)
    print *, 'ipar(11) - preconditioning used?                         ', (this%ipar(11) /= 0)
    print *, 'ipar(12) - automatic zero-norm test?                     ', (this%ipar(12) /= 0)
    if      (this%ipar(13) < 0) then
      print *, 'ipar(13) - x and b not updated'
    else if (this%ipar(13) == 0) then
      print *, 'ipar(13) - writing solution into x'
    else
      print *, 'ipar(13) - writing solution into b'
    end if
    print *, 'ipar(15) - number of non-restarted iterations:           ', this%ipar(15)

    print *, 'dpar(1) - relative residual tolerance:       ', this%dpar(1)
    print *, 'dpar(2) - absolute residual tolerance:       ', this%dpar(2)
    print *, 'dpar(3) - initial residual:                  ', this%dpar(3)
    print *, 'dpar(4) - total residual tolerance:          ', this%dpar(4)
    print *, 'dpar(5) - current residual:                  ', this%dpar(5)
    print *, 'dpar(7) - norm of currently generated vector:', this%dpar(7)
    print *, 'dpar(8) - zero-norm tolerance:               ', this%dpar(8)
  end subroutine

  subroutine gmres_c(opts, b, mulvec, x, precon)
    !! wrapper around mkl's dfgmres for real(!) arrays.
    !!
    !! solves Ax=b and optionally, uses a preconditioner P.
    !!
    !! complex arrays are upscaled to real arrays of double the length.

    type(gmres_options), intent(inout)                :: opts
      !! parameters ipar, dpar.
      !! inout: gmres calls change some of them, e.g. work variables and iteration count.
    complex,             intent(in)                   :: b(:)
      !! rhs
    class(matop_cmplx),  intent(in),           target :: mulvec
      !! matrix vector operation: x \mapsto A*x
      !! reason for class: either use already defined single_matop or derive your own.
    complex,             intent(inout)                :: x(:)
      !! on input:  inital guess x_0
      !! on output: solution x
    class(matop_cmplx),  intent(in), optional, target :: precon
      !! matrix vector operation: x \mapsto P*x
      !! reason for class: either use already defined single_matop or derive your own, e.g. derive a matop for mkl's ILUt.

    real, allocatable :: b_c2r(:), x_c2r(:)
    type(matop_c2r)   :: mulvec_c2r, precon_c2r

    ASSERT(size(x) == size(b))

    call mulvec_c2r%init(mulvec)

    b_c2r = [real(b), aimag(b)]
    x_c2r = [real(x), aimag(x)]

    if (present(precon)) then
      call precon_c2r%init(precon)
      call gmres(opts, b_c2r, mulvec_c2r, x_c2r, precon=precon_c2r)

    else
      call gmres(opts, b_c2r, mulvec_c2r, x_c2r)
    end if

    x = cmplx(x_c2r(:size(x)), x_c2r(size(x)+1:))
  end subroutine

  subroutine gmres_r(opts, b, mulvec, x, precon)
    !! wrapper around mkl's dfgmres.
    !!
    !! solves Ax=b and optionally, uses a preconditioner P.

    type(gmres_options), intent(inout)           :: opts
      !! parameters ipar, dpar.
      !! inout: gmres calls change some of them, e.g. work variables and iteration count.
    real,                intent(in)              :: b(:)
      !! rhs
    class(matop_real),   intent(inout)           :: mulvec
      !! matrix vector operation: x \mapsto A*x
      !! reason for class: either use already defined single_matop or derive your own.
    real,                intent(inout)           :: x(:)
      !! on input:  inital guess x_0
      !! on output: solution x
    class(matop_real),   intent(inout), optional :: precon
      !! matrix vector operation: x \mapsto P*x
      !! reason for class: either use already defined single_matop or derive your own, e.g. derive a matop for mkl's ILUt.

    integer           :: rci_request, itercount, n
    real, allocatable :: tmp(:)

    n = size(x)
    ASSERT(n == size(b))

    ! (ipar(13) > 0) == (write solution into b)
    !     not implemented as b has intent(in)
    if (opts%ipar(13) > 0) call program_error('writing solution into b is not supported.')

    ! check precon flag
    if (present(precon)) then
      if (opts%ipar(11) == 0) print *, '[gmres] WARN! Preconditioner supplied to gmres but not turned on according to options?!'
    else
      if (opts%ipar(11) /= 0) call program_error('precon requested by ipar but not supplied!')
    end if

    ! init vars: dfgmres_init sets first elements to 0. we do that here manually.
    allocate (tmp(tmp_size(n, opts%ipar(15))), source=0.0)

    ! dfgmres_init, dfgmres_check return rci_req==0 if no errors occured
    rci_request = 0

    itercount = 0

    ! gmres solution process
    LOOP: do
      ! Compute the solution by RCI (P)FGMRES solver
      call dfgmres(n, x, b, rci_request, opts%ipar, opts%dpar, tmp)

      select case (rci_request)
        ! If RCI_REQUEST==0, then the solution was found with the required precision
        case (0)
          exit

        ! If RCI_REQUEST==1, then compute the vector A*TMP(IPAR(22))
        case (1)
          associate (v_in  => tmp(opts%ipar(22):opts%ipar(22)+(n-1)), &
            &        v_out => tmp(opts%ipar(23):opts%ipar(23)+(n-1))  )

            call mulvec%exec(v_in, v_out)
          end associate
          itercount = itercount + 1

        ! stopping test
        case (2)
          call program_error('user-defined stopping test not supplied but requested!!')

        ! apply preconditioner
        case (3)
          if (.not. present(precon)) call program_error('preconditioner rountine not supplied but requested!!')
          associate (bb => tmp(opts%ipar(22):opts%ipar(22)+(n-1)), &
            &        xx => tmp(opts%ipar(23):opts%ipar(23)+(n-1))  )

            call precon%exec(bb, xx)
          end associate

        ! vector zero test
        case (4)
          call program_error('user-defined zero vector test not supplied but requested!!')

        ! If RCI_REQUEST=anything else, then DFGMRES subroutine failed
        case default
          print *, '[gmres] gmres loop failed. rci_request: ' // int2str(rci_request)
          exit LOOP
      end select
    end do LOOP

    ! get solution
    call dfgmres_get(n, x, b, rci_request, opts%ipar, opts%dpar, tmp, itercount)
    if (rci_request /= 0) call program_error('gmres get failed. rci_request: ' // int2str(rci_request))
  end subroutine

  integer function tmp_size(n, ipar15)
    !! mkl's gmres needs a work variable "tmp". this function computes its size according to mkl documentation.

    integer, intent(in) :: n
      !! system size
    integer, intent(in) :: ipar15
      !! ipar(15). see mkl documentation

    tmp_size = n*(2*ipar15+1) + (ipar15*(ipar15+9))/2 + 1
  end function

end module
