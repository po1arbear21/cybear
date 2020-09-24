#include "../macro.f90.inc"

module gmres_m

  use error_m
  use matop_m, only: matop_real
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
    procedure :: init  => gmres_options_init
    procedure :: check => gmres_options_check
    procedure :: print => gmres_options_print
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

contains

  subroutine gmres_options_init(this, x, b)
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

  subroutine gmres_options_check(this, x, b)
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

  subroutine gmres(opts, b, mulvec, x, precon)
    !! wrapper around mkl's dfgmres.
    !!
    !! solves Ax=b and optionally, uses a preconditioner P.

    type(gmres_options), intent(inout)        :: opts
      !! parameters ipar, dpar.
      !! inout: gmres calls change some of them, e.g. work variables and iteration count.
    real,                intent(in)           :: b(:)
      !! rhs
    class(matop_real),   intent(in)           :: mulvec
      !! matrix vector operation: x \mapsto A*x
      !! reason for class: either use already defined single_matop or derive your own.
    real,                intent(inout)        :: x(:)
      !! on input:  inital guess x_0
      !! on output: solution x
    class(matop_real),   intent(in), optional :: precon
      !! matrix vector operation: x \mapsto P*x
      !! reason for class: either use already defined single_matop or derive your own, e.g. derive a matop for mkl's ILUt.

    integer           :: rci_request, itercount, n
    real, allocatable :: tmp(:)

    n = size(x)
    ASSERT(n == size(b))

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
