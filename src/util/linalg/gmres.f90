m4_include(../macro.f90.inc)

module gmres_m

  use error_m,         only: assert_failed, program_error
  use iso_c_binding,   only: c_loc, c_f_pointer
  use iso_fortran_env, only: real64
  use matop_m,         only: matop_real, matop_cmplx, matop_c2r
  use mkl_ilu_m,       only: mkl_ilu
  use util_m,          only: int2str
  use vector_m,        only: vector_real

  implicit none

  private
  public gmres
  public gmres_options

  type gmres_options
    integer :: max_it    = 150
      !! maximum number of iterations - ipar(5)

    real    :: rtol      = 1e-12
      !! relative residual tolerance - dpar(1)
    real    :: atol      = 0.0
      !! absolute residual tolerance - dpar(2)

    logical :: print_msg = .false.
      !! print messages?
  contains
    procedure :: apply => gmres_options_apply
  end type

  interface gmres
    module procedure :: gmres_r
    module procedure :: gmres_c
  end interface

  interface
    ! - copied from /opt/exports/intel20/compilers_and_libraries_2020.1.217/linux/mkl/include/mkl_rci.fi
    ! - added kind information and changed formatting
    ! - including interfaces by '#include "mkl_rci.fi"' doesnt work for gfortran b.c. real data type is in interfaces
    !   defined as double precision which is real(kind=16) when using "-fdefault-real-8" compiler flag

    subroutine dfgmres(n, x, b, rci_request, ipar, dpar, tmp)
      import real64
      integer      :: n, rci_request, ipar(*)
      real(real64) :: x(*), b(*), dpar(*), tmp(*)
    end subroutine

    subroutine dfgmres_init(n, x, b, rci_request, ipar, dpar, tmp)
      import real64
      integer      :: n, rci_request, ipar(*)
      real(real64) :: x(*), b(*), dpar(*), tmp(*)
    end subroutine

    subroutine dfgmres_check(n, x, b, rci_request, ipar, dpar, tmp)
      import real64
      integer      :: n, rci_request, ipar(*)
      real(real64) :: x(*), b(*), dpar(*), tmp(*)
    end subroutine

    subroutine dfgmres_get(n, x, b, rci_request, ipar, dpar, tmp, itercount)
      import real64
      integer      :: n, rci_request, itercount, ipar(*)
      real(real64) :: x(*), b(*), dpar(*), tmp(*)
    end subroutine

  end interface

contains

  subroutine gmres_options_apply(this, ipar, dpar)
    !! apply this' values to ipar and dpar as well as setting useful default values.

    class(gmres_options), intent(in)    :: this
    integer,              intent(inout) :: ipar(128)
      !! integer parameters. see mkl doc
    real,                 intent(inout) :: dpar(128)
      !! real parameters. see mkl doc

    ipar( 5) = this%max_it
    ipar( 9) = 1              ! performs stopping test: dpar(5) <= dpar(4) aka residual norm smaller abs+tol thresholds
    ipar(10) = 0              ! user-defined stopping test is not requested
    ipar(12) = 1              ! zero-norm test is performed: dpar(7) <= dpar(8)

    dpar(1) = this%rtol
    dpar(2) = this%atol
  end subroutine

  subroutine gmres_c(b, mulvec, x, opts, precon, itercount, residual)
    !! wrapper around mkl's dfgmres for real(!) arrays.
    !!
    !! solves Ax=b and optionally, uses a preconditioner P.
    !!
    !! complex arrays are upscaled to real arrays of double the length.

    complex,             intent(in),    target           :: b(:)
      !! rhs
      !! reason for target: we will create real pointer pointing to its real and imag part.
    class(matop_cmplx),  intent(in),    target           :: mulvec
      !! matrix vector operation: x \mapsto A*x
      !! reason for target: matop_c2r will point to this matop.
    complex,             intent(inout), target           :: x(:)
      !! on input:  inital guess x_0
      !! on output: solution x
      !! reason for target: we will create real pointer pointing to its real and imag part.
    class(matop_cmplx),  intent(inout),    target, optional :: precon
      !! matrix vector operation: x \mapsto P*x
      !! reason for target: matop_c2r will point to this matop.
      !! inout: preconditioner will need to be factorized.
    type(gmres_options), intent(in),            optional :: opts
      !! parameters ipar, dpar.
    integer,             intent(out),           optional :: itercount
      !! iteration count
    real, allocatable,   intent(out),           optional :: residual(:)
      !! residuals over iteration

    real, pointer   :: x_ptr(:), b_ptr(:)
    type(matop_c2r) :: mulvec_c2r, precon_c2r

    m4_assert(size(x) == size(b))

    call c_f_pointer(c_loc(x), x_ptr, shape=shape(x)*2)
    call c_f_pointer(c_loc(b), b_ptr, shape=shape(b)*2)

    call mulvec_c2r%init(mulvec)

    if (present(precon)) then
      call precon_c2r%init(precon)
      call gmres(b_ptr, mulvec_c2r, x_ptr, opts=opts, precon=precon_c2r, itercount=itercount, residual=residual)

    else
      call gmres(b_ptr, mulvec_c2r, x_ptr, opts=opts,                    itercount=itercount, residual=residual)
    end if
  end subroutine

  subroutine gmres_r(b, mulvec, x, opts, precon, itercount, residual)
    !! wrapper around mkl's dfgmres.
    !!
    !! solves Ax=b and optionally, uses a preconditioner P.

    real,                intent(in)              :: b(:)
      !! rhs
    class(matop_real),   intent(in)              :: mulvec
      !! matrix vector operation: x \mapsto A*x
    real,                intent(inout)           :: x(:)
      !! on input:  inital guess x_0
      !! on output: solution x
    type(gmres_options), intent(in),    optional :: opts
      !! parameters ipar, dpar.
    class(matop_real),   intent(inout), optional :: precon
      !! matrix vector operation: x \mapsto P*x
      !! inout: preconditioner will need to be factorized.
    integer,             intent(out),   optional :: itercount
      !! iteration count
    real, allocatable,   intent(out),   optional :: residual(:)
      !! residuals over iteration

    integer                   :: rci_request, itercount_, n, ipar(128), ipar15_default
    real                      :: dpar(128)
    real, allocatable, target :: tmp(:)
    type(vector_real)         :: res

    n = size(x)
    m4_assert(n == size(b))

    ! init ipar, dpar by default values
    ipar15_default = min(150, n)  ! default value. needs to be set before allocating tmp?
    allocate (tmp(tmp_size(n, ipar15_default)))
    call dfgmres_init(n, x, b, rci_request, ipar, dpar, tmp)
    if (rci_request /= 0) call program_error('gmres init failed. rci_request: ' // int2str(rci_request))

    ! set values from options into ipar, dpar
    block
      type(gmres_options) :: opts_

      if (present(opts)) then
        call opts%apply( ipar, dpar)    ! use options supplied via argument
      else
        call opts_%apply(ipar, dpar)    ! use default options
      end if
    end block

    ! turn on/off preconditioning according to optional "precon" argument
    ipar(11) = merge(1, 0, present(precon))

    ! set ipar/dpar wrt preconditioner
    if (present(precon)) then
      select type (p => precon)
        class is (mkl_ilu)
          call p%apply(ipar, dpar)
      end select
    end if

    ! check ipar, dpar
    call dfgmres_check(n, x, b, rci_request, ipar, dpar, tmp)
    if (rci_request /= 0) call program_error('gmres check params failed. rci_request: ' // int2str(rci_request))

    ! (ipar(13) > 0) == (write solution into b)
    !     not implemented as b has intent(in)
    if (ipar(13) > 0) call program_error('writing solution into b is not supported.')

    ! allocated tmp array by default length. we might need to reallocate if it got changed.
    if (ipar(15) /= ipar15_default) then
      deallocate (tmp)
      allocate (tmp(tmp_size(n, ipar(15))), source=0.0)   ! init sets first elements to 0. we do it here manually
    end if

    ! factorize preconditioner
    if (present(precon)) then
      select type (p => precon)
        class is (mkl_ilu)
          block
            integer :: t(2)
            logical :: print_msg_
            real    :: cr

            print_msg_ = .false.
            if (present(opts)) print_msg_ = opts%print_msg

            if (print_msg_) call system_clock(count=t(1), count_rate=cr)

            call p%factorize(ipar, dpar)

            if (print_msg_) then
              call system_clock(count=t(2))
              print '(A, ES10.3, A)', " [gmres] precon factorized. time:", (t(2)-t(1))/cr, "s."
            end if
          end block
      end select
    end if

    itercount_ = 0
    if (present(residual)) call res%init(0, c=30)

    ! gmres solution process
    LOOP: do
      ! Compute the solution by RCI (P)FGMRES solver
      call dfgmres(n, x, b, rci_request, ipar, dpar, tmp)

      ! extract current residual
      if (present(residual)) call res%push(dpar(5))

      select case (rci_request)
        ! If RCI_REQUEST==0, then the solution was found with the required precision
        case (0)
          exit

        ! If RCI_REQUEST==1, then compute the vector A*TMP(IPAR(22))
        case (1)
          associate (v_in  => tmp(ipar(22):ipar(22)+(n-1)), &
            &        v_out => tmp(ipar(23):ipar(23)+(n-1))  )

            call mulvec%exec(v_in, v_out)
          end associate
          itercount_ = itercount_ + 1

        ! stopping test
        case (2)
          call program_error('user-defined stopping test not supplied but requested!!')

        ! apply preconditioner
        case (3)
          if (.not. present(precon)) call program_error('preconditioner rountine not supplied but requested!!')
          associate (bb => tmp(ipar(22):ipar(22)+(n-1)), &
            &        xx => tmp(ipar(23):ipar(23)+(n-1))  )

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
    call dfgmres_get(n, x, b, rci_request, ipar, dpar, tmp, itercount_)
    if (rci_request /= 0) call program_error('gmres get failed. rci_request: ' // int2str(rci_request))

    ! set output variables
    if (present(itercount)) itercount = itercount_
    if (present(residual )) then
      residual = res%to_array()
      if (res%n > 1) residual = residual(2:)
    end if
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
