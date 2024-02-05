m4_include(macro.f90.inc)

module newton_m

  use error_m,         only: assert_failed, program_error
  use gmres_m,         only: gmres_options, gmres
  use ieee_arithmetic, only: ieee_value, ieee_is_finite, ieee_negative_inf, ieee_positive_inf
  use logging_m,       only: logging
  use matop_m,         only: single_matop_real
  use matrix_m,        only: matrix_real
  use util_m,          only: int2str, real2str

  implicit none

  private
  public newton1D, newton1D_opt
  public newton,   newton_opt

  type newton1D_opt
    !! options for 1D newton iteration

    real    :: atol
      !! absolute tolerance
    real    :: rtol
      !! relative tolerance
    real    :: xmin
      !! lower solution bound (can be -inf)
    real    :: xmax
      !! upper solution bound (can be +inf)
    real    :: dx_lim
      !! limit solution update
    integer :: min_it
      !! minimum number of iterations
    integer :: max_it
      !! maximum number of iterations

    logical                   :: log
      !! enable/disable logging
    character(:), allocatable :: msg
      !! message to print each iteration
  contains
    procedure :: init => newton1D_opt_init
  end type

  type newton_opt
    !! options for multidimensional newton iteration

    real, allocatable :: atol(:)
      !! absolute solution tolerance
    real, allocatable :: rtol(:)
      !! relative solution tolerance
    real, allocatable :: ftol(:)
      !! absolute residual tolerance
    real, allocatable :: dx_lim(:)
      !! limit solution update
    real, allocatable :: xmin(:), xmax(:)
      !! lower/upper solution bound (can be -inf/+inf)
    integer :: min_it
      !! minimum number of iterations
    integer           :: max_it
      !! maximum number of iterations

    logical                   :: it_solver
      !! direct solver or iterative solver
    logical                   :: log
      !! enable/disable logging
    logical                   :: error_if_not_converged
      !! stop program if convergence failed
    character(:), allocatable :: msg
      !! message to print each iteration
  contains
    procedure :: init => newton_opt_init
  end type

  abstract interface
    subroutine newton1D_fun(x, p, f, dfdx, dfdp)
      real,              intent(in)  :: x
        !! argument
      real,              intent(in)  :: p(:)
        !! parameters
      real,              intent(out) :: f
        !! output function value
      real,    optional, intent(out) :: dfdx
        !! optional output derivative of f wrt x
      real,    optional, intent(out) :: dfdp(:)
        !! optional output derivatives of f wrt p
    end subroutine

    subroutine newton_fun(x, p, f, dfdx, dfdx_prec, dfdp)
      import matrix_real
      real,                                  intent(in)  :: x(:)
        !! arguments
      real,                                  intent(in)  :: p(:)
        !! parameters
      real,               optional,          intent(out) :: f(:)
        !! output function values
      class(matrix_real), optional, pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
        !! must be factorized if preconditioner not used
      class(matrix_real), optional, pointer, intent(out) :: dfdx_prec
        !! optional output pointer to preconditioner jacobian of f wrt x
        !! must be factorized
      real,               optional,          intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p
    end subroutine
  end interface

contains

  subroutine newton1D(fun, p, opt, x0, x, dxdp)
    !! get root of 1D function by newton iteration with bisection stabilization
    procedure(newton1D_fun)         :: fun
      !! pointer to function
    real,               intent(in)  :: p(:)
      !! function parameters (can be empty array if not needed)
    type(newton1D_opt), intent(in)  :: opt
      !! iteration options
    real,               intent(in)  :: x0
      !! first guess for solution
    real,               intent(out) :: x
      !! output solution
    real,    optional,  intent(out) :: dxdp(:)
      !! optional output derivative of x wrt parameters

    ! local variables
    integer :: it
    real    :: err, err0, xmin, xmax, fmin, fmax
    real    :: f, dfdx, dfdp(size(p)), dx
    logical :: bounded

    ! init iteration params
    it   = 0
    err  = 0.5*huge(err)
    err0 = huge(err0)
    xmin = opt%xmin
    xmax = opt%xmax

    ! start value
    x = x0

    ! evaluate function at bounds
    bounded = .false.
    if (ieee_is_finite(xmin) .and. ieee_is_finite(xmax)) then
      bounded = .true.
      call fun(xmin, p, fmin)
      call fun(xmax, p, fmax)
      if ((fmin == 0) .and. (fmax == 0)) then
        x = 0.5 * (xmin + xmax)
      elseif (fmin == 0) then
        x = xmin
      elseif (fmax == 0) then
        x = xmax
      else
        if (sign(1.0, fmin) == sign(1.0, fmax)) call program_error("solution bounds are invalid, no sign change")
      end if
    end if

    ! newton iteration with bisection stabilization
    if (opt%log) m4_info(opt%msg  // " ITER      ABS_ERROR")
    do while ((it < opt%min_it) .or. ((err > opt%atol) .and. (err > abs(x) * opt%rtol)))
      it = it + 1

      ! check for maximum number of iterations
      if (it > opt%max_it) then
        m4_info("atol   = " // real2str(opt%atol))
        m4_info("rtol   = " // real2str(opt%rtol))
        m4_info("xmin   = " // real2str(xmin))
        m4_info("xmax   = " // real2str(xmax))
        m4_info("min_it = " // int2str(it))
        m4_info("max_it = " // int2str(it))
        m4_info("x      = " // real2str(x))
        m4_info("f      = " // real2str(f))
        m4_info("err    = " // real2str(err))
        call program_error("solution could not be found within maximum number of iterations")
      end if

      ! bisection
      if ((x < xmin) .or. (x > xmax) .or. (err0 <= err)) then
        ! at least one of the bounds is definitely finite, since it >= 2
        if      (.not. ieee_is_finite(xmin)) then
          x = 2 * x - xmax
        else if (.not. ieee_is_finite(xmax)) then
          x = 2 * x - xmin
        else
          x = 0.5 * (xmin + xmax)
        end if
      end if

      ! evaluate function
      call fun(x, p, f, dfdx=dfdx)

      ! calculate newton update and new error
      dx   = f / dfdx
      err0 = err
      err  = abs(dx)

      ! update bounds
      if (bounded) then
        if (sign(1.0, f) * sign(1.0, fmax) > 0) then
          xmax = x
          fmax = f
        else
          xmin = x
          fmin = f
        end if
      else
        if (sign(1.0, f) * sign(1.0, dfdx) > 0) then
          xmax = x
        else
          xmin = x
        end if
      end if

      ! limit update
      if (abs(dx) > opt%dx_lim) then
        dx = dx * opt%dx_lim / abs(dx)
      end if

      ! update solution
      x = x - dx

      if (opt%log) m4_info(opt%msg  // " " // int2str(it, "(I4)") // " "  // real2str(err, "(ES14.6E3)"))

      ! exit if close to solution
      if (ieee_is_finite(xmin) .and. ieee_is_finite(xmax)) then
        if (((xmax - xmin) < 0.5 * abs(xmax + xmin) * opt%rtol) .or. (0.5 * (xmax - xmin) < opt%atol)) then
          x = 0.5 * (xmax + xmin)
          exit
        end if
      end if
    end do

    ! calculate derivatives of solution wrt params by implicit differentiation
    if (present(dxdp)) then
      m4_assert(size(dxdp) == size(p))
      call fun(x, p, f, dfdx=dfdx, dfdp=dfdp)
      dxdp = - dfdp / dfdx
    end if
  end subroutine

  subroutine newton(fun, p, opt, x0, x, gmres_opt, dxdp)
    !! get root of multidimensional function by newton-raphson iteration
    procedure(newton_fun)                      :: fun
      !! pointer to function
    real,                          intent(in)  :: p(:)
      !! function parameters (can be empty array if not needed)
    type(newton_opt),              intent(in)  :: opt
      !! iteration options
    real,                          intent(in)  :: x0(:)
      !! first guess for solution
    real,                          intent(out) :: x(:)
      !! output solution
    type(gmres_options), optional, intent(in)  :: gmres_opt
      !! options for gmres iterative solver. only used when opt%it_sol is set.
    real,                optional, intent(out) :: dxdp(:,:)
      !! optional output derivatives of x wrt p

    ! local variables
    integer                         :: i, it
    real                            :: err, abs_err, dx0
    real,               allocatable :: f(:), dfdp(:,:), dx(:), xold(:)
    class(matrix_real), pointer     :: dfdx, dfdx_prec
    type(gmres_options)             :: gmres_opt_
    type(single_matop_real)         :: mulvec, precon

    m4_assert(size(x0      ) == size(x))
    m4_assert(size(opt%atol) == size(x))

    ! optional args
    gmres_opt_%atol = minval(opt%atol, dim = 1)
    gmres_opt_%rtol = minval(opt%rtol, dim = 1)
    if (present(gmres_opt)) gmres_opt_ = gmres_opt

    ! allocate memory
    allocate (f(   size(x)        ), source = huge(err))
    allocate (dfdp(size(x),size(p)), source = 0.0 )
    allocate (dx(  size(x)        ), source = 0.0 )
    allocate (xold(size(x)        ), source = 0.0 )
    nullify (dfdx, dfdx_prec)

    ! init iteration params
    it  = 0
    err = huge(err)

    ! start values
    x = x0
    where(x < opt%xmin) x = opt%xmin
    where(x > opt%xmax) x = opt%xmax

    ! newton-raphson iteration
    if (opt%log) m4_info(opt%msg // " ITER      REL_ERROR      ABS_ERROR       RESIDUAL")
    do while ((it < opt%min_it) .or. (any(err > opt%rtol) .and. any(abs(f) > opt%ftol)))
      it = it + 1

      ! check for maximum number of iterations
      if (it > opt%max_it) then
        m4_info("it      = " // int2str(it))
        m4_info("err     = " // real2str(err))
        m4_info("abs_err = " // real2str(abs_err))
        if (opt%error_if_not_converged) then
          m4_error("solution could not be found within maximum number of iterations")
        else
          m4_info("solution could not be found within maximum number of iterations")
          return
        end if
      end if

      ! determine update dx (either directly or by iterative solver)
      if (opt%it_solver) then
        ! initial solution
        dx = 0

        ! evaluate function: compute residuals, jacobian, preconditioner
        if (associated(dfdx     )) call dfdx%reset()
        if (associated(dfdx_prec)) call dfdx_prec%reset()
        call fun(x, p, f = f, dfdx = dfdx, dfdx_prec = dfdx_prec)
        call dfdx_prec%factorize()

        ! set matops + solve by gmres
        call mulvec%init(dfdx)
        call precon%init(dfdx_prec, inv = .true.)
        call gmres(f, mulvec, dx, opts = gmres_opt_, precon = precon)
      else
        ! evaluate function and solve
        if (associated(dfdx)) call dfdx%reset()
        call fun(x, p, f = f, dfdx = dfdx)
        call dfdx%factorize()
        call dfdx%solve_vec(f, dx)
      end if

      ! limit update
      dx0 = maxval(abs(dx) / opt%dx_lim, dim=1)
      if (dx0 > 1) dx = dx / dx0

      ! update solution
      xold = x
      x    = x - dx

      ! limit solution to bounds if necessary
      where(x < opt%xmin) x = opt%xmin
      where(x > opt%xmax) x = opt%xmax

      ! calculate new error
      err     = maxval(abs(x-xold) / (abs(xold) + opt%atol / opt%rtol))
      abs_err = maxval(abs(x-xold))

      if (opt%log) m4_info(opt%msg // " " // int2str(it, "(I4)") // " " // real2str(err, "(ES14.6E3)") // " " // real2str(abs_err, "(ES14.6E3)") // " " // real2str(maxval(abs(f)), "(ES14.6E3)"))
    end do

    ! calculate derivatives of solution wrt params by implicit differentiation
    if (present(dxdp)) then
      m4_assert(all(shape(dxdp) == [size(x), size(p)]))

      ! set dxdp by solving: dfdx * dxdp = -dfdp
      call fun(x, p, dfdp = dfdp)
      if (opt%it_solver) then
        dxdp = 0  ! init sols
        do i = 1, size(p)
          call gmres(-dfdp(:,i), mulvec, dxdp(:,i), opts = gmres_opt_, precon = precon)
        end do
      else
        call dfdx%solve_mat(-dfdp, dxdp)
      end if
    end if
  end subroutine

  subroutine newton1D_opt_init(this, atol, rtol, xmin, xmax, dx_lim, min_it, max_it, log, msg)
    !! initialize 1D newton iteration options
    class(newton1D_opt),    intent(out) :: this
    real,         optional, intent(in)  :: atol
      !! absolute tolerance (default: 1e-16)
    real,         optional, intent(in)  :: rtol
      !! relative tolerance (default: 1e-12)
    real,         optional, intent(in)  :: xmin
      !! lower solution bound (default: -inf)
    real,         optional, intent(in)  :: xmax
      !! upper solution bound (default: +inf)
    real,         optional, intent(in)  :: dx_lim
      !! update limit (default: huge)
    integer,      optional, intent(in)  :: min_it
      !! minimum number of iterations (default: 1)
    integer,      optional, intent(in)  :: max_it
      !! maximum number of iterations (default: huge)
    logical,      optional, intent(in)  :: log
      !! enable logging (default .false.)
    character(*), optional, intent(in)  :: msg
      !! message to print each iteration if logging is enabled (default: "")

    this%atol = 1e-16
    if (present(atol)) this%atol = atol
    this%rtol = 1e-12
    if (present(rtol)) this%rtol = rtol
    this%xmin = ieee_value(1.0, ieee_negative_inf)
    if (present(xmin)) this%xmin = xmin
    this%xmax = ieee_value(1.0, ieee_positive_inf)
    if (present(xmax)) this%xmax = xmax
    this%dx_lim = huge(1.0)
    if (present(dx_lim)) this%dx_lim = dx_lim
    this%min_it = 1
    if (present(min_it)) this%min_it = min_it
    this%max_it = huge(0)
    if (present(max_it)) this%max_it = max_it
    this%log = .false.
    if (present(log)) this%log = log
    this%msg = ""
    if (present(msg)) this%msg = msg
  end subroutine

  subroutine newton_opt_init(this, n, atol, rtol, ftol, dx_lim, xmin, xmax, min_it, max_it, it_solver, log, error_if_not_converged, msg)
    !! initialize multidimensional newton iteration options
    class(newton_opt),      intent(out) :: this
    integer,                intent(in)  :: n
      !! system size
    real,         optional, intent(in)  :: atol
      !! absolute solution tolerance (default: 1e-16)
    real,         optional, intent(in)  :: rtol
      !! relative solution tolerance (default: 1e-12)
    real,         optional, intent(in)  :: ftol
      !! absolute residual tolerance (default: 0.0)
    real,         optional, intent(in)  :: dx_lim
      !! limit for newton update (default: inf)
    real,         optional, intent(in)  :: xmin, xmax
      !! lower/upper solution bound (default: -inf/+inf)
    integer,      optional, intent(in)  :: min_it
      !! minimum number of iterations (default: 1)
    integer,      optional, intent(in)  :: max_it
      !! maximum number of iterations (default: huge)
    logical,      optional, intent(in)  :: it_solver
      !! direct solver or iterative solver (default: false == direct)
    logical,      optional, intent(in)  :: log
      !! enable logging (default .false.)
    logical,      optional, intent(in)  :: error_if_not_converged
      !! stop program if no convergence (default .true.)
    character(*), optional, intent(in)  :: msg
      !! message to print each iteration if logging is enabled (default: "")

    allocate (this%atol(n), this%rtol(n), this%ftol(n), this%dx_lim(n), this%xmin(n), this%xmax(n))

    this%atol = 1e-16
    if (present(atol)) this%atol = atol
    this%rtol = 1e-12
    if (present(rtol)) this%rtol = rtol
    this%ftol = 0.0
    if (present(ftol)) this%ftol = ftol
    this%dx_lim = ieee_value(1.0, ieee_positive_inf)
    if (present(dx_lim)) this%dx_lim = dx_lim
    this%xmin = ieee_value(1.0, ieee_negative_inf)
    if (present(xmin)) this%xmin = xmin
    this%xmax = ieee_value(1.0, ieee_positive_inf)
    if(present(xmax)) this%xmax = xmax
    this%min_it = 1
    if (present(min_it)) this%min_it = min_it
    this%max_it = huge(0)
    if (present(max_it)) this%max_it = max_it
    this%it_solver = .false.
    if (present(it_solver)) this%it_solver = it_solver
    this%log = .false.
    if (present(log)) this%log = log
    this%error_if_not_converged = .true.
    if (present(error_if_not_converged)) this%error_if_not_converged = error_if_not_converged
    this%msg = ""
    if (present(msg)) this%msg = msg
  end subroutine

end module
