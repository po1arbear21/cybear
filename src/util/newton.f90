module newton_m
  use error_m
  use ieee_arithmetic
  use matrix_m

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
    integer           :: max_it
      !! maximum number of iterations

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
    subroutine newton1D_fun(x, p, f, ipar, dfdx, dfdp)
      real,              intent(in)  :: x
        !! argument
      real,              intent(in)  :: p(:)
        !! parameters
      real,              intent(out) :: f
        !! output function value
      integer, optional, intent(in)  :: ipar(:)
        !! integer parameters
      real,    optional, intent(out) :: dfdx
        !! optional output derivative of f wrt x
      real,    optional, intent(out) :: dfdp(:)
        !! optional output derivatives of f wrt p
    end subroutine

    subroutine newton_fun(x, p, f, dfdx, dfdp)
      import matrix_real
      real,                        intent(in)  :: x(:)
        !! arguments
      real,                        intent(in)  :: p(:)
        !! parameters
      real,                        intent(out) :: f(:)
        !! output function values
      class(matrix_real), pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
      real, optional,              intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p
    end subroutine
  end interface

contains

  subroutine newton1D(fun, p, opt, x0, x, ipar, dxdp)
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
    integer, optional,  intent(in)  :: ipar(:)
      !! integer function parameters
    real,    optional,  intent(out) :: dxdp(:)
      !! optional output derivative of x wrt parameters

    ! local variables
    integer :: i, it
    real    :: err, err0, xmin, xmax, fmin, fmax
    real    :: f, dfdx, dfdp(size(p)), dx
    logical :: bounded

    ! init iteration params
    it   = 0
    err  = 1e99
    err0 = 2e99
    xmin = opt%xmin
    xmax = opt%xmax

    ! start value
    x = x0

    ! evaluate function at bounds
    bounded = .false.
    if (ieee_is_finite(xmin) .and. ieee_is_finite(xmax)) then
      bounded = .true.
      call fun(xmin, p, fmin, ipar=ipar)
      call fun(xmax, p, fmax, ipar=ipar)
      if (sign(1.0, fmin) == sign(1.0, fmax)) then
        call program_error("solution bounds are invalid, no sign change")
      end if
    end if

    ! newton iteration with bisection stabilization
    do while ((it < 1) .or. ((err > opt%atol) .and. (err > abs(x) * opt%rtol)))
      it = it + 1

      ! check for maximum number of iterations
      if (it > opt%max_it) then
        print "(A,ES24.16)", "atol   = ", opt%atol
        print "(A,ES24.16)", "rtol   = ", opt%rtol
        print "(A,ES24.16)", "xmin   = ", xmin
        print "(A,ES24.16)", "xmax   = ", xmax
        print "(A,I0)",      "max_it = ", opt%max_it
        print "(A,ES24.16)", "x      = ", x
        print "(A,ES24.16)", "f      = ", f
        print "(A,ES24.16)", "err    = ", err
        call program_error("solution could not be found within maximum number of iterations")
      end if

      ! bisection
      if ((x < xmin) .or. (x > xmax) .or. (err0 <= err)) then
        ! at least one of the bounds is definitely finite, since it >= 2
        if (.not. ieee_is_finite(xmin)) then
          x = 2 * x - xmax
        elseif (.not. ieee_is_finite(xmax)) then
          x = 2 * x - xmin
        else
          x = 0.5 * (xmin + xmax)
        end if
      end if

      ! evaluate function
      call fun(x, p, f, dfdx=dfdx, ipar=ipar)

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

      ! update solution
      x = x - dx

      if (opt%log) print "(A,I0,ES24.16)", opt%msg, it, err

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
      call fun(x, p, f, dfdx=dfdx, dfdp=dfdp, ipar=ipar)
      do i = 1, size(p)
        dxdp(i) = - dfdp(i) / dfdx
      end do
    end if
  end subroutine

  subroutine newton(fun, p, opt, x0, x, dxdp)
    !! get root of multidimensional function by newton-raphson iteration
    procedure(newton_fun)         :: fun
      !! pointer to function
    real,             intent(in)  :: p(:)
      !! function parameters (can be empty array if not needed)
    type(newton_opt), intent(in)  :: opt
      !! iteration options
    real,             intent(in)  :: x0(:)
      !! first guess for solution
    real,             intent(out) :: x(:)
      !! output solution
    real, optional,   intent(out) :: dxdp(:,:)
      !! optional output derivatives of x wrt p

    ! local variables
    integer                         :: i, it
    real                            :: err, abs_err
    real,               allocatable :: f(:), dfdp(:,:), dx(:)
    class(matrix_real), pointer     :: dfdx

    ! allocate memory
    allocate (f(   size(x)        ), source = 1e99)
    allocate (dfdp(size(x),size(p)), source = 0.0)
    allocate (dx(  size(x)        ), source = 0.0)
    nullify (dfdx)

    ! init iteration params
    it  = 0
    err = 1e99

    ! start values
    x = x0

    ! newton-raphson iteration
    do while (any(err > opt%rtol) .and. any(abs(f) > opt%ftol))
      it = it + 1

      ! check for maximum number of iterations
      if (it > opt%max_it) then
        print "(A,I0)",      "it      = ", it
        print "(A,ES24.16)", "err     = ", err
        print "(A,ES24.16)", "abs_err = ", abs_err
        if (opt%error_if_not_converged) then
          call program_error("solution could not be found within maximum number of iterations")
        else
          print *, "solution could not be found within maximum number of iterations"
          return
        end if
      end if

      ! evaluate function
      call fun(x, p, f, dfdx)

      ! factorize and solve
      call dfdx%factorize()
      call dfdx%solve_vec(f, dx)
      call dfdx%destruct()

      ! calculate new error
      err     = maxval(abs(dx) / (abs(x) + opt%atol / opt%rtol))
      abs_err = maxval(abs(dx))

      ! limit update
      do i = 1, size(dx)
        if (abs(dx(i)) > opt%dx_lim(i)) then
          dx = dx * opt%dx_lim(i) / abs(dx(i))
        end if
      end do

      ! update solution
      x = x - dx

      if (opt%log) then
        print "(A,I0,3ES25.16)", opt%msg, it, err, abs_err, maxval(abs(f))
      end if
    end do

    ! calculate derivatives of solution wrt params by implicit differentiation
    if (present(dxdp)) then
      call fun(x, p, f, dfdx, dfdp = dfdp)

      call dfdx%factorize()
      call dfdx%solve_mat(-dfdp, dxdp)
    end if
  end subroutine

  subroutine newton1D_opt_init(this, atol, rtol, xmin, xmax, max_it, log, msg)
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
    this%max_it = huge(0)
    if (present(max_it)) this%max_it = max_it
    this%log = .false.
    if (present(log)) this%log = log
    this%msg = ""
    if (present(msg)) this%msg = msg
  end subroutine

  subroutine newton_opt_init(this, n, atol, rtol, ftol, dx_lim, max_it, log, error_if_not_converged, msg)
    !! initialize multidimensional newton iteration options
    class(newton_opt),      intent(out) :: this
    integer,                intent(in)  :: n
      !! system size
    real,         optional, intent(in)  :: atol(:)
      !! absolute solution tolerance (default: 1e-16)
    real,         optional, intent(in)  :: rtol(:)
      !! relative solution tolerance (default: 1e-12)
    real,         optional, intent(in)  :: ftol(:)
      !! absolute residual tolerance (default: 0.0)
    real,         optional, intent(in)  :: dx_lim(:)
      !! limit for newton update (default: inf)
    integer,      optional, intent(in)  :: max_it
      !! maximum number of iterations (default: huge)
    logical,      optional, intent(in)  :: log
      !! enable logging (default .false.)
    logical,      optional, intent(in)  :: error_if_not_converged
      !! stop program if no convergence (default .true.)
    character(*), optional, intent(in)  :: msg
      !! message to print each iteration if logging is enabled (default: "")

    allocate (this%atol(n), this%rtol(n), this%ftol(n), this%dx_lim(n))

    this%atol = 1e-16
    if (present(atol)) this%atol = atol
    this%rtol = 1e-12
    if (present(rtol)) this%rtol = rtol
    this%ftol = 0.0
    if (present(ftol)) this%ftol = ftol
    this%dx_lim = ieee_value(1.0, ieee_positive_inf)
    if (present(dx_lim)) this%dx_lim = dx_lim
    this%max_it = huge(0)
    if (present(max_it)) this%max_it = max_it
    this%log = .false.
    if (present(log)) this%log = log
    this%error_if_not_converged = .true.
    if (present(error_if_not_converged)) this%error_if_not_converged = error_if_not_converged
    this%msg = ""
    if (present(msg)) this%msg = msg
  end subroutine

end module
