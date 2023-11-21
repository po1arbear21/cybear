m4_include(../macro.f90.inc)

module ode_m
  use error_m
  implicit none

  private
  public ode_fun
  public ode_options
  public ode_result
  public ode_solve

  type ode_options
    !! options for ode solver
    real, allocatable :: atol(:)
      !! absolute tolerance for each component
    real, allocatable :: rtol(:)
      !! relative tolerance for each component
    real              :: newton_atol
      !! absolute tolerance for newton iteration
    real              :: newton_rtol
      !! relative tolerance for newton iteration
    integer           :: newton_max_it
      !! maximum number of allowed newton iterations
    integer           :: max_rejected
      !! maximum number of allowed rejected steps
    real              :: min_rx
      !! minimum relative stepsize
    real              :: initial_rx
      !! initial relative stepsize
  contains
    procedure :: init => ode_options_init ! initialize
  end type

  type ode_result
    !! results returned by ode solver

    real, allocatable :: Usmp(:,:)
      !! return values at sample points (nU x size(xs))
    real, allocatable :: dUsmpdU0(:,:,:)
      !! return derivatives of sample values wrt U0 (nU x nU x size(xs))
    real, allocatable :: dUsmpdP(:,:,:)
      !! return derivatives of sample values wrt P (nU x nP x size(xs))

    integer :: nsteps
      !! return number of steps taken
  end type

  interface
    subroutine ode_fun(x, U, P, status, f, dfdU, dfdP)
      real,           intent(in)  :: x
        !! x coordinate
      real,           intent(in)  :: U(:)
        !! state (nU)
      real,           intent(in)  :: P(:)
        !! parameters (nP)
      logical,        intent(out) :: status
        !! success/fail of dU/dx calculation
      real, optional, intent(out) :: f(:)
        !! output dU/dx (nU)
      real, optional, intent(out) :: dfdU(:,:)
        !! output derivatives of f wrt U (nU,nU)
      real, optional, intent(out) :: dfdP(:,:)
        !! output derivatives of f wrt P (nU,nP)
    end subroutine

    subroutine ode_kernel(fun, xold, Uold, x, dxk, Uk, dUkdQ, fk, dfkdUk, dfkdP, polyk, P, opt, &
      &                                       dxn, Un, dUndQ, fn, dfndUn, dfndP, polyn, dpolyndQ, err, status)
      import ode_options, ode_fun
      procedure(ode_fun)               :: fun
        !! function to integrate
      real,              intent(in)    :: xold
        !! initial x position of previous step
      real,              intent(in)    :: Uold(:)
        !! state at beginning of previous step
      real,              intent(in)    :: x
        !! initial x position of step
      real,              intent(in)    :: dxk
        !! stepsize
      real,              intent(in)    :: Uk(:)
        !! state at beginning of step
      real,              intent(in)    :: dUkdQ(:,:)
        !! derivatives of Uk wrt Q = [U0; P]
      real,              intent(in)    :: fk(:)
        !! dUk/dx
      real,              intent(in)    :: dfkdUk(:,:)
        !! derivatives of fk wrt Uk
      real,              intent(in)    :: dfkdP(:,:)
        !! derivatives of fk wrt P
      real,              intent(in)    :: polyk(:,:)
        !! interpolation polynomial (from previous step)
      real,              intent(in)    :: P(:)
        !! parameters
      type(ode_options), intent(in)    :: opt
        !! options
      real,              intent(out)   :: dxn
        !! output new stepsize
      real,              intent(out)   :: Un(:)
        !! output state at end of step
      real,              intent(out)   :: dUndQ(:,:)
        !! output derivatives of Un wrt Q = [U0; P]
      real,              intent(out)   :: fn(:)
        !! output dUn/dx
      real,              intent(out)   :: dfndUn(:,:)
        !! output derivatives of fn wrt Un
      real,              intent(out)   :: dfndP(:,:)
        !! output derivatives of fn wrt P
      real,              intent(out)   :: polyn(:,:)
        !! output interpolation polynomial for this step
      real,              intent(out)   :: dpolyndQ(:,:,:)
        !! output derivatives of polyn wrt Q = [U0; P]
      real,              intent(inout) :: err
        !! output scalar error estimate
      logical,           intent(inout) :: status
        !! output success (true) or fail (false)
    end subroutine
  end interface

contains

  subroutine ode_solve(kernel, nS, fun, x0, x1, xsmp, U0, P, opt, status, res)
    procedure(ode_kernel)          :: kernel
      !! pointer to ode solver kernel
    integer,           intent(in)  :: nS
      !! number of stages of kernel
    procedure(ode_fun)             :: fun
      !! pointer to function to integrate
    real,              intent(in)  :: x0
      !! initial x
    real,              intent(in)  :: x1
      !! final x
    real,              intent(in)  :: xsmp(:)
      !! x sample points
    real,              intent(in)  :: U0(:)
      !! initial state
    real,              intent(in)  :: P(:)
      !! parameters
    type(ode_options), intent(in)  :: opt
      !! solver options
    logical,           intent(out) :: status
      !! true if ode was successfully solved, false if not
    type(ode_result),  intent(out) :: res
      !! output result object

    integer :: nU, nP
    logical :: final_step

    ! get system size and number of parameters
    nU = size(U0)
    nP = size(P)

    block
      integer :: i, j, rejection_counter
      integer :: ismp, nsmp, dsmp, ismpmax
      real    :: x, xold, dxk, dxn, err
      real    :: Uold(nU)
      real    :: Uk(nU), dUkdQ(nU,nU+nP), fk(nU), dfkdUk(nU,nU), dfkdP(nU,nP), polyk(nU,nS)
      real    :: Un(nU), dUndQ(nU,nU+nP), fn(nU), dfndUn(nU,nU), dfndP(nU,nP), polyn(nU,nS), dpolyndQ(nU,nU+nP,nS)

      ! initial x and dx
      x    = x0
      xold = x0
      dxk = (x1 - x0) * opt%initial_rx
      final_step = (sign(1.0, dxk) * (x + dxk - x1) > - x1 * opt%min_rx)
      if (final_step) dxk = x1 - x

      ! reset step and rejection counter
      res%nsteps = 0
      rejection_counter = 0

      ! initial state
      Uk    = U0
      Uold  = U0
      dUkdQ = 0
      do i = 1, nU
        dUkdQ(i,i) = 1.0
      end do

      ! number of samples
      nsmp = size(xsmp)

      m4_divert(m4_ifdef({m4_debug},0,-1))
        ! check if requested sample points are sorted
        do ismp = 2, nsmp
          if (xsmp(ismp) < xsmp(ismp-1)) call program_error("sample points xsmp must be sorted in ascending order")
        end do

        ! check if requested sample points are inside of interval
        if ((xsmp(1) < min(x0, x1)) .or. (xsmp(nsmp) > max(x0, x1))) then
          call program_error("sample points xsmp must lie inside of interval [min(x0, x1), max(x0, x1)]")
        end if
      m4_divert(0)

      ! init samples
      allocate(res%Usmp(    nU,   nsmp), source = 0.0)
      allocate(res%dUsmpdU0(nU,nU,nsmp), source = 0.0)
      allocate(res%dUsmpdP (nU,nP,nsmp), source = 0.0)

      ! go forward or backward through sample points
      if (x1 > x0) then
        ismp    = 1
        dsmp    = 1
        ismpmax = nsmp
      else
        ismp    = nsmp
        dsmp    = -1
        ismpmax = 1
      end if

      ! eval f and derivatives at U0
      call fun(x0, U0, P, status, f = fk, dfdU = dfkdUk)
      if (.not. status) return

      ! init interpolation polynomial (linear extrapolation)
      polyk      = 0
      polyk(:,1) = fk

      ! initial values for error control
      status = .false.
      err    = 1e99

      ! main loop
      do while (.true.)
        ! kernel
        call kernel(fun, xold, Uold, x, dxk, Uk, dUkdQ, fk, dfkdUk, dfkdP, polyk, P, opt, &
                                        dxn, Un, dUndQ, fn, dfndUn, dfndP, polyn, dpolyndQ, err, status)

        if (.not. status) then
          ! reject step
          rejection_counter = rejection_counter + 1
          if (rejection_counter > opt%max_rejected) return
        else
          ! accept step, reset rejection counter
          rejection_counter = 0

          ! advance x
          xold = x
          x = x + dxk
          if (final_step) x = x1

          ! get samples
          if (ismp /= ismpmax+dsmp) then
            do while ((min(xold, x) <= xsmp(ismp)) .and. (max(xold, x) >= xsmp(ismp)))
              res%Usmp(    :  ,ismp) = Uk
              res%dUsmpdU0(:,:,ismp) = dUkdq(:,1:nU)
              res%dUsmpdP( :,:,ismp) = dUkdq(:,(nU+1):(nU+nP))
              do j = 1, nS
                res%Usmp(    :  ,ismp) = res%Usmp(    :,  ismp) + polyn(   :               ,j) * (xsmp(ismp) - xold)**j
                res%dUsmpdU0(:,:,ismp) = res%dUsmpdU0(:,:,ismp) + dpolyndQ(:,          1:nU,j) * (xsmp(ismp) - xold)**j
                res%dUsmpdP( :,:,ismp) = res%dUsmpdP( :,:,ismp) + dpolyndQ(:,(nU+1):(nU+nP),j) * (xsmp(ismp) - xold)**j
              end do

              ! go to next sample point
              ismp = ismp + dsmp
              if (ismp == ismpmax + dsmp) exit
            end do
          end if

          ! copy n state to k state
          Uold   = Uk
          Uk     = Un
          dUkdQ  = dUndQ
          fk     = fn
          dfkdUk = dfndUn
          polyk  = polyn

          ! update step counter
          res%nsteps = res%nsteps + 1

          ! finished?
          if (final_step) exit
        end if

        ! adjust step size
        dxk = dxn
        final_step = (sign(1.0, dxk) * (x + dxk - x1) > - x1 * opt%min_rx)
        if (final_step) dxk = x1 - x
      end do
    end block
  end subroutine

  subroutine ode_options_init(this, nU, atol, rtol, newton_atol, newton_rtol, newton_max_it, max_rejected, min_rx, initial_rx)
    !! initialize ode solver options
    class(ode_options), intent(out) :: this
    integer,            intent(in)  :: nU
      !! system size
    real,    optional,  intent(in)  :: atol(:)
      !! absolute tolerance for each component (default: 1e-16)
    real,    optional,  intent(in)  :: rtol(:)
      !! relative tolerance for each component (default: 1e-12)
    real,    optional,  intent(in)  :: newton_atol
      !! absolute tolerance for newton iteration (default: 1e-16)
    real,    optional,  intent(in)  :: newton_rtol
      !! relative tolerance for newton iteration (default: 1e-14)
    integer, optional,  intent(in)  :: newton_max_it
      !! maximum number of allowed newton iterations (default: 5)
    integer, optional,  intent(in)  :: max_rejected
      !! maximum number of allowed rejected steps (default: 10)
    real,    optional,  intent(in)  :: min_rx
      !! minimum relative stepsize (default: 1e-8)
    real,    optional,  intent(in)  :: initial_rx
      !! initial stepsize, normalized to total interval size (default: 0.125)

    allocate (this%atol(nU), source = 1e-16)
    if (present(atol)) this%atol = atol
    allocate (this%rtol(nU), source = 1e-12)
    if (present(rtol)) this%rtol = rtol
    this%newton_atol = 1e-16
    if (present(newton_atol)) this%newton_atol = newton_atol
    this%newton_rtol = 1e-14
    if (present(newton_rtol)) this%newton_rtol = newton_rtol
    this%newton_max_it = 5
    if (present(newton_max_it)) this%newton_max_it = newton_max_it
    this%max_rejected = 10
    if (present(max_rejected)) this%max_rejected = max_rejected
    this%min_rx = 1e-8
    if (present(min_rx)) this%min_rx = min_rx
    this%initial_rx = 0.125
    if (present(initial_rx)) this%initial_rx = initial_rx
  end subroutine

end module
