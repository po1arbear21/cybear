module ode_m
  use error_m
  use high_precision_m
  implicit none

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
  contains
    procedure :: init => ode_options_init ! initialize
  end type

  type ode_result
    !! results returned by ode solver

    real, allocatable :: Us(:,:)
      !! return values at sample points (nU x size(xs))

    real, allocatable :: U1(:)
      !! return final state
    real, allocatable :: dU1dU0(:,:)
      !! return derivatives of U1 wrt U0
    real, allocatable :: dU1dP(:,:)
      !! return derivatives of U1 wrt P

    real, allocatable :: UA(:)
      !! return average state in interval
    real, allocatable :: dUAdU0(:,:)
      !! return derivatives of UA wrt U0
    real, allocatable :: dUAdP(:,:)
      !! return derivatives of UA wrt P

    integer :: nsteps
      !! return number of steps taken
  end type

  interface
    subroutine ode_fun(x, U, P, f, dfdU, dfdP)
      real,              intent(in)  :: x
        !! x coordinate
      real,              intent(in)  :: U(:)
        !! state (nU)
      real,              intent(in)  :: P(:)
        !! parameters (nP)
      real, optional,    intent(out) :: f(:)
        !! output dU/dx (nU)
      real, optional,    intent(out) :: dfdU(:,:)
        !! output derivatives of f wrt U (nU,nU)
      real, optional,    intent(out) :: dfdP(:,:)
        !! output derivatives of f wrt P (nU,nP)
    end subroutine

    subroutine ode_kernel(fun, xold, x, dxk, Uk, dUkdQ, fk, dfkdUk, polyk, &
      P, opt, dxn, Un, dUndQ, fn, dfndUn, polyn, dpolyndQ, err, status)
      import ode_options, ode_fun
      procedure(ode_fun), pointer, intent(in)    :: fun
        !! function to integrate
      real,                        intent(in)    :: xold
        !! initial x position of previous step
      real,                        intent(in)    :: x
        !! initial x position of step
      real,                        intent(in)    :: dxk
        !! stepsize
      real,                        intent(in)    :: Uk(:)
        !! state at beginning of step
      real,                        intent(in)    :: dUkdQ(:,:)
        !! derivatives of Uk wrt Q = [U0; P]
      real,                        intent(in)    :: fk(:)
        !! dUk/dx
      real,                        intent(in)    :: dfkdUk(:,:)
        !! derivatives of fk wrt Uk
      real,                        intent(in)    :: polyk(:,:)
        !! interpolation polynomial (from previous step)
      real,                        intent(in)    :: P(:)
        !! parameters
      type(ode_options),           intent(in)    :: opt
        !! options
      real,                        intent(out)   :: dxn
        !! output new stepsize
      real,                        intent(out)   :: Un(:)
        !! output state at end of step
      real,                        intent(out)   :: dUndQ(:,:)
        !! output derivatives of Un wrt Q = [U0; P]
      real,                        intent(out)   :: fn(:)
        !! output dUn/dx
      real,                        intent(out)   :: dfndUn(:,:)
        !! output derivatives of fn wrt Un
      real,                        intent(out)   :: polyn(:,:)
        !! output interpolation polynomial for this step
      real,                        intent(out)   :: dpolyndQ(:,:,:)
        !! output derivatives of polyn wrt Q = [U0; P]
      real,                        intent(inout) :: err
        !! output scalar error estimate
      logical,                     intent(inout) :: status
        !! output success (true) or fail (false)
    end subroutine
  end interface

contains

  subroutine ode_solve(kernel, nS, fun, x0, x1, U0, P, opt, res, xs)
    procedure(ode_kernel), pointer, intent(in)  :: kernel
      !! pointer to ode solver kernel
    integer,                        intent(in)  :: nS
      !! number of stages of kernel
    procedure(ode_fun),    pointer, intent(in)  :: fun
      !! pointer to function to integrate
    real,                           intent(in)  :: x0
      !! initial x
    real,                           intent(in)  :: x1
      !! final x
    real,                           intent(in)  :: U0(:)
      !! initial state
    real,                           intent(in)  :: P(:)
      !! parameters
    type(ode_options),              intent(in)  :: opt
      !! solver options
    type(ode_result),               intent(out) :: res
      !! output result object
    real, optional,                 intent(in)  :: xs(:)
      !! x sample points

    ! local variables
    integer :: nU, nP

    ! get system size and number of parameters
    nU = size(U0)
    nP = size(P)

    block
      integer       :: i, j, is, rejection_counter
      real          :: x, xold, dxk, dxn, err
      real          :: Uk(nU), dUkdQ(nU,nU+nP), fk(nU), dfkdUk(nU,nU), polyk(nU,nS)
      real          :: Un(nU), dUndQ(nU,nU+nP), fn(nU), dfndUn(nU,nU), polyn(nU,nS), dpolyndQ(nU,nU+nP,nS)
      type(hp_real) :: hp_UA(nU), hp_dUAdQ(nU,nU+nP)
      logical       :: status

      ! initial x and dx
      x    = x0
      xold = x0
      dxk = (x1 - x0) / 8

      ! reset step and rejection counter
      res%nsteps = 0
      rejection_counter = 0

      ! initial state
      Uk    = U0
      dUkdQ = 0
      do i = 1, nU
        dUkdQ(i,i) = 1.0
      end do

      ! init samples
      if (present(xs)) then
        allocate(res%Us(nU,size(xs)), source = 0d0)
        is = 1
      end if

      ! eval f and derivatives at U0
      call fun(x0, U0, P, f = fk, dfdU = dfkdUk)

      ! init interpolation polynomial (linear extrapolation)
      polyk      = 0
      polyk(:,1) = fk

      ! init average state
      hp_UA    = real_to_hp(0.0)
      hp_dUAdQ = real_to_hp(0.0)

      ! initial values for error control
      status = .false.
      err    = 1e99

      ! main loop
      do while (abs(x - x0) < abs(x1 - x0))
        ! kernel
        call kernel(fun, xold, x, dxk, Uk, dUkdQ, fk, dfkdUk, polyk, P, opt, &
                                  dxn, Un, dUndQ, fn, dfndUn, polyn, dpolyndQ, err, status)

        if (.not. status) then
          ! reject step
          rejection_counter = rejection_counter + 1
          if (rejection_counter > opt%max_rejected) then
            print "(1A, 1E64.56)", "x0 = ", x0
            print "(1A, 1E64.56)", "x1 = ", x1
            do i = 1, nU
              print "(1A, 1I1, 1A, 1E64.56)", "U0(", i, ") = ", U0(i)
            end do
            do i = 1, nP
              print "(1A, 1I1, 1A, 1E64.56)", "P(", i, ") = ", P(i)
            end do

            call program_error("Can not solve ode, limit for rejected steps reached!")
          end if
        else
          ! accept step, reset rejection counter
          rejection_counter = 0

          ! update average state by integrating the interpolation polynomial
          hp_UA    = hp_UA    + dxk * Uk
          hp_dUAdQ = hp_dUAdQ + dxk * dUkdQ
          do j = 1, nS
            hp_UA    = hp_UA    + dxk**(j+1) / real(j+1) * polyn(:,j)
            hp_dUAdQ = hp_dUAdQ + dxk**(j+1) / real(j+1) * dpolyndQ(:,:,j)
          end do

          ! sample
          if (present(xs)) then
            if (is <= size(xs)) then
              do while (xs(is) <= x + dxk)
                res%Us(:,is) = Uk + polyn(:,1) * (xs(is) - x) + polyn(:,2) * (xs(is) - x)**2 + polyn(:,3) * (xs(is) - x)**3
                is = is + 1
                if (is > size(xs)) exit
              end do
            end if
          end if

          ! advance x
          xold = x
          x    = x + dxk

          ! copy n state to k state
          Uk     = Un
          dUkdQ  = dUndQ
          fk     = fn
          dfkdUk = dfndUn
          polyk  = polyn

          ! update step counter
          res%nsteps = res%nsteps + 1
        end if

        ! adjust step size
        dxk = dxn
        if (abs(x + dxk - x0) > abs(x1 - x0)) dxk = x1 - x
      end do

      ! set end state and extract derivatives from dUkdQ
      res%U1 = Uk
      res%dU1dU0 = dUkdQ(:,1:nU)
      res%dU1dP  = dUkdQ(:,(nU+1):(nU+nP))

      ! average state
      res%UA     = hp_to_real(hp_UA                     ) / (x1 - x0)
      res%dUAdU0 = hp_to_real(hp_dUAdQ(:,1:nU)          ) / (x1 - x0)
      res%dUAdP  = hp_to_real(hp_dUAdQ(:,(nU+1):(nU+nP))) / (x1 - x0)
    end block
  end subroutine

  subroutine ode_options_init(this, nU, atol, rtol, newton_atol, newton_rtol, newton_max_it, max_rejected, min_rx)
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

    allocate (this%atol(nU), source = 1e-16)
    allocate (this%rtol(nU), source = 1e-12)

    if (present(atol)) then
      this%atol = atol
    end if
    if (present(rtol)) then
      this%rtol = rtol
    end if
    if (present(newton_atol)) then
      this%newton_atol = newton_atol
    else
      this%newton_atol = 1e-16
    end if
    if (present(newton_rtol)) then
      this%newton_rtol = newton_rtol
    else
      this%newton_rtol = 1e-14
    end if
    if (present(newton_max_it)) then
      this%newton_max_it = newton_max_it
    else
      this%newton_max_it = 5
    end if
    if (present(max_rejected)) then
      this%max_rejected = max_rejected
    else
      this%max_rejected = 10
    end if
    if (present(min_rx)) then
      this%min_rx = min_rx
    else
      this%min_rx = 1e-8
    end if
  end subroutine

end module