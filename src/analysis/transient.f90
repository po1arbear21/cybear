m4_include(../util/macro.f90.inc)

module transient_m

  use bin_search_m,    only: bin_search, BS_LESS
  use error_m,         only: assert_failed, program_error
  use esystem_m,       only: esystem
  use gmres_m,         only: gmres_options
  use input_src_m,     only: input_src, periodic_src
  use math_m,          only: interp1, PI
  use matrix_m,        only: matrix_real, matrix_add, sparse_real
  use newton_m,        only: newton, newton_opt
  use normalization_m, only: denorm

  implicit none

  !FIXME: implement for dense system

  private
  public backward_euler, bdf2, mbdf2, trapezoidal, tr_bdf2, transient

  type, abstract :: transient
    !! Base transient type

    type(esystem), pointer    :: sys => null()
      !! pointer to corresponding equation system
    real,         allocatable :: x(:,:)
      !! data for each time step
    real                      :: t_0
      !! starting time
    real                      :: t_1
      !! ending time
    real                      :: delta_t
      !! initial / overall time-step
    integer                   :: n
      !! number of time steps
    real,         allocatable :: t(:)
      !! corresponding time array
    character(:), allocatable :: name
      !! transient methods name
    type(newton_opt)          :: nopt
      !! options for the newton solver
    type(gmres_options)       :: gopt
      !! options for the iterative solver in each newton iteration (only used when nopt%it_solver == true)
    type(sparse_real)         :: df, dft, dfp, mat
      !! matrices
    logical                   :: log
      !! enable/disable logging (default: false)
  contains
    procedure :: init       => transient_init
    procedure :: run        => transient_run
    procedure :: next_time  => transient_next_time
    procedure :: time_alloc => transient_time_alloc
    generic   :: select     => transient_select_ind, transient_select_time

    procedure, private :: transient_select_ind, transient_select_time

    procedure(transient_solve), deferred :: solve
  end type

  abstract interface
    subroutine transient_solve(this, i, input)
      import transient, input_src
      class(transient), target,   intent(inout) :: this
      integer,                    intent(in)    :: i
        !! current iteration step
      class(input_src), optional, intent(in)    :: input
        !! optional input parameter
    end subroutine

  end interface

  type, extends(transient) :: backward_euler
    !! Backward-Euler method
  contains
    procedure :: solve => backward_euler_solve
  end type

  type, extends(transient) :: trapezoidal
    !! Trapezoidal method

    real, allocatable :: f_old(:), f_new(:)
      !! residual work arrays
  contains
    procedure :: solve => trapezoidal_solve
  end type

  type, abstract, extends(transient) :: bdf2_base
    !! BDF2 and MBDF2 method base

    real :: a0, a1, a2
      !! coefficients
  contains
    procedure :: solve => bdf2_base_solve
  end type

  type, extends(bdf2_base) :: bdf2
    !! BDF2 method
  end type

  type, extends(bdf2_base) :: mbdf2
    !! MBDF2 method

    real :: z
      !! frequency dependent parameter
  end type

  type, extends(transient) :: tr_bdf2
    !! TRBDF2 method (combination of Trapezoidal and BDF2)

    logical           :: adaptive_timestep
      !! adaptive (true) or constant (false) time-step

    real, allocatable :: f_old(:), f_new(:), f_gam(:), x_gam(:)
      !! residual and solution work arrays
  contains
    procedure :: solve     => tr_bdf2_solve
    procedure :: next_time => tr_bdf2_next_time
  end type

contains

  subroutine transient_init(this, sys, t_0, t_1, delta_t, n, nopt, gopt, adaptive, log)
    !! initialize the transient type
    class(transient),              intent(out)   :: this
    type(esystem),       target,   intent(inout) :: sys
      !! equation system
    real,                          intent(in)    :: t_0
      !! starting time
    real,                          intent(in)    :: t_1
      !! ending time
    real,                optional, intent(in)    :: delta_t
      !! either provide the time-step of the iteration
    integer,             optional, intent(in)    :: n
      !! or the number of iterations
    type(newton_opt),    optional, intent(in)    :: nopt
      !! options for the newton solver
    type(gmres_options), optional, intent(in)    :: gopt
      !! options for the iterative solver in each newton iteration (only used when nopt%it_solver == true)
    logical,             optional, intent(in)    :: adaptive
      !! const (default, false) or adaptive (true) time-step (only used for tr-bdf2)
    logical,             optional, intent(in)    :: log
    !! enable/disable logging (default: false)

    ! set starting time
    this%t_0 = t_0

    ! set ending time
    this%t_1 = t_1

    ! initial timestep configuration
    if (present(n) .and. .not. present(delta_t)) then
      this%delta_t = (t_1 - t_0) / n
      this%n       = n
    else if (.not. present(n) .and. present(delta_t)) then
      this%delta_t = delta_t
      this%n       = int((t_1 - t_0) / delta_t)
    else
      m4_assert(.false.)
    end if

    ! set logging
    this%log = .false.
    if (present(log)) this%log = log

    ! save pointer to equation system for later use
    this%sys => sys

    ! parse optional newton options
    if (present(nopt)) then
      m4_assert(size(nopt%atol) == sys%n)
      this%nopt = nopt
    else
      call this%nopt%init(sys%n)
    end if
    if (present(gopt)) this%gopt = gopt

    ! set method dependent parameters
    select type (this)
    type is (backward_euler)
      this%name = "Backward_Euler"
      allocate (this%t(0:this%n))
      allocate (this%x(sys%n, 0:this%n))
    type is (bdf2)
      this%name = "BDF2"
      this%a0 =  1.5
      this%a1 = -2.0
      this%a2 =  0.5
      allocate (this%t(-1:this%n))
      allocate (this%x(sys%n, -1:this%n))
    type is (mbdf2)
      this%name = "MBDF2"
      this%z  = 2 * PI / this%n
      this%a0 = (this%z * cos(2*this%z) - this%z * cos(this%z)) / (sin(2*this%z) - 2 * sin(this%z))
      this%a1 = (this%z * sin(this%z)) / (cos(this%z) - 1)
      this%a2 = this%z / (2 * sin(this%z))
      allocate (this%t(-1:this%n))
      allocate (this%x(sys%n, -1:this%n))
    type is (trapezoidal)
      this%name = "Trapezoidal"
      allocate (this%t(0:this%n))
      allocate (this%x(sys%n, 0:this%n))
      allocate (this%f_old(this%sys%n), this%f_new(this%sys%n))
    type is (tr_bdf2)
      this%name = "TR-BDF2"
      this%adaptive_timestep = .false.
      if (present(adaptive)) this%adaptive_timestep = adaptive
      allocate (this%t(0:this%n))
      allocate (this%x(sys%n, 0:this%n))
      allocate (this%f_old(this%sys%n), this%f_new(this%sys%n), this%f_gam(this%sys%n), this%x_gam(this%sys%n))
    end select
  end subroutine

  subroutine transient_run(this, input)
    !! run the transient iteration
    class(transient),           intent(inout) :: this
    class(input_src), optional, intent(in)    :: input
      !! optional input parameter

    integer           :: i
    real, allocatable :: tmp_t(:), tmp_x(:,:)

    ! check input parameters
    if (present(input)) then
      m4_assert(input%n == this%sys%ninput)
      call this%sys%set_input(input%get(this%t_0))
    else
      m4_assert(this%sys%ninput == 0)
    end if

    ! initialize first iteration step
    this%t(  0) = this%t_0
    this%x(:,0) = this%sys%get_x()
    select type (this)
    type is (bdf2)
      this%t(  -1) = this%t_0 - this%delta_t
      this%x(:,-1) = this%sys%get_x()
    type is (mbdf2)
      this%t(  -1) = this%t_0 - this%delta_t
      this%x(:,-1) = this%sys%get_x()
    type is (trapezoidal)
      call this%sys%eval(f = this%f_old)
    type is (tr_bdf2)
      call this%sys%eval(f = this%f_old)
    end select

    ! jacobian wrt time derivative (constant)
    call this%sys%get_dft(this%dft)

    ! init matrix
    call this%mat%init(this%sys%n)

    ! time loop
    i = 0
    do while (this%t(i) < this%t_1)
      ! next timestep
      call this%next_time(i)
      if (this%log) print "(I0, ES25.16, x, A)", i, denorm(this%t(i), "fs"), "fs"

      ! solve system at current time
      call this%solve(i, input = input)
    end do

    ! truncate the arrays
    allocate(tmp_t(0:i))
    tmp_t(0:i) = this%t(0:i)
    call move_alloc(tmp_t, this%t)
    allocate (tmp_x(this%sys%n, 0:i))
    tmp_x(:,0:i) = this%x(:,0:i)
    call move_alloc(tmp_x, this%x)

    call this%dft%destruct()
  end subroutine

  subroutine transient_next_time(this, i)
    class(transient), intent(inout) :: this
    integer,          intent(inout) :: i
      !! current iterator of the loop

    ! advance time for fixed delta_t
    i = i + 1
    call this%time_alloc(i)
    this%t(i) = this%t(i-1) + this%delta_t
    if ((abs(this%t(i) - this%t_1) < 1e-10 * abs(this%t_1 - this%t_0)) .or. (this%t(i) > this%t_1)) this%t(i) = this%t_1
  end subroutine

  subroutine transient_time_alloc(this, n)
    class(transient), intent(inout) :: this
    integer,          intent(in)    :: n
      !! minimum new upper bound

    integer           :: l, u
    real, allocatable :: tmp_t(:), tmp_x(:,:)

    l = lbound(this%t, dim = 1)
    u = ubound(this%t, dim = 1)
    if (n > u) then
      allocate (tmp_t(l:2*n))
      tmp_t(l:u) = this%t(l:u)
      call move_alloc(tmp_t, this%t)
      allocate (tmp_x(this%sys%n, l:2*n))
      tmp_x(:,l:u) = this%x(:,l:u)
      call move_alloc(tmp_x, this%x)
    end if
  end subroutine

  subroutine tr_bdf2_next_time(this, i)
    class(tr_bdf2), intent(inout) :: this
    integer,        intent(inout) :: i
      !! current iterator of the loop

    if (i == 0) then
      call transient_next_time(this, i)
    elseif (this%adaptive_timestep) then
      call time_step_control()
    else
      call transient_next_time(this, i)
      this%f_old = this%f_new
    end if

  contains

    subroutine time_step_control()
      real, parameter   :: p = 2.0, erel = 1e-6, eabs = 1e-16, gamma = 2 - sqrt(2.0)
      real, parameter   :: k_g = (-3*gamma**2 + 4*gamma - 2) / (12*(2 - gamma))
      real, allocatable :: err(:)
      real              :: r

      ! estimate error
      allocate (err(this%sys%n))
      err = 2 * k_g * this%delta_t * ((1/gamma)*this%f_old +(1/(gamma*(gamma-1)))*this%f_gam + (1/(1-gamma))*this%f_new)
      r = norm2(err / (erel * this%x(:,i) + eabs)) / sqrt(real(this%sys%n))

      if (r <= 2) then
        ! accept timestep
        i = i + 1
        call this%time_alloc(i)
        this%f_old = this%f_new
      else
        ! reject timestep
        this%t(i) = this%t(i-1)
      end if

      ! adjust timestep
      this%delta_t = this%delta_t / (r ** (1 / (p + 1)))
      this%t(i) = this%t(i-1) + this%delta_t
      if ((abs(this%t(i) - this%t_1) < 1e-10 * abs(this%t_1 - this%t_0)) .or. (this%t(i) > this%t_1)) then
        this%t(i) = this%t_1
        this%delta_t = this%t(i) - this%t(i-1)
      end if
    end subroutine

  end subroutine

  subroutine transient_select_ind(this, i)
    !! save trasient results at index i (overwrite esystem variables)
    class(transient), intent(in) :: this
    integer,          intent(in) :: i

    ! checking if output is in the simulation range
    m4_assert(i >= lbound(this%x, 2))
    m4_assert(i <= ubound(this%x, 2))

    call this%sys%set_x(this%x(:,i))
  end subroutine

  subroutine transient_select_time(this, t)
    !! save transient interpolated results at time t (overwrite esystem variables)
    class(transient), intent(in) :: this
    real,             intent(in) :: t

    integer           :: ind
    real, allocatable :: x(:)

    ! check if output is in the simulation range
    m4_assert(t >= this%t_0)
    m4_assert(t <= this%t_1)

    ! allocate interpolated input array
    allocate(x(this%sys%n))

    ! interpolate x between points
    if (t == this%t_0) then
      x = this%x(:,lbound(this%x, dim = 2))
    elseif (t == this%t_1) then
      x = this%x(:,ubound(this%x, dim = 2))
    else
      ind = bin_search(this%t, t, BS_LESS) - 1 + lbound(this%t, dim = 1)
      x = this%x(:,ind) + (t - this%t(ind)) / (this%t(ind+1) - this%t(ind)) * (this%x(:,ind+1) - this%x(:,ind))
    end if

    call this%sys%set_x(x)
  end subroutine

  subroutine backward_euler_solve(this, i, input)
    !! solve Backward-Euler timestep
    class(backward_euler), target, intent(inout) :: this
    integer,                       intent(in)    :: i
      !! current iteration step
    class(input_src), optional,    intent(in)    :: input
      !! optional input parameter

    real :: p(0)

    ! set input
    if (present(input)) then
      call this%sys%set_input(input%get(this%t(i)))
    end if

    ! newton iteration
    call newton(fun, p, this%nopt, this%x(:,i-1), this%x(:,i), gmres_opt = this%gopt)
    call this%mat%reset()

  contains

    subroutine fun(x, p, f, dfdx, dfdx_prec, dfdp)
      real,                                  intent(in)  :: x(:)
        !! arguments
      real,                                  intent(in)  :: p(:)
        !! parameters
      real,               optional,          intent(out) :: f(:)
        !! output function values
      class(matrix_real), optional, pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
      class(matrix_real), optional, pointer, intent(out) :: dfdx_prec
        !! optional output pointer to preconditioner jacobian of f wrt x
      real,               optional,          intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p

      ! params not needed
      m4_ignore(p)
      m4_assert(present(dfdx))
      m4_assert(.not. present(dfdp))
      m4_ignore(dfdp)

      ! set x in sys
      call this%sys%set_x(x)

      ! compute residue, jacobian, and possible preconditioner
      if (this%nopt%it_solver) then
        m4_assert(present(dfdx_prec))
        call this%sys%eval(f = f, df = this%df, dfp = this%dfp)
        dfdx_prec => this%dfp
      else
        m4_assert(.not. present(dfdx_prec))
        call this%sys%eval(f = f, df = this%df)
      end if

      ! construct f
      call this%dft%mul_vec((x - this%x(:,i-1)) / this%delta_t, f, fact_y = 1.0)

      ! construct Jacobian of f
      call matrix_add(this%df, this%dft, this%mat, fact2 = (1 / this%delta_t))
      dfdx => this%mat
    end subroutine

  end subroutine

  subroutine trapezoidal_solve(this, i, input)
    !! solve Trapezoidal timestep
    class(trapezoidal), target, intent(inout) :: this
    integer,                    intent(in)    :: i
      !! current iteration step
    class(input_src), optional, intent(in)    :: input
      !! optional input parameter

    real :: p(0)

    ! set input
    if (present(input)) then
      call this%sys%set_input(input%get(this%t(i)))
    end if

    ! newton iteration
    call newton(fun, p, this%nopt, this%x(:,i-1), this%x(:,i), gmres_opt = this%gopt)
    call this%mat%reset()

    ! save f
    this%f_old = this%f_new

  contains

    subroutine fun(x, p, f, dfdx, dfdx_prec, dfdp)
      real,                                  intent(in)  :: x(:)
        !! arguments
      real,                                  intent(in)  :: p(:)
        !! parameters
      real,               optional,          intent(out) :: f(:)
        !! output function values
      class(matrix_real), optional, pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
      class(matrix_real), optional, pointer, intent(out) :: dfdx_prec
        !! optional output pointer to preconditioner jacobian of f wrt x
      real,               optional,          intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p

      ! params not needed
      m4_assert(size(p) == 0)
      m4_ignore(p)
      m4_assert(present(dfdx))
      m4_assert(.not. present(dfdp))
      m4_ignore(dfdp)

      ! set x in sys
      call this%sys%set_x(x)

      ! compute residue, jacobian, and possible preconditioner
      if (this%nopt%it_solver) then
        m4_assert(present(dfdx_prec))
        call this%sys%eval(f = this%f_new, df = this%df, dfp = this%dfp)
        dfdx_prec => this%dfp
      else
        m4_assert(.not. present(dfdx_prec))
        call this%sys%eval(f = this%f_new, df = this%df)
      end if

      ! construct f
      f = this%f_new + this%f_old
      call this%dft%mul_vec((x - this%x(:,i-1)) / this%delta_t, f, fact_y = 0.5)

      ! construct Jacobian of f
      call matrix_add(this%df, this%dft, this%mat, fact1 = 0.5, fact2 = (1 / this%delta_t))

      dfdx => this%mat
    end subroutine

  end subroutine

  subroutine bdf2_base_solve(this, i, input)
    !! solve BDF2 or MBDF2 timestep
    class(bdf2_base), target,   intent(inout) :: this
    integer,                    intent(in)    :: i
      !! current iteration step
    class(input_src), optional, intent(in)    :: input
      !! optional input parameter

    real :: p(0)

    ! set input
    if (present(input)) then
      call this%sys%set_input(input%get(this%t(i)))
    end if

    ! newton iteration
    call newton(fun, p, this%nopt, this%x(:,i-1), this%x(:,i), gmres_opt = this%gopt)
    call this%mat%reset()

  contains

    subroutine fun(x, p, f, dfdx, dfdx_prec, dfdp)
      real,                                  intent(in)  :: x(:)
        !! arguments
      real,                                  intent(in)  :: p(:)
        !! parameters
      real,               optional,          intent(out) :: f(:)
        !! output function values
      class(matrix_real), optional, pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
      class(matrix_real), optional, pointer, intent(out) :: dfdx_prec
        !! optional output pointer to preconditioner jacobian of f wrt x
      real,               optional,          intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p

      ! params not needed
      m4_assert(size(p) == 0)
      m4_ignore(p)
      m4_assert(present(dfdx))
      m4_assert(.not. present(dfdp))
      m4_ignore(dfdp)

      ! save input variable
      call this%sys%set_x(x)

      ! compute residue, jacobian, and possible preconditioner
      if (this%nopt%it_solver) then
        m4_assert(present(dfdx_prec))
        call this%sys%eval(f = f, df = this%df, dfp = this%dfp)
        dfdx_prec => this%dfp
      else
        m4_assert(.not. present(dfdx_prec))
        call this%sys%eval(f = f, df = this%df)
      end if

      ! construct f
      call this%dft%mul_vec((this%a0 * x + this%a1 * this%x(:,i-1) + this%a2 * this%x(:,i-2)) / this%delta_t, f, fact_y = 1.0)

      ! construct Jacobian of f
      call matrix_add(this%df, this%dft, this%mat, fact1 = 1.0, fact2 = (this%a0 / this%delta_t))

      dfdx => this%mat
    end subroutine

  end subroutine

  subroutine tr_bdf2_solve(this, i, input)
    !! solve TRBDF2 timestep
    class(tr_bdf2), target,     intent(inout) :: this
    integer,                    intent(in)    :: i
      !! current iteration step
    class(input_src), optional, intent(in)    :: input
      !! optional input parameter

    real, parameter :: gamma = 2.0 - sqrt(2.0)
    real            :: p(0)

    ! set input at gamma step
    if (present(input)) then
      call this%sys%set_input(input%get(this%t(i-1) + this%delta_t * gamma))
    end if

    ! solve trapezoidal with newton iteration
    call newton(tr_fun, p, this%nopt, this%x(:,i-1), this%x_gam, gmres_opt = this%gopt)
    call this%mat%reset()

    ! set input at new time
    if (present(input)) then
      call this%sys%set_input(input%get(this%t(i)))
    end if

    ! solve bdf2 with newton iteration
    call newton(bdf2_fun, p, this%nopt, this%x_gam, this%x(:,i), gmres_opt = this%gopt)
    call this%mat%reset()

  contains

    subroutine tr_fun(x, p, f, dfdx, dfdx_prec, dfdp)
      real,                                  intent(in)  :: x(:)
        !! arguments
      real,                                  intent(in)  :: p(:)
        !! parameters
      real,               optional,          intent(out) :: f(:)
        !! output function values
      class(matrix_real), optional, pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
      class(matrix_real), optional, pointer, intent(out) :: dfdx_prec
        !! optional output pointer to preconditioner jacobian of f wrt x
      real,               optional,          intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p

      ! params not needed
      m4_assert(size(p) == 0)
      m4_ignore(p)
      m4_assert(present(dfdx))
      m4_assert(.not. present(dfdp))
      m4_ignore(dfdp)

      ! set x in sys
      call this%sys%set_x(x)

      ! compute residue, jacobian, and possible preconditioner
      if (this%nopt%it_solver) then
        m4_assert(present(dfdx_prec))
        call this%sys%eval(f = this%f_gam, df = this%df, dfp = this%dfp)
        dfdx_prec => this%dfp
      else
        m4_assert(.not. present(dfdx_prec))
        call this%sys%eval(f = this%f_gam, df = this%df)
      end if

      ! construct f
      f = this%f_gam + this%f_old
      call this%dft%mul_vec((x - this%x(:,i-1)) / (this%delta_t * gamma), f, fact_y = 0.5)

      ! construct Jacobian of f
      call matrix_add(this%df, this%dft, this%mat, fact1 = 0.5, fact2 = 1 / (this%delta_t * gamma))

      dfdx => this%mat
    end subroutine

    subroutine bdf2_fun(x, p, f, dfdx, dfdx_prec, dfdp)
      real,                                  intent(in)  :: x(:)
        !! arguments
      real,                                  intent(in)  :: p(:)
        !! parameters
      real,               optional,          intent(out) :: f(:)
        !! output function values
      class(matrix_real), optional, pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
      class(matrix_real), optional, pointer, intent(out) :: dfdx_prec
        !! optional output pointer to preconditioner jacobian of f wrt x
      real,               optional,          intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p

      ! params not needed
      m4_assert(size(p) == 0)
      m4_ignore(p)
      m4_assert(present(dfdx))
      m4_assert(.not. present(dfdp))
      m4_ignore(dfdp)

      ! save input variable
      call this%sys%set_x(x)

      ! compute residue, jacobian, and possible preconditioner
      if (this%nopt%it_solver) then
        m4_assert(present(dfdx_prec))
        call this%sys%eval(f = this%f_new, df = this%df, dfp = this%dfp)
        dfdx_prec => this%dfp
      else
        m4_assert(.not. present(dfdx_prec))
        call this%sys%eval(f = this%f_new, df = this%df)
      end if

      ! construct f
      f = this%f_new
      call this%dft%mul_vec(((2-gamma)/(1-gamma))*x + (1/(gamma*(gamma-1)))*this%x_gam + ((1-gamma)/gamma)*this%x(:,i-1), f, fact_y = this%delta_t)

      ! construct Jacobian of f
      call matrix_add(this%df, this%dft, this%mat, fact1 = this%delta_t, fact2 = (2-gamma)/(1-gamma))

      dfdx => this%mat
    end subroutine

  end subroutine

end module
