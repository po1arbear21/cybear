m4_include(../util/macro.f90.inc)

module transient_m

  use bin_search_m,    only: bin_search
  use color_m,         only: COL_DEFAULT, COL_MAGENTA
  use container_m,     only: container, STORAGE_WRITE, STORAGE_READ, DYNAMIC_APP, ERR_ALREADY_EXISTS
  use error_m,         only: assert_failed, program_error
  use esystem_m,       only: esystem
  use gmres_m,         only: gmres_options, gmres
  use ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_positive_inf, ieee_negative_inf
  use input_src_m,     only: input_src
  use logging_m,       only: logging
  use matop_m,         only: single_matop_real
  use matrix_m,        only: SPSOLVER_PARDISO, SOLVER_GMRES, block_real, matrix_add
  use normalization_m, only: norm, denorm
  use string_m,        only: string, new_string
  use util_m,          only: int2str, real2str
  use variable_m,      only: variable_ptr, variable_real

  implicit none

  ! Todo: better handling of jump conditions within transient simulation

  private
  public transient
  public TRANS_BE, TRANS_TRAPZ, TRANS_BDF2, TRANS_MBDF2, TRANS_TRBDF2

  !! methods
  integer, parameter :: TRANS_BE     = 1
  integer, parameter :: TRANS_TRAPZ  = 2
  integer, parameter :: TRANS_BDF2   = 3
  integer, parameter :: TRANS_MBDF2  = 4
  integer, parameter :: TRANS_TRBDF2 = 5
  !! bdf2 and mbdf2 parameters
  ! Todo: add BE, BDF2 and MBDF2
  !! tr-bdf2 parameters
  real, parameter :: GAM = 2 - sqrt(2.0)
  ! parameters for error erstimation
  real, parameter :: POW = 1.0/3.0, C1 = (GAM-1)/3.0, C2 = 1.0/3.0, C3 = -GAM/3.0
  ! parameter for stepsize regulation (Matlab : 0.7, Bank (1985): 1.0)
  real, parameter :: THETA_1 = 1.0

  type transient
    !! transient simulation -> solve ODE from esystem: M*dx/dt + f(t,x) = 0

    ! general data
    type(esystem), pointer    :: sys => null()
      !! pointer to corresponding equation system
    real,         allocatable :: t(:)
      !! time points
    integer,      allocatable :: n_steps(:)
      !! number of time steps for each section of transient simulation, not including the start point
    real,         allocatable :: x(:,:)
      !! data for each time step (sys%n, size(this%t)), only used if results are stored in RAM memory
    logical                   :: use_ram
      !! save data also in RAM memory (true, default) or only on disk (false)
    logical                   :: log
      !! enable/disable logging (default false)
    character(:), allocatable :: msg
      !! message to print each Newton iteration if logging is enabled (default "Transient: ")

    ! ODE solver parameters
    integer           :: method
      !! ODE solver
    logical           :: adaptive
      !! use adaptive time steps (true) or not (false, default)
    real              :: reduction_factor
      !! adaptive: time step reduction factor if Newton does not converge (default 0.5)
    real, allocatable :: eabs(:)
      !! adaptive: absolute error tolerance for error estimate (default 1e-16)
    real, allocatable :: erel(:)
      !! adaptive: relative error tolerance for error estimate (default 1e-6)

    ! data for Newton
    integer             :: solver
      !! matrix solver (default SPSOLVER_PARDISO)
    type(gmres_options) :: gopt
      !! options for gmres in case an iterative solver is used (default atol 1e-16, default rtol 1e-12)
    real, allocatable   :: atol(:)
      !! absolute solution tolerance (default 1e-16)
    real, allocatable   :: rtol(:)
      !! relative solution tolerance (default 1e-12)
    real, allocatable   :: ftol(:)
      !! absolute residual tolerance (default 0.0)
    real, allocatable   :: dx_lim(:)
      !! limit for Newton update (default inf)
    real, allocatable   :: dx_lim_rel(:)
      !! limit for Newton update, relative to current value of variable (default inf)
    real, allocatable   :: xmin(:)
      !! lower solution bound (default -inf)
    real, allocatable   :: xmax(:)
      !! upper solution bound (default +inf)
    integer             :: min_it
      !! minimum number of Newton iterations (default 0)
    integer             :: max_it
      !! maximum number of Newton iterations (default huge)

    ! data for output (default no output)
    type(string), allocatable :: output_vars(:)
      !! list of variables for output
    integer                   :: vars_delta_it
      !! write variables at each vars_delta_it-th step...
    real,         allocatable :: vars_t(:)
      !! ... or write variables at these exact time points (relative tolerance 1e-10)
    character(:), allocatable :: varfile
      !! output file for variables
    integer                   :: cache_delta_it
      !! save this%system%get_x() at each cache_delta_it-th step...
    real,         allocatable :: cache_t(:)
      !! ... or save this%system%get_x() at these exact time points (relative tolerance 1e-10)
    character(:), allocatable :: cachefile
      !! output file for cache
  contains
    procedure :: init              => transient_init
    procedure :: set_newton_params => transient_set_newton_params
    procedure :: set_var_params    => transient_set_var_params
    procedure :: init_output       => transient_init_output
    procedure :: init_cache        => transient_init_cache
    procedure :: run               => transient_run
    generic   :: select_result     => transient_select_i, transient_select_t

    procedure, private :: transient_select_i, transient_select_t
    procedure, private :: solve_trapz  => transient_solve_trapz
    procedure, private :: solve_trbdf2 => transient_solve_trbdf2
    procedure, private :: newton_x     => transient_newton_x
    procedure, private :: newton_z     => transient_newton_z
    procedure, private :: time_alloc   => transient_time_alloc
    procedure, private :: output_grids => transient_output_grids
    procedure, private :: output       => transient_output
  end type

  abstract interface
    subroutine get_update(var, dvar, f)
      real, intent(in)  :: var(:)
      real, intent(out) :: dvar(:)
      real, intent(out) :: f(:)
    end subroutine

    function get_val(var1) result(var2)
      real, intent(in)  :: var1(:)
      real, allocatable :: var2(:)
    end function
  end interface

contains

  subroutine transient_init(this, sys, method, adaptive, reduction_factor, eabs, erel, use_ram, log, msg)
    !! initialize transient simulation
    class(transient),       intent(out) :: this
    type(esystem), target,  intent(in)  :: sys
      !! equation system
    integer,      optional, intent(in)  :: method
      !! ODE solver used for simulation (TRANS_BE, TRANS_TRAPZ, TRANS_BDF2, TRANS_MBDF2, TRANS_TRBDF2)
    logical,      optional, intent(in)  :: adaptive
      !! adaptive time step size (default false)
    real,         optional, intent(in)  :: reduction_factor
      !! adaptive: time step reduction factor if Newton does not converge (default 0.5)
    real,         optional, intent(in)  :: eabs
      !! adaptive: absolute error tolerance for error estimate (default 1e-16)
    real,         optional, intent(in)  :: erel
      !! adaptive: relative error tolerance for error estimate (default 1e-6)
    logical,      optional, intent(in)  :: use_ram
      !! save data also in RAM memory (true, default) or only on disk (false)
    logical,      optional, intent(in)  :: log
      !! enable/disable logging (default false)
    character(*), optional, intent(in)  :: msg
      !! message to print each Newton iteration if logging is enabled (default "Transient: ")

    this%sys => sys
    this%method = method
    this%adaptive = .false.
    if (present(adaptive)) this%adaptive = adaptive
    if (this%adaptive) then
      this%reduction_factor = 0.5
      if (present(reduction_factor)) this%reduction_factor = reduction_factor
      m4_assert(this%reduction_factor > 0.0 .and. this%reduction_factor < 1.0)
      allocate(this%eabs(sys%n), source = 1e-16)
      if (present(eabs)) this%eabs = eabs
      allocate(this%erel(sys%n), source = 1e-6)
      if (present(erel)) this%erel = erel
    else
      m4_assert(.not. present(reduction_factor))
      m4_assert(.not. present(eabs))
      m4_assert(.not. present(erel))
    end if
    this%use_ram = .true.
    if (present(use_ram)) this%use_ram = use_ram
    this%log = .false.
    if (present(log)) this%log = log
    this%msg = "Transient: "
    if (present(msg)) deallocate(this%msg)
    if (present(msg)) this%msg = msg

    this%solver = SPSOLVER_PARDISO
    allocate(this%atol(sys%n), source = 1e-16)
    allocate(this%rtol(sys%n), source = 1e-12)
    allocate(this%ftol(sys%n), source = 0.0)
    allocate(this%dx_lim(sys%n), source = ieee_value(1.0, ieee_positive_inf))
    allocate(this%dx_lim_rel(sys%n), source = ieee_value(1.0, ieee_positive_inf))
    allocate(this%xmin(sys%n), source = ieee_value(1.0, ieee_negative_inf))
    allocate(this%xmax(sys%n), source = ieee_value(1.0, ieee_positive_inf))
    this%min_it = 0
    this%max_it = huge(this%max_it)
  end subroutine

  subroutine transient_set_newton_params(this, solver, gopt, atol, rtol, ftol, dx_lim, dx_lim_rel, xmin, xmax, min_it, max_it)
    !! set parameters for Newton iteration, overwriting the whole array
    class(transient),              intent(inout) :: this
    integer,             optional, intent(in)    :: solver
      !! matrix solver
    type(gmres_options), optional, intent(in)    :: gopt
      !! options for gmres in case an iterative solver is used
    real,                optional, intent(in)    :: atol
      !! absolute solution tolerance
    real,                optional, intent(in)    :: rtol
      !! relative solution tolerance
    real,                optional, intent(in)    :: ftol
      !! absolute residual tolerance
    real,                optional, intent(in)    :: dx_lim
      !! limit for Newton update
    real,                optional, intent(in)    :: dx_lim_rel
      !! limit for Newton update, relative to current value of variable
    real,                optional, intent(in)    :: xmin
      !! lower solution bound
    real,                optional, intent(in)    :: xmax
      !! upper solution bound
    integer,             optional, intent(in)    :: min_it
      !! minimum number of Newton iterations
    integer,             optional, intent(in)    :: max_it
      !! maximum number of Newton iterations

    if (present(solver)) this%solver = solver
    if (present(gopt)) this%gopt = gopt
    if (present(atol)) this%atol = atol
    if (present(rtol)) this%rtol = rtol
    if (present(ftol)) this%ftol = ftol
    if (present(dx_lim)) this%dx_lim = dx_lim
    if (present(dx_lim_rel)) this%dx_lim_rel = dx_lim_rel
    if (present(xmin)) this%xmin = xmin
    if (present(xmax)) this%xmax = xmax
    if (present(min_it)) this%min_it = min_it
    if (present(max_it)) this%max_it = max_it
  end subroutine

  subroutine transient_set_var_params(this, vs_name, eabs, erel, atol, rtol, ftol, dx_lim, dx_lim_rel, xmin, xmax)
    !! set parameters for Newton iteration for specific vselector, leaving the rest of the array unchanged
    class(transient), intent(inout) :: this
    character(*),     intent(in)    :: vs_name
      !! name of vselector
    real, optional,   intent(in)    :: eabs
      !! adaptive: absolute error tolerance for error estimate
    real, optional,   intent(in)    :: erel
      !! adaptive: relative error tolerance for error estimate
    real, optional,   intent(in)    :: atol
      !! absolute solution tolerance
    real, optional,   intent(in)    :: rtol
      !! relative solution tolerance
    real, optional,   intent(in)    :: ftol
      !! absolute residual tolerance
    real, optional,   intent(in)    :: dx_lim
      !! limit for Newton update
    real, optional,   intent(in)    :: dx_lim_rel
      !! limit for Newton update, relative to current value of variable
    real, optional,   intent(in)    :: xmin
      !! lower solution bound
    real, optional,   intent(in)    :: xmax
      !! upper solution bound

    integer :: ivar, itab, ibl, i0, i1

    ivar = this%sys%search_main_var(vs_name)
    do itab = 1, size(this%sys%res2block(ivar)%d)
      ibl = this%sys%res2block(ivar)%d(itab)
      i0  = this%sys%i0(ibl)
      i1  = this%sys%i1(ibl)
      if (present(eabs)) this%eabs(i0:i1) = eabs
      if (present(erel)) this%erel(i0:i1) = erel
      if (present(atol)) this%atol(i0:i1) = atol
      if (present(rtol)) this%rtol(i0:i1) = rtol
      if (present(ftol)) this%ftol(i0:i1) = ftol
      if (present(dx_lim)) this%dx_lim(i0:i1) = dx_lim
      if (present(dx_lim_rel)) this%dx_lim_rel(i0:i1) = dx_lim_rel
      if (present(xmin)) this%xmin(i0:i1) = xmin
      if (present(xmax)) this%xmax(i0:i1) = xmax
    end do
  end subroutine

  subroutine transient_init_output(this, vars, varfile, delta_it, t)
    !! add variable to output
    class(transient),  intent(inout) :: this
    type(string),      intent(in)    :: vars(:)
      !! names of the variables
    character(*),      intent(in)    :: varfile
      !! file to which the variables are written
    integer, optional, intent(in)    :: delta_it
      !! write variables at every delta_it-th step. 1: all steps (default)
    real,    optional, intent(in)    :: t(:)
      !! write variables at these specific time points

    integer            :: i
    type(variable_ptr) :: ptr

    ! check that input is valid
    m4_assert(size(vars) /= 0)
    do i = 1, size(vars)
      ptr = this%sys%search_var(vars(i)%s)
    end do

    ! delete old output settings
    if (allocated(this%output_vars)) deallocate(this%output_vars)
    if (allocated(this%varfile)) deallocate(this%varfile)

    ! new output settings
    this%output_vars = vars
    this%varfile = varfile
    this%vars_delta_it = 1
    if (present(delta_it)) this%vars_delta_it = delta_it
    m4_assert(this%vars_delta_it > 0)
    if (allocated(this%vars_t)) deallocate(this%vars_t)
    allocate (this%vars_t(0))
    if (present(t)) this%vars_t = t
  end subroutine

  subroutine transient_init_cache(this, cachefile, delta_it, t)
    !! init output for solution vector this%sys%get_x()
    class(transient),    intent(inout) :: this
    character(*),        intent(in)    :: cachefile
      !! file to which the variables are written
    integer, optional,   intent(in)    :: delta_it
      !! write cache at every delta_it-th step. 1: all steps (default)
    real,    optional,   intent(in)    :: t(:)
      !! write cache at these specific time points

    integer :: unit, stat

    ! delete old output settings
    if (allocated(this%cachefile)) deallocate(this%cachefile)

    ! new output settings
    this%cachefile = cachefile
    this%cache_delta_it = 1
    if (present(delta_it)) this%cache_delta_it = delta_it
    m4_assert(this%cache_delta_it > 0)
    if (allocated(this%cache_t)) deallocate(this%cache_t)
    allocate (this%cache_t(0))
    if (present(t)) this%cache_t = t
    ! if cachefile already exists, delete it
    open(newunit = unit, file = this%cachefile, iostat = stat, status = "replace")
    close(unit, status="delete")
  end subroutine

  subroutine transient_run(this, times, n_steps, dt0, input, eval_before_output, start_steady_state, start_file, start_i)
    !! perform transient analysis
    class(transient),           intent(inout) :: this
    real,                       intent(in)    :: times(:)
      !! sorted array of time points to be hit exactly, including start and end time
    integer,          optional, intent(in)    :: n_steps(:)
      !! number of time steps for each section of transient simulation, not including the start point (only pass for fixed time steps)
    real,             optional, intent(in)    :: dt0
      !! adaptive: initial time step (default 1fs)
    class(input_src), optional, intent(in)    :: input
      !! optional input source
    logical,          optional, intent(in)    :: eval_before_output
      !! call sys%eval before output to make all variables consistent with the last Newton update (default true)
    logical,          optional, intent(in)    :: start_steady_state
      !! assume that the start values are a steady-state solution (default false)
    character(*),     optional, intent(in)    :: start_file
      !! path to cachefile of old simulation which shall be used as start value. Otherwise this%sys%get_x() is used
    integer,          optional, intent(in)    :: start_i
      !! time step number of desired start solution in start_file

    if (allocated(this%t))       deallocate(this%t)
    if (allocated(this%x))       deallocate(this%x)
    if (allocated(this%n_steps)) deallocate(this%n_steps)

    ! check input parameters
    m4_assert(size(times) > 1)
    m4_assert(all(times(2:size(times))-times(1:size(times)-1) > 0))
    if (this%adaptive) then
      m4_assert(.not. present(n_steps))
      if (present(dt0)) then
        m4_assert(dt0 > 0)
      end if
      allocate(this%n_steps(size(times) - 1))
    else
      m4_assert(present(n_steps))
      m4_assert(size(n_steps) > 0)
      m4_assert(all(n_steps > 0))
      this%n_steps = n_steps
      m4_assert(.not. present(dt0))
      if (size(times) - 1 /= size(n_steps)) call program_error("number of time points and time steps does not match")
    end if
    if (present(input)) then
      m4_assert(input%n == this%sys%ninput)
    else
      m4_assert(this%sys%ninput == 0)
    end if
    if (present(start_file) .neqv. present(start_i)) call program_error("when using an old solution as start, both a file and iteration number have to be provided")

    ! call method-specific subroutines
    select case(this%method)
    case (TRANS_BE)
      call program_error("not implemented yet")
    case (TRANS_TRAPZ)
      call this%solve_trapz(times, dt0 = dt0, input = input, eval_before_output = eval_before_output, &
                           &start_steady_state = start_steady_state, start_file = start_file, start_i = start_i)
    case (TRANS_BDF2)
      call program_error("not implemented yet")
    case (TRANS_MBDF2)
      call program_error("not implemented yet")
    case (TRANS_TRBDF2)
      call this%solve_trbdf2(times, dt0 = dt0, input = input, eval_before_output = eval_before_output, &
                            &start_steady_state = start_steady_state, start_file = start_file, start_i = start_i)
    case default
      call program_error("invalid transient method")
    end select
  end subroutine

  subroutine transient_solve_trapz(this, times, dt0, input, eval_before_output, start_steady_state, start_file, start_i)
    !! trapezoidal method
    class(transient),           intent(inout) :: this
    real,                       intent(in)    :: times(:)
      !! time points which shall be hit exactly, including start and end point
    real,             optional, intent(in)    :: dt0
      !! adaptive: initial time step (default 1fs)
    class(input_src), optional, intent(in)    :: input
      !! optional input source
    logical,          optional, intent(in)    :: eval_before_output
      !! call sys%eval before output to make all variables consistent with the last Newton update (default true)
    logical,          optional, intent(in)    :: start_steady_state
      !! assume that the start values are a steady-state solution (default false)
    character(*),     optional, intent(in)    :: start_file
      !! path to cachefile of old simulation which shall be used as start value. Otherwise this%sys%get_x() is used
    integer,          optional, intent(in)    :: start_i
      !! time step number of desired start solution in start_file

    integer           :: i, p_cntr
    logical           :: convergence, start_steady_state_
    real              :: delta_t
    real, allocatable :: x_old(:), x_new(:), f_old(:), f_new(:), tmp_t(:)
    type(block_real)  :: mat, mat_prec
    type(container)   :: c_input

    ! allocate arrays
    if (this%adaptive) then
      allocate(this%t(0 : 500))
      if (this%use_ram) allocate(this%x(this%sys%n, 0 : 500))
    else
      allocate(this%t(0 : sum(this%n_steps)))
      if (this%use_ram) allocate(this%x(this%sys%n, 0 : sum(this%n_steps)))
    end if
    allocate(x_old(this%sys%n), x_new(this%sys%n), f_old(this%sys%n), f_new(this%sys%n))

    ! set start input
    call this%sys%set_input(input%get(times(1)))
    ! set first delta_t
    if (this%adaptive) then
      delta_t = norm(1.0, "fs")
      if (present(dt0)) delta_t = dt0
    else
      delta_t = (times(2)-times(1)) / this%n_steps(1)
    end if
    ! read start values from file if present, otherwise take current values from esystem
    if (present(start_file)) then
      call c_input%open(start_file, flag = STORAGE_READ)
      call c_input%read("transient/cache/"//int2str(start_i), x_old)
      call c_input%read("transient/t", tmp_t)
      call c_input%close()
      tmp_t = norm(tmp_t, "s")
      if (abs(tmp_t(start_i+1)-times(1))/abs(tmp_t(start_i+1)) > minval(this%rtol)) call program_error("given start time is not the same as time from start solution")
      deallocate (tmp_t)
      call this%sys%set_x(x_old)
    else
      x_old = this%sys%get_x()
    end if
    start_steady_state_ = .false.
    if (present(start_steady_state)) start_steady_state_ = start_steady_state
    if (start_steady_state_) then
      ! f_old = 0.0 if start values are steady-state
      f_old = 0.0
    else
      call this%sys%eval(f = f_old)
    end if
    this%t(0) = times(1)
    if (this%use_ram) this%x(:, 0) = x_old

    ! output start values
    call this%output_grids()
    call this%output(0, eval_before_output)

    ! time loop
    i = 1
    p_cntr = 2
    this%t(1) = this%t(0) + delta_t
    do
      if (this%log) print "(I0, ES25.16, x, A)", i, denorm(this%t(i), "fs"), "fs"

      ! set input at t(i)
      if (present(input)) call this%sys%set_input(input%get(this%t(i)))

      ! trapezoidal rule with previous solution as starting guess
      call this%newton_x(trapz_get_dx, x_old, x_new, convergence)

      if (this%adaptive) then
        call program_error("adaptive time steps not implemented for trapezoidal rule")
      else
        if (.not. convergence) then
          ! Newton not converged: error
          call program_error("Newton not converged for fixed time step")
        else
          ! accept time step
          if (this%use_ram) this%x(:, i) = x_new
          call this%output(i, eval_before_output)
          ! check if this was one of the given time points
          if (abs(this%t(i)-times(p_cntr))/this%t(i) < 1e-10) then
            ! check if end time is reached
            if (p_cntr == size(times)) exit
            p_cntr = p_cntr + 1
          end if
          ! adjust step size
          if (abs(this%t(i)-times(p_cntr-1))/this%t(i) < 1e-10) delta_t = (times(p_cntr)-times(p_cntr-1)) / this%n_steps(p_cntr-1)
          ! advance to next step
          i = i + 1
          f_old = f_new
          x_old = x_new
        end if
      end if
      this%t(i) = this%t(i-1) + delta_t
    end do
    call mat_prec%destruct()
    call mat%destruct()

  contains

    subroutine trapz_get_dx(x, dx, f)
      !! get dx and f
      real, intent(in)  :: x(:)
        !! current solution
      real, intent(out) :: dx(:)
        !! Newton update
      real, intent(out) :: f(:)
        !! residual

      type(block_real), pointer :: dfdx, dfdx_prec, dft
      type(single_matop_real) :: mulvec, precon

      dft => this%sys%dft
      call mat_prec%destruct()
      call mat%destruct()

      if (this%solver == SOLVER_GMRES) then
        ! initial solution
        dx = 0
        ! evaluate function, construct f, Jacobian and preconditioner for trapezoidal rule
        call this%sys%eval(f = f_new, df = dfdx, dfp = dfdx_prec)
        f = 0.5 * delta_t * (f_old + f_new)
        call dft%mul_vec(x - x_old, f, fact_y = 1.0)
        call matrix_add(dft, dfdx,      mat,      fact2 = 0.5 * delta_t)
        call matrix_add(dft, dfdx_prec, mat_prec, fact2 = 0.5 * delta_t)
        ! init gmres
        call mat_prec%factorize(solver = this%gopt%solver)
        call mulvec%init(mat)
        call precon%init(mat_prec, inv = .true.)
        ! solve for dx by gmres
        call gmres(f, mulvec, dx, opts = this%gopt, precon = precon)
      else
        ! evaluate function, construct f and Jacobian for trapezoidal rule
        call this%sys%eval(f = f_new, df = dfdx)
        f = 0.5 * delta_t * (f_old + f_new)
        call dft%mul_vec(x - x_old, f, fact_y = 1.0)
        call matrix_add(dft, dfdx, mat, fact2 = 0.5 * delta_t)
        ! solve for dx
        call mat%factorize(solver = this%solver)
        call mat%solve_vec(f, dx)
      end if
    end subroutine
  end subroutine

  subroutine transient_solve_trbdf2(this, times, dt0, input, eval_before_output, start_steady_state, start_file, start_i)
    !! tr-bdf2
    class(transient),           intent(inout) :: this
    real,                       intent(in)    :: times(:)
      !! time points which shall be hit exactly, including start and end point
    real,             optional, intent(in)    :: dt0
      !! adaptive: initial time step (default 1fs)
    class(input_src), optional, intent(in)    :: input
      !! optional input source
    logical,          optional, intent(in)    :: eval_before_output
      !! call sys%eval before output to make all variables consistent with the last Newton update (default true)
    logical,          optional, intent(in)    :: start_steady_state
      !! assume that the start values are a steady-state solution (default false)
    character(*),     optional, intent(in)    :: start_file
      !! path to cachefile of old simulation which shall be used as start value. Otherwise this%sys%get_x() is used
    integer,          optional, intent(in)    :: start_i
      !! time step number of desired start solution in start_file

    integer           :: istop, i, p_cntr
    logical           :: convergence, start_steady_state_
    real              :: delta_t, r, fact
    real, allocatable :: z_old(:), z_gam(:), z_new(:), x_old(:), x_gam(:), x_new(:), z0(:), tmp_f(:), tmp_t(:), tmp_x(:,:)
    type(block_real)  :: mat, mat_prec
    type(container)   :: c_input

    ! allocate arrays
    if (this%adaptive) then
      allocate(this%t(0 : 500))
      if (this%use_ram) allocate(this%x(this%sys%n, 0 : 500))
    else
      allocate(this%t(0 : sum(this%n_steps)))
      if (this%use_ram) allocate(this%x(this%sys%n, 0 : sum(this%n_steps)))
    end if
    allocate(z_old(this%sys%n), z_gam(this%sys%n), z_new(this%sys%n), x_old(this%sys%n), x_gam(this%sys%n), x_new(this%sys%n))

    ! set start input
    call this%sys%set_input(input%get(times(1)))
    ! set first delta_t
    if (this%adaptive) then
      delta_t = norm(1.0, "fs")
      if (present(dt0)) delta_t = dt0
    else
      delta_t = (times(2)-times(1)) / this%n_steps(1)
    end if
    ! read start values from file if present, otherwise take current values from esystem
    if (present(start_file)) then
      call c_input%open(start_file, flag = STORAGE_READ)
      call c_input%read("transient/cache/"//int2str(start_i), x_old)
      call c_input%read("transient/t", tmp_t)
      call c_input%close()
      tmp_t = norm(tmp_t, "s")
      if (abs(tmp_t(start_i+1)-times(1))/abs(tmp_t(start_i+1)) > minval(this%rtol)) call program_error("given start time is not the same as time from start solution")
      deallocate (tmp_t)
      call this%sys%set_x(x_old)
    else
      x_old = this%sys%get_x()
    end if
    start_steady_state_ = .false.
    if (present(start_steady_state)) start_steady_state_ = start_steady_state
    if (start_steady_state_) then
      ! z_old = 0.0 if start values are steady-state
      z_old = 0.0
    else
      ! Solve M * z_old = delta_t * f_old in the least squares sense. Assume that M is sparse.
      ! For consistent start values, the (possibly underdetermined) system has at least one solution
      allocate (tmp_f(this%sys%n))
      call this%sys%eval(f = tmp_f)
      call this%sys%dft%solve_lsqr(delta_t*tmp_f, 0.0, z_old, istop, atol = 0.0, btol = minval(this%rtol))
      ! istop = 1: solution is found -> norm(b-A*x) <= btol*norm(b) + atol*norm(A*x)
      if (istop > 1) call program_error("No start solution could be found. Probably the given start values are inconsistent")
      ! istop = 0: z = 0 is solution (-> steady-state), but we still need to make sure the system is consistent
      if (istop == 0 .and. any(tmp_f > this%ftol)) call program_error("No start solution could be found. Probably the given start values are inconsistent")
    end if
    this%t(0) = times(1)
    if (this%use_ram) this%x(:, 0) = x_old

    ! output start values
    call this%output_grids()
    call this%output(0, eval_before_output)

    ! time loop
    i = 1
    p_cntr = 2
    this%t(1) = this%t(0) + delta_t
    do
      if (this%log) print "(I0, ES25.16, x, A)", i, denorm(this%t(i), "fs"), "fs"

      ! set input at t_gam
      if (present(input)) call this%sys%set_input(input%get(this%t(i-1) + delta_t * GAM))

      ! trapezoidal stage
      ! initial guess such that x = xold for stability
      z0 = -z_old
      ! M*z_gam - delta_t*f(t,x(z_gam)) = 0 with x(z_gam) = x_old - GAM/2 * (z_old + z_gam)
      call this%newton_z(trbdf2_get_dz, trapz_get_x, trapz_get_z, -GAM/2.0, z0, z_gam, convergence)
      x_gam = trapz_get_x(z_gam)

      if (convergence) then
        ! set input at t(i)
        if (present(input)) call this%sys%set_input(input%get(this%t(i)))

        ! bdf2 stage
        ! initial guess such that x = xgam for stability
        z0 = (GAM-1.0)/2.0 * (z_old + z_gam)
        ! M*z_new - delta_t*f(t,x(z_new)) = 0 with x(z_new) = x_old - 1/(2*(2-GAM)) * (z_old + z_gam) - GAM/2 * z_new
        call this%newton_z(trbdf2_get_dz, bdf2_get_x, bdf2_get_z, -GAM/2.0, z0, z_new, convergence)
        x_new = bdf2_get_x(z_new)
      end if

      if (this%adaptive) then
        if(.not. convergence) then
          ! Newton not converged: reject time step
          fact = this%reduction_factor
        else
          r = trbdf2_error_est()
          if (r > 2) then
            ! error too large: reject time step
            fact = max(0.1, THETA_1 / (r ** POW))
          else
            ! accept time step
            if (this%use_ram) this%x(:, i) = x_new
            call this%output(i, eval_before_output)
            ! check if this was one of the given time points
            if (abs(this%t(i)-times(p_cntr))/this%t(i) < 1e-10) then
              this%n_steps(p_cntr-1) = i - sum(this%n_steps(1:p_cntr-2))
              ! check if end time is reached
              if (p_cntr == size(times)) exit
              p_cntr = p_cntr + 1
            end if
            ! adjust step size
            fact = min(1.2, THETA_1 / (r ** POW))
            ! further adjustment if next point can be hit with a max 10% increase of delta_t
            if (this%t(i) + 1.1 * delta_t * fact > times(p_cntr)) fact = (times(p_cntr) - this%t(i)) / delta_t
            ! advance to next step
            i = i + 1
            call this%time_alloc(i)
            z_old = z_new
            x_old = x_new
          end if
        end if
      else
        if (.not. convergence) then
          ! Newton not converged: error
          call program_error("Newton not converged for fixed time step")
        else
          ! accept time step
          if (this%use_ram) this%x(:, i) = x_new
          call this%output(i, eval_before_output)
          ! check if this was one of the given time points
          if (abs(this%t(i)-times(p_cntr))/this%t(i) < 1e-10) then
            ! check if end time is reached
            if (p_cntr == size(times)) exit
            p_cntr = p_cntr + 1
          end if
          ! adjust step size
          fact = 1.0
          if (abs(this%t(i)-times(p_cntr-1))/this%t(i) < 1e-10) fact = (times(p_cntr)-times(p_cntr-1)) / this%n_steps(p_cntr-1) / delta_t
          ! advance to next step
          i = i + 1
          z_old = z_new
          x_old = x_new
        end if
      end if
      delta_t = delta_t * fact
      z_old = z_old * fact
      this%t(i) = this%t(i-1) + delta_t
    end do
    call mat%destruct()
    call mat_prec%destruct()

    if (this%adaptive) then
      ! truncate the arrays
      allocate(tmp_t(0:i))
      tmp_t(0:i) = this%t(0:i)
      call move_alloc(tmp_t, this%t)
      if (this%use_ram) then
        allocate (tmp_x(this%sys%n, 0:i))
        tmp_x(:,0:i) = this%x(:,0:i)
        call move_alloc(tmp_x, this%x)
      end if
    end if
    m4_assert(size(this%t) == sum(this%n_steps)+1)

  contains

    subroutine trbdf2_get_dz(z, dz, f)
      !! get dz and residual for M*z - delta_t*f(t,x(z)) = 0 with dx/dz = -GAM/2.0
      real, intent(in)  :: z(:)
        !! current solution
      real, intent(out) :: dz(:)
        !! Newton update
      real, intent(out) :: f(:)
        !! residual

      type(block_real), pointer :: dfdx, dfdx_prec, dft
      type(single_matop_real)   :: mulvec, precon

      dft => this%sys%dft
      call mat%destruct()
      call mat_prec%destruct()

      if (this%solver == SOLVER_GMRES) then
        ! initial solution
        dz = 0
        ! evaluate function, construct f, Jacobian and preconditioner for trapezoidal rule
        call this%sys%eval(f = f, df = dfdx, dfp = dfdx_prec)
        call dft%mul_vec(z, f, fact_y = -delta_t)
        call matrix_add(dft, dfdx,      mat,      fact2 = 0.5 * delta_t * GAM)
        call matrix_add(dft, dfdx_prec, mat_prec, fact2 = 0.5 * delta_t * GAM)
        ! init gmres
        call mat_prec%factorize(solver = this%gopt%solver)
        call mulvec%init(mat)
        call precon%init(mat_prec, inv = .true.)
        ! solve for dz by gmres
        call gmres(f, mulvec, dz, opts = this%gopt, precon = precon)
      else
        ! evaluate function, construct f and Jacobian for trapezoidal rule
        call this%sys%eval(f = f, df = dfdx)
        call dft%mul_vec(z, f, fact_y = -delta_t)
        call matrix_add(dft, dfdx, mat, fact2 = 0.5 * delta_t * GAM)
        ! solve for dz
        call mat%factorize(solver = this%solver)
        call mat%solve_vec(f, dz)
      end if
    end subroutine

    function trapz_get_x(z) result(x)
      !! get x from z
      real, intent(in)  :: z(:)
      real, allocatable :: x(:)

      x = x_old - GAM/2.0 * (z_old + z)
    end function

    function trapz_get_z(x) result(z)
      !! get z from x
      real, intent(in)  :: x(:)
      real, allocatable :: z(:)

      z = -z_old + 2.0/GAM * (x_old - x)
    end function

    function bdf2_get_x(z) result(x)
      !! get x from z
      real, intent(in)  :: z(:)
      real, allocatable :: x(:)

      x = x_old - 1.0/(2.0*(2.0-GAM)) * (z_old + z_gam) - GAM/2.0 * z
    end function

    function bdf2_get_z(x) result(z)
      !! get z from x
      real, intent(in)  :: x(:)
      real, allocatable :: z(:)

      z = -1.0/(GAM*(2-GAM)) * (z_old + z_gam) + 2.0/GAM * (x_old - x)
    end function

    function trbdf2_error_est() result(est)
      !! generalized error estimate by Bank (1985) and Hosea (1996): (M+delta_t*gamma/2*J) * err = M * (y* - y)
      real                    :: est

      real, allocatable       :: xmax(:), denom(:), err(:), tmp1(:), tmp2(:)
      type(single_matop_real) :: mulvec, precon

      allocate(err(this%sys%n))
      allocate(tmp2(this%sys%n))

      xmax = max(abs(x_old), abs(x_gam), abs(x_new))
      denom = max(xmax * this%erel, this%eabs)
      tmp1 = -(C1 * z_old + C2 * z_gam + C3 * z_new)
      call this%sys%dft%mul_vec(tmp1, tmp2)
      if (this%solver == SOLVER_GMRES) then
        ! init gmres
        m4_assert(mat_prec%factorized)
        call mulvec%init(mat)
        call precon%init(mat_prec, inv = .true.)
        ! solve for err by gmres
        call gmres(tmp2, mulvec, err, opts = this%gopt, precon = precon)
      else
        ! solve directly
        m4_assert(mat%factorized)
        call mat%solve_vec(tmp2, err)
      end if
      ! use inf norm instead of 2 norm (Matlab)
      est = maxval(abs(err / denom))
    end function

  end subroutine

  subroutine transient_newton_x(this, get_dx, x0, x, convergence)
    !! Newton iteration for solution vector x
    class(transient), intent(inout) :: this
    procedure(get_update)           :: get_dx
      !! procedure to calculate dx and f in Newton iteration, depends on transient method used
    real,             intent(in)    :: x0(:)
      !! initial guess
    real,             intent(out)   :: x(:)
      !! solution
    logical,          intent(out)   :: convergence
      !! convergence status

    character(:), allocatable :: msg_lim
    integer                   :: it
    logical                   :: limmin, limmax
    real                      :: abs_err, dx0
    real,         allocatable :: f(:), xold(:), dx(:), err(:)

    ! allocate temporary arrays for Newton iteration
    allocate (f  (this%sys%n))
    allocate (dx (this%sys%n))
    allocate (err(this%sys%n))

    ! start values for Newton iteration
    it = 0
    err = huge(err)
    f = huge(err)
    x = x0
    convergence = .true.
    ! make sure start solution is within xmin and xmax
    limmin = any(x <= this%xmin)
    limmax = any(x >= this%xmax)
    x = min(this%xmax, max(this%xmin, x))
    call this%sys%set_x(x)

    ! Newton iteration
    if (this%log) m4_info(this%msg // " ITER      REL_ERROR      ABS_ERROR       RESIDUAL")
    if (limmin .and. this%log) m4_info(this%msg // " " // int2str(it, "(I4)") // COL_MAGENTA // ' Solution limited to min' // COL_DEFAULT)
    if (limmax .and. this%log) m4_info(this%msg // " " // int2str(it, "(I4)") // COL_MAGENTA // ' Solution limited to max' // COL_DEFAULT)
    do while ((it < this%min_it) .or. (any(err > this%rtol) .and. any(abs(f) > this%ftol)) .or. limmin .or. limmax)
      it = it + 1

      ! check for maximum number of iterations
      if (it > this%max_it) then
        m4_info("it      = " // int2str(it))
        m4_info("err     = " // real2str(maxval(err)))
        m4_info("abs_err = " // real2str(abs_err))
        if (.not. this%adaptive) then
          call program_error("solution could not be found within maximum number of iterations")
        else
          m4_info("solution could not be found within maximum number of iterations")
          convergence = .false.
          exit
        end if
      end if

      ! calculation of dx and f depends on the transient method
      call get_dx(x, dx, f)

      ! limit update by absolute dx limit
      dx0 = maxval(abs(dx) / this%dx_lim, dim=1)
      if (dx0 > 1) dx = dx / dx0

      ! limit update relative to current variable values
      dx0 = maxval(abs(dx / (x * this%dx_lim_rel)), dim=1)
      ! if x=0 everywhere and dx_lim_rel is not set:
      if (ieee_is_nan(dx0)) dx0 = 0.0
      ! if dx_lim_rel is set at a point where variable is 0:
      if (.not. ieee_is_finite(dx0)) call program_error("dx_lim_rel is set where a variable is 0")
      if (dx0 > 1) dx = dx / dx0

      ! update solution
      xold = x
      x = x - dx

      ! check if new solution violates bounds
      limmin = any(x <= this%xmin)
      limmax = any(x >= this%xmax)
      x = min(this%xmax, max(this%xmin, x))

      ! update esystem
      call this%sys%set_x(x)

      ! calculate new error
      err     = abs(x-xold) / (abs(xold) + this%atol / this%rtol)
      abs_err = maxval(abs(dx))

      if (limmin .and. limmax) then
        msg_lim = COL_MAGENTA // 'Solution limited to min and max' // COL_DEFAULT
      elseif (limmin) then
        msg_lim = COL_MAGENTA // 'Solution limited to min' // COL_DEFAULT
      elseif (limmax) then
        msg_lim = COL_MAGENTA // 'Solution limited to max' // COL_DEFAULT
      else
        msg_lim = ''
      end if
      if (this%log) m4_info(this%msg // " " // int2str(it, "(I4)") // " " // real2str(maxval(err), "(ES14.6E3)") // " " // &
        &                  real2str(abs_err, "(ES14.6E3)") // " " // real2str(maxval(abs(f)), "(ES14.6E3)") // " " // msg_lim)
    end do
  end subroutine

  subroutine transient_newton_z(this, get_dz, get_x, get_z, dxdz, z0, z, convergence)
    !! Newton iteration for scaled derivative z
    class(transient), intent(inout) :: this
    procedure(get_update)           :: get_dz
      !! subroutine to calculate dz and f in Newton iteration, depends on transient method used
    procedure(get_val)              :: get_x
      !! function to update x from current value of z, depends on transient method used
    procedure(get_val)              :: get_z
      !! function to update z from current value of x, depends on transient method used
    real,             intent(in)    :: dxdz
      !! constant derivative of solution x w.r.t. scaled derivative z, depends on transient method used
    real,             intent(in)    :: z0(:)
      !! initial guess
    real,             intent(out)   :: z(:)
      !! solution
    logical,          intent(out)   :: convergence
      !! convergence status

    character(:), allocatable :: msg_lim
    integer                   :: it
    logical                   :: limmin, limmax
    real                      :: abs_err, dz0
    real,         allocatable :: f(:), x(:), zold(:), dz(:), err(:)

    ! allocate temporary arrays for Newton iteration
    allocate (f  (this%sys%n))
    allocate (x  (this%sys%n))
    allocate (dz (this%sys%n))
    allocate (err(this%sys%n))

    ! start values for Newton iteration
    it = 0
    err = huge(err)
    f = huge(err)
    z = z0
    convergence = .true.
    ! make sure start solution is within xmin and xmax
    x = get_x(z)
    limmin = any(x <= this%xmin)
    limmax = any(x >= this%xmax)
    x = min(this%xmax, max(this%xmin, x))
    if (limmin .or. limmax) z = get_z(x)
    call this%sys%set_x(x)

    ! Newton iteration
    if (this%log) m4_info(this%msg // " ITER      REL_ERROR      ABS_ERROR       RESIDUAL")
    if (limmin .and. this%log) m4_info(this%msg // " " // int2str(it, "(I4)") // COL_MAGENTA // ' Solution limited to min' // COL_DEFAULT)
    if (limmax .and. this%log) m4_info(this%msg // " " // int2str(it, "(I4)") // COL_MAGENTA // ' Solution limited to max' // COL_DEFAULT)
    do while ((it < this%min_it) .or. (any(err > this%rtol) .and. any(abs(f) > this%ftol)) .or. limmin .or. limmax)
      it = it + 1

      ! check for maximum number of iterations
      if (it > this%max_it) then
        m4_info("it      = " // int2str(it))
        m4_info("err     = " // real2str(maxval(err)))
        m4_info("abs_err = " // real2str(abs_err))
        if (.not. this%adaptive) then
          call program_error("solution could not be found within maximum number of iterations")
        else
          m4_info("solution could not be found within maximum number of iterations")
          convergence = .false.
          exit
        end if
      end if

      ! calculation of dz and f depends on the transient method
      call get_dz(z, dz, f)

      ! limit update by absolute dz limit (derived from dx limit via dxdz)
      dz0 = maxval(abs(dz*dxdz) / this%dx_lim, dim=1)
      if (dz0 > 1) dz = dz / dz0

      ! limit update relative to current variable values
      dz0 = maxval(abs(dz*dxdz / (x * this%dx_lim_rel)), dim=1)
      ! if x=0 everywhere and dx_lim_rel is not set:
      if (ieee_is_nan(dz0)) dz0 = 0.0
      ! if dx_lim_rel is set at a point where variable is 0:
      if (.not. ieee_is_finite(dz0)) call program_error("dx_lim_rel is set where a variable is 0")
      if (dz0 > 1) dz = dz / dz0

      ! update solution
      zold = z
      z = z - dz

      ! check if new solution violates bounds
      x = get_x(z)
      limmin = any(x <= this%xmin)
      limmax = any(x >= this%xmax)
      x = min(this%xmax, max(this%xmin, x))
      if (limmin .or. limmax) z = get_z(x)

      ! update esystem
      call this%sys%set_x(x)

      ! calculate new error
      err     = abs(z-zold) / (abs(zold) + this%atol / this%rtol)
      abs_err = maxval(abs(dz))

      if (limmin .and. limmax) then
        msg_lim = COL_MAGENTA // 'Solution limited to min and max' // COL_DEFAULT
      elseif (limmin) then
        msg_lim = COL_MAGENTA // 'Solution limited to min' // COL_DEFAULT
      elseif (limmax) then
        msg_lim = COL_MAGENTA // 'Solution limited to max' // COL_DEFAULT
      else
        msg_lim = ''
      end if
      if (this%log) m4_info(this%msg // " " // int2str(it, "(I4)") // " " // real2str(maxval(err), "(ES14.6E3)") // " " // &
        &                  real2str(abs_err, "(ES14.6E3)") // " " // real2str(maxval(abs(f)), "(ES14.6E3)") // " " // msg_lim)
    end do
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
      if (this%use_ram) then
        allocate (tmp_x(this%sys%n, l:2*n))
        tmp_x(:,l:u) = this%x(:,l:u)
        call move_alloc(tmp_x, this%x)
      end if
    end if
  end subroutine

  subroutine transient_output_grids(this)
    !! output grids into varfile
    class(transient), intent(inout) :: this

    integer            :: i, stat
    type(container)    :: c_vars
    type(variable_ptr) :: ptr

    if (allocated(this%varfile)) then
      call c_vars%open(this%varfile, flag = STORAGE_WRITE)
      ! output the grids on which the output variables are defined
      do i = 1, size(this%output_vars)
        ptr = this%sys%search_var(this%output_vars(i)%s)
        ! save grid if it is not already saved
        call c_vars%save(ptr%p%g, stat = stat)
        if (stat /= 0 .and. stat /= ERR_ALREADY_EXISTS) call program_error("Error when saving grids")
      end do
      call c_vars%close()
    end if
  end subroutine

  subroutine transient_output(this, i, eval_before_output)
    !! output values at every delta_it-th iteration or if the current time matches any of the output times exactly
    class(transient),  intent(inout) :: this
    integer,           intent(in)    :: i
      !! time point index
    logical, optional, intent(in)    :: eval_before_output
      !! call sys%eval before output to make all variables consistent with the last Newton update (default true)

    integer            :: j
    logical            :: eval_before_output_
    type(container)    :: c_vars, c_cache
    type(variable_ptr) :: ptr

    eval_before_output_ = .true.
    if (present(eval_before_output)) eval_before_output_ = eval_before_output

    if (allocated(this%varfile)) then
      if (mod(i, this%vars_delta_it)==0 .or. any((this%vars_t-this%t(i))/this%t(i) < 1e-10)) then
        call c_vars%open(this%varfile, flag = STORAGE_WRITE)
        call c_vars%write("transient/t", [this%t(i)], unit = "s", dynamic = DYNAMIC_APP)
        if (eval_before_output_) call this%sys%eval()
        do j = 1, size(this%output_vars)
          ptr = this%sys%search_var(this%output_vars(j)%s)
          select type(pp => ptr%p)
          class is(variable_real)
            call c_vars%save(pp, "transient", dynamic = .true.)
          end select
        end do
        call c_vars%close()
      end if
    end if
    ! output solution vector
    if (allocated(this%cachefile)) then
      if (mod(i, this%cache_delta_it)==0 .or. any((this%cache_t-this%t(i))/this%t(i) < 1e-10)) then
        call c_cache%open(this%cachefile, flag = STORAGE_WRITE)
        call c_cache%write("transient/t", [this%t(i)], unit = "s", dynamic = DYNAMIC_APP)
        call c_cache%write("transient/cache/"//int2str(i), this%sys%get_x())
        call c_cache%close()
      end if
    end if
  end subroutine

  subroutine transient_select_i(this, i, cachefile, eval)
    !! save steady-state result at step i (overwrite esystem variables)
    class(transient),       intent(inout) :: this
    integer,                intent(in)    :: i
      !! step index
    character(*), optional, intent(in)    :: cachefile
      !! file where the system is stored. If not present the result is taken from this%x
    logical,      optional, intent(in)    :: eval
      !! evaluate the esystem after setting the main variables (default false)

    real, allocatable :: x(:)
    type(container)   :: c_cache

    ! store solution vector in esystem
    if (present(cachefile)) then
      call c_cache%open(cachefile, flag = STORAGE_READ)
      call c_cache%read("transient/cache/"//int2str(i), x)
      call c_cache%close()
      call this%sys%set_x(x)
    elseif (allocated(this%x)) then
      call this%sys%set_x(this%x(:,i))
    elseif (allocated(this%cachefile)) then
      call c_cache%open(this%cachefile, flag = STORAGE_READ)
      call c_cache%read("transient/cache/"//int2str(i), x)
      call c_cache%close()
      call this%sys%set_x(x)
    else
      call program_error("Neither this%x nor this%cachefile are allocated")
    end if

    ! evaluate esystem to update its other variables
    if (present(eval)) then
      if (eval) call this%sys%eval()
    end if
  end subroutine

  subroutine transient_select_t(this, t, cachefile, eval)
    !! save steady-state result at step i (overwrite esystem variables)
    class(transient),       intent(inout) :: this
    real,                   intent(in)    :: t
      !! step index
    character(*), optional, intent(in)    :: cachefile
      !! file where the system is stored. If not present the result is taken from this%x
    logical,      optional, intent(in)    :: eval
      !! evaluate the esystem after setting the main variables (default false)

    integer           :: i
    real, allocatable :: times(:), x(:)
    type(container)   :: c_cache

    ! store solution vector in esystem
    if (present(cachefile)) then
      call c_cache%open(cachefile, flag = STORAGE_READ)
      call c_cache%read("transient/t", times)
      i = bin_search(times, t)
      call c_cache%read("transient/cache/"//int2str(i), x)
      call c_cache%close()
      call this%sys%set_x(x)
    elseif (allocated(this%x)) then
      i = bin_search(this%t, t)
      call this%sys%set_x(this%x(:,i))
    elseif (allocated(this%cachefile)) then
      call c_cache%open(this%cachefile, flag = STORAGE_READ)
      call c_cache%read("transient/t", times)
      i = bin_search(times, t)
      call c_cache%read("transient/cache/"//int2str(i), x)
      call c_cache%close()
      call this%sys%set_x(x)
    else
      call program_error("Neither this%x nor this%cachefile are allocated")
    end if

    ! evaluate esystem to update its other variables
    if (present(eval)) then
      if (eval) call this%sys%eval()
    end if
  end subroutine

end module