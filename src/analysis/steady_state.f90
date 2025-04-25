m4_include(../util/macro.f90.inc)

module steady_state_m

  use color_m,         only: COL_DEFAULT, COL_MAGENTA
  use container_m,     only: container, STORAGE_READ, STORAGE_WRITE, ERR_ALREADY_EXISTS
  use error_m,         only: assert_failed, program_error
  use esystem_m,       only: esystem
  use gmres_m,         only: gmres_options, gmres
  use ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_positive_inf, ieee_negative_inf
  use input_src_m,     only: input_src
  use logging_m,       only: logging
  use matop_m,         only: single_matop_real
  use matrix_m,        only: block_real, SOLVER_GMRES, SPSOLVER_PARDISO
  use string_m,        only: string, new_string
  use util_m,          only: int2str, real2str
  use variable_m,      only: variable_ptr, variable_real

  implicit none

  private
  public steady_state

  type steady_state
    !! steady-state simulation -> solve stationary equation from esystem: f(x) = 0

    ! general data
    type(esystem), pointer    :: sys => null()
      !! pointer to corresponding equation system
    real, allocatable         :: x(:,:)
      !! steady-state data (sys%n, size(t_input)), only used if results are stored in RAM memory
    logical                   :: use_ram
      !! save data also in RAM memory (true, default) or only on disk (false)
    logical                   :: log
      !! enable/disable logging (default false)
    character(:), allocatable :: msg
      !! message to print each Newton iteration if logging is enabled (default "Steady-State: ")
    
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
    logical             :: error_if_not_converged
      !! error if one step does not converge (default true)
    logical             :: converged_when_lim
      !! convergence can be achieved even if the solution is limited to xmin or xmax (default false)

    ! data for output (default no output)
    type(string), allocatable :: output_vars(:)
      !! list of variables for output
    integer                   :: vars_delta_it
      !! write variables at each vars_delta_it-th step
    character(:), allocatable :: varfile
      !! output file for variables
    integer                   :: cache_delta_it
      !! save this%system%get_x() at each cache_delta_it-th step
    character(:), allocatable :: cachefile
      !! output file for cache
  contains
    procedure :: init              => steady_state_init
    procedure :: set_newton_params => steady_state_set_newton_params
    procedure :: set_var_params    => steady_state_set_var_params
    procedure :: init_output       => steady_state_init_output
    procedure :: init_cache        => steady_state_init_cache
    procedure :: run               => steady_state_run
    procedure :: select            => steady_state_select
  end type

  abstract interface
    subroutine gummel_proc()
    end subroutine
  end interface

contains

  subroutine steady_state_init(this, sys, use_ram, log, msg)
    !! initialize steady-state simulation
    class(steady_state),    intent(out) :: this
    type(esystem), target,  intent(in)  :: sys
      !! equation system
    logical,      optional, intent(in)  :: use_ram
      !! save data also in RAM memory (true, default) or only on disk (false)
    logical,      optional, intent(in)  :: log
      !! enable/disable logging (default false)
    character(*), optional, intent(in)  :: msg
      !! message to print each Newton iteration if logging is enabled (default "Steady-State: ")

    this%sys => sys
    this%use_ram = .true.
    if (present(use_ram)) this%use_ram = use_ram
    this%log = .false.
    if (present(log)) this%log = log
    this%msg = "Steady-State: "
    if (present(msg)) deallocate (this%msg)
    if (present(msg)) this%msg = msg

    this%solver = SPSOLVER_PARDISO
    allocate (this%atol(sys%n), source = 1e-16)
    allocate (this%rtol(sys%n), source = 1e-12)
    allocate (this%ftol(sys%n), source = 0.0)
    allocate (this%dx_lim(sys%n), source = ieee_value(1.0, ieee_positive_inf))
    allocate (this%dx_lim_rel(sys%n), source = ieee_value(1.0, ieee_positive_inf))
    allocate (this%xmin(sys%n), source = ieee_value(1.0, ieee_negative_inf))
    allocate (this%xmax(sys%n), source = ieee_value(1.0, ieee_positive_inf))
    this%min_it = 0
    this%max_it = huge(this%max_it)
    this%error_if_not_converged = .true.
    this%converged_when_lim = .false.
  end subroutine

  subroutine steady_state_set_newton_params(this, solver, gopt, atol, rtol, ftol, dx_lim, dx_lim_rel, xmin, xmax, min_it, max_it, error_if_not_converged, converged_when_lim)
    !! set parameters for Newton iteration, overwriting the whole array
    class(steady_state),           intent(inout) :: this
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
    logical,             optional, intent(in)    :: error_if_not_converged
      !! error if one step does not converge
    logical,             optional, intent(in)    :: converged_when_lim
      !! convergence can be achieved even if the solution is limited to xmin or xmax (default false)

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
    if (present(error_if_not_converged)) this%error_if_not_converged = error_if_not_converged
    if (present(converged_when_lim)) this%converged_when_lim = converged_when_lim
  end subroutine

  subroutine steady_state_set_var_params(this, vs_name, atol, rtol, ftol, dx_lim, dx_lim_rel, xmin, xmax)
    !! set parameters for Newton iteration for specific vselector, leaving the rest of the array unchanged
    class(steady_state), intent(inout) :: this
    character(*),        intent(in)    :: vs_name
      !! name of vselector
    real, optional,      intent(in)    :: atol
      !! absolute solution tolerance
    real, optional,      intent(in)    :: rtol
      !! relative solution tolerance
    real, optional,      intent(in)    :: ftol
      !! absolute residual tolerance
    real, optional,      intent(in)    :: dx_lim
      !! limit for Newton update
    real, optional,      intent(in)    :: dx_lim_rel
      !! limit for Newton update, relative to current value of variable
    real, optional,      intent(in)    :: xmin
      !! lower solution bound
    real, optional,      intent(in)    :: xmax
      !! upper solution bound

    integer :: ivar, itab, ibl, i0, i1

    ivar = this%sys%search_main_var(vs_name)
    do itab = 1, size(this%sys%res2block(ivar)%d)
      ibl = this%sys%res2block(ivar)%d(itab)
      i0  = this%sys%i0(ibl)
      i1  = this%sys%i1(ibl)
      if (present(atol)) this%atol(i0:i1) = atol
      if (present(rtol)) this%rtol(i0:i1) = rtol
      if (present(ftol)) this%ftol(i0:i1) = ftol
      if (present(dx_lim)) this%dx_lim(i0:i1) = dx_lim
      if (present(dx_lim_rel)) this%dx_lim_rel(i0:i1) = dx_lim_rel
      if (present(xmin)) this%xmin(i0:i1) = xmin
      if (present(xmax)) this%xmax(i0:i1) = xmax
    end do
  end subroutine

  subroutine steady_state_init_output(this, vars, varfile, delta_it)
    !! add variable to output
    class(steady_state), intent(inout) :: this
    type(string),        intent(in)    :: vars(:)
      !! names of the variables
    character(*),        intent(in)    :: varfile
      !! file to which the variables are written
    integer, optional,   intent(in)    :: delta_it
      !! write variables at every delta_it-th step. 1: all steps (default), 0: all steps including start values, -1: last step

    integer            :: i
    type(variable_ptr) :: ptr

    ! check that input is valid
    m4_assert(size(vars) /= 0)
    do i = 1, size(vars)
      ptr = this%sys%search_var(vars(i)%s)
    end do

    ! delete old output settings
    if (allocated(this%output_vars)) deallocate (this%output_vars)
    if (allocated(this%varfile)) deallocate (this%varfile)

    ! new output settings
    this%output_vars = vars
    this%varfile = varfile
    this%vars_delta_it = 1
    if (present(delta_it)) this%vars_delta_it = delta_it
    m4_assert(this%vars_delta_it >= -1)
  end subroutine

  subroutine steady_state_init_cache(this, cachefile, delta_it)
    !! init output for solution vector this%sys%get_x()
    class(steady_state), intent(inout) :: this
    character(*),        intent(in)    :: cachefile
      !! file to which the variables are written
    integer, optional,   intent(in)    :: delta_it
      !! write cache at every delta_it-th step. 1: all steps (default), 0: all steps including start values, -1: last step

    ! delete old output settings
    if (allocated(this%cachefile)) deallocate (this%cachefile)

    ! new output settings
    this%cachefile = cachefile
    this%cache_delta_it = 1
    if (present(delta_it)) this%cache_delta_it = delta_it
    m4_assert(this%cache_delta_it >= -1)
  end subroutine

  subroutine steady_state_run(this, input, t_input, gummel, eval_before_output)
    !! perform steady-state analysis for one or multiple sets of input parameters
    class(steady_state),              intent(inout) :: this
    class(input_src),       optional, intent(in)    :: input
      !! optional input source
    real,                   optional, intent(in)    :: t_input(:)
      !! perform quasi-stationary simulations for these time-points (default: [0.0])
    procedure(gummel_proc), optional                :: gummel
      !! optional procedure to run before each newton iteration
    logical,                optional, intent(in)    :: eval_before_output
      !! call sys%eval before output to make all variables consistent with the last Newton update (default true)

    character(:), allocatable :: msg_lim
    integer                   :: i, j, it, nt, stat
    logical                   :: eval_before_output_, limmin, limmax
    real                      :: abs_err, dx0
    real,         allocatable :: t_input_(:), x(:), dx(:), xold(:), f(:), err(:)
    type(block_real), pointer :: dfdx, dfdx_prec
    type(container)           :: c_vars, c_cache
    type(single_matop_real)   :: mulvec, precon
    type(variable_ptr)        :: ptr

    eval_before_output_ = .true.
    if (present(eval_before_output)) eval_before_output_ = eval_before_output

    if (allocated(this%x)) deallocate (this%x)

    ! input
    if (present(input)) then
      m4_assert(input%n == this%sys%ninput)
      if (present(t_input)) then
        t_input_ = t_input
      else
        t_input_ = [0.0]
      end if
      nt = size(t_input_)
    else
      m4_assert(this%sys%ninput == 0)
      m4_assert(.not. present(t_input))
      nt = 1
    end if

    ! output start values
    if (allocated(this%varfile)) then
      call c_vars%open(this%varfile, flag = STORAGE_WRITE)
      ! output the grids on which the output variables are defined
      do j = 1, size(this%output_vars)
        ptr = this%sys%search_var(this%output_vars(j)%s)
        ! save grid if it is not already saved
        call c_vars%save(ptr%p%g, stat = stat)
        if (stat /= 0 .and. stat /= ERR_ALREADY_EXISTS) call program_error("Error when saving grids")
      end do
      ! output variables
      if (this%vars_delta_it == 0) then
        if (eval_before_output_) call this%sys%eval()
        do j = 1, size(this%output_vars)
          ptr = this%sys%search_var(this%output_vars(j)%s)
          select type(pp => ptr%p)
          class is(variable_real)
            call c_vars%save(pp, "steady-state", dynamic = .true.)
          end select
        end do
      end if
      call c_vars%close()
    end if
    ! output solution vector
    if (allocated(this%cachefile)) then
      if (this%cache_delta_it == 0) then
        call c_cache%open(this%cachefile, flag = STORAGE_WRITE)
        call c_cache%write("steady-state/cache/0", this%sys%get_x())
        call c_cache%close()
      end if
    end if

    ! allocate result array if data are stored in RAM
    if (this%use_ram) allocate (this%x(this%sys%n, nt))

    ! allocate temporary arrays for Newton iteration
    allocate (f  (this%sys%n))
    allocate (dx (this%sys%n))
    allocate (err(this%sys%n))

    ! solve steady-state for each input time
    do i = 1, nt
      if (this%log) print *, "steady-state step " // int2str(i) // " of " // int2str(nt)
      
      if (present(input))  call this%sys%set_input(input%get(t_input_(i)))
      if (present(gummel)) call gummel()

      ! start values for Newton iteration
      it = 0
      err = huge(err)
      f = huge(err)
      x = this%sys%get_x()
      ! make sure start solution is within xmin and xmax
      limmin = any(x <= this%xmin)
      limmax = any(x >= this%xmax)
      x = min(this%xmax, max(this%xmin, x))
      call this%sys%set_x(x)

      ! Newton iteration
      if (this%log) m4_info(this%msg // " ITER      REL_ERROR      ABS_ERROR       RESIDUAL")
      if (limmin .and. this%log) m4_info(this%msg // " " // int2str(it, "(I4)") // COL_MAGENTA // ' Solution limited to min' // COL_DEFAULT)
      if (limmax .and. this%log) m4_info(this%msg // " " // int2str(it, "(I4)") // COL_MAGENTA // ' Solution limited to max' // COL_DEFAULT)
      do while ((it < this%min_it) .or. (any(err > this%rtol) .and. any(abs(f) > this%ftol)) .or. ((.not. this%converged_when_lim) .and. (limmin .or. limmax)))
        it = it + 1

        ! check for maximum number of iterations
        if (it > this%max_it) then
          m4_info("it      = " // int2str(it))
          m4_info("err     = " // real2str(maxval(err)))
          m4_info("abs_err = " // real2str(abs_err))
          if (this%error_if_not_converged) then
            call program_error("solution could not be found within maximum number of iterations")
          else
            m4_info("solution could not be found within maximum number of iterations")
            exit
          end if
        end if

        if (this%solver == SOLVER_GMRES) then
          ! initial solution
          dx = 0
          ! evaluate function, construct f, Jacobian and preconditioner
          call this%sys%eval(f = f, df = dfdx, dfp = dfdx_prec)
          ! init gmres
          call dfdx_prec%factorize(solver = this%gopt%solver)
          call mulvec%init(dfdx)
          call precon%init(dfdx_prec, inv = .true.)
          ! solve for dx by gmres
          call gmres(f, mulvec, dx, opts = this%gopt, precon = precon)
          call dfdx_prec%reset(only_factorization = .true.)
        else
          ! evaluate function, construct f and Jacobian
          call this%sys%eval(f = f, df = dfdx)
          ! solve for dx
          call dfdx%factorize(solver = this%solver)
          call dfdx%solve_vec(f, dx)
          call dfdx%reset(only_factorization = .true.)
        end if

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
        if (this%use_ram) this%x(:,i) = x

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

      ! output variables
      if (allocated(this%varfile)) then
        out = .false.
        if (this%vars_delta_it > 0) out = mod(i, this%vars_delta_it)==0
        if (this%vars_delta_it == 0 .or. i == nt) out = .true.
        if (out) then
          call c_vars%open(this%varfile, flag = STORAGE_WRITE)
          if (eval_before_output_) call this%sys%eval()
          do j = 1, size(this%output_vars)
            ptr = this%sys%search_var(this%output_vars(j)%s)
            select type(pp => ptr%p)
            class is(variable_real)
              call c_vars%save(pp, "steady-state", dynamic = .true.)
            end select
          end do
          call c_vars%close()
        end if
      end if
      ! output solution vector
      if (allocated(this%cachefile)) then
        out = .false.
        if (this%vars_delta_it > 0) out = mod(i, this%vars_delta_it)==0
        if (this%vars_delta_it == 0 .or. i == nt) out = .true.
        if (out) then
          call c_cache%open(this%cachefile, flag = STORAGE_WRITE)
          call c_cache%write("steady-state/cache/"//int2str(i), this%sys%get_x())
          call c_cache%close()
        end if
      end if

    end do
  end subroutine

  subroutine steady_state_select(this, i, cachefile, eval)
    !! save steady-state result at step i (overwrite esystem variables)
    class(steady_state),    intent(inout) :: this
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
      call c_cache%read("steady-state/cache/"//int2str(i), x)
      call c_cache%close()
      call this%sys%set_x(x)
    elseif (allocated(this%x)) then
      call this%sys%set_x(this%x(:,i))
    elseif (allocated(this%cachefile)) then
      call c_cache%open(this%cachefile, flag = STORAGE_READ)
      call c_cache%read("steady-state/cache/"//int2str(i), x)
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
