program efield_test
  !! Electric Field Verification Test (1D)
  !!
  !! Tests the electric_field module by comparing computed E-field
  !! against analytical solution for a simple 1D device with uniform field.

  use approx_m,          only: approx_imref, approx_potential
  use block_m,           only: block_real
  use cl_options_m,      only: cl_option_descriptor, cl_option, get_cl_options
  use device_m,          only: dev
  use electric_field_m,  only: electric_field, calc_efield
  use error_m,           only: program_error
  use input_m,           only: input_file, input_section
  use input_src_m,       only: polygon_src
  use math_m,            only: linspace
  use normalization_m,   only: init_normconst, norm, denorm
  use semiconductor_m,   only: CR_NAME
  use solver_base_m,     only: solver_real
  use solver_m,          only: default_solver_params, init_solver_real
  use steady_state_m,    only: steady_state
  use string_m,          only: string
  use util_m,            only: get_hostname

  implicit none

  logical          :: gummel_restart, gummel_once, gummel_enabled
  real             :: temperature
  type(input_file) :: runfile

  ! Electric field variables (local to test program)
  type(electric_field), allocatable :: efield(:)
  type(calc_efield),    allocatable :: calc_ef(:)

  print "(A)", "Start simulation on " // get_hostname()
  print *, "Electric Field Verification Test (1D)"

  ! parse command line arguments
  call command_line()

  ! initialize electric field after device is ready
  call init_electric_field()

  ! solve
  call solve_steady_state()

contains

  subroutine command_line()
    !! parse command line arguments, init normalization, device and run file

    integer                      :: idesc, i, j
    integer,         allocatable :: iclopt(:), jclopt(:)
    type(cl_option), allocatable :: clopt(:)
    type(cl_option_descriptor)   :: desc(3) = [ &
          cl_option_descriptor('T', "temperature", .true., .false., .true., .true.), &
          cl_option_descriptor('r', "run",         .true., .false., .true., .true.), &
          cl_option_descriptor('d', "device",      .true., .false., .true., .true.)]

    call get_cl_options(desc, clopt, iclopt, jclopt)

    temperature = 300 ! default

    do idesc = 1, size(desc)
      do i = iclopt(idesc), iclopt(idesc + 1) - 1
        j = jclopt(i)
        select case (clopt(j)%short)
      case ('T')
        read (clopt(j)%arg, *) temperature
        call init_normconst(temperature)
        print "(A,ES25.16E3,A)", "T = ", temperature, " K"

      case ('d')
        call dev%init(clopt(j)%arg, temperature)

      case ('r')
        call runfile%init(clopt(j)%arg)

      end select
      end do
    end do
  end subroutine

  subroutine init_electric_field()
    !! Initialize electric field variables and equations
    integer :: dir

    ! Allocate for all spatial dimensions
    allocate(efield(dev%par%g%dim))
    allocate(calc_ef(dev%par%g%dim))

    ! Initialize each component
    do dir = 1, dev%par%g%dim
      call efield(dir)%init(dev%par, dir)
      call calc_ef(dir)%init(dev%par, dev%pot, efield(dir))
    end do

    print "(A,I0,A)", "Initialized electric field for ", dev%par%g%dim, " dimension(s)"
  end subroutine

  subroutine voltage_input_ss(sid, t_inp, V, t_sim)
    !! read in voltage configuration of contacts for steady-state from runfile
    integer,           intent(in)  :: sid
    !! section id for runfile
    real, allocatable, intent(out) :: t_inp(:)
    !! pseudo time points for polygon input source
    real, allocatable, intent(out) :: V(:,:)
    !! voltages at each pseudo time point (dev%par%nct, size(t_inp))
    real, allocatable, intent(out) :: t_sim(:)
    !! pseudo time points at which steady-state simulations are performed

    integer              :: ict, ict_sweep, i
    integer, allocatable :: nsweep_ct(:), nsweep(:)
    logical              :: status
    real,    allocatable :: Vbounds(:), tmp(:)

    ! find out which, if any, contact is swept and how many points are used
    ict_sweep = 0
    do ict = 1, dev%par%nct
      call runfile%get(sid, "N_"//dev%par%contacts(ict)%name, nsweep_ct, status = status)
      if (status) status = (size(nsweep_ct) > 1 .or. nsweep_ct(1) > 1)
      if (status) then
        if (ict_sweep /= 0) call program_error("only one voltage sweep per steady-state section allowed")
        ict_sweep = ict
        nsweep = nsweep_ct
      end if
    end do
    if (.not. allocated(nsweep)) allocate (nsweep(0))

    ! t_inp
    t_inp = linspace(0.0, real(size(nsweep)), size(nsweep)+1)

    ! contact voltages
    allocate(V(dev%par%nct, size(nsweep)+1))
    do ict = 1, dev%par%nct
      call runfile%get(sid, "V_"//dev%par%contacts(ict)%name, Vbounds)
      if (ict == ict_sweep) then
        V(ict,:) = Vbounds
      else
        if (size(Vbounds) > 1) call program_error("Voltage bounds given for a contact without sweep")
        V(ict,:) = Vbounds(1)
      end if
    end do

    ! t_sim
    if (ict_sweep == 0) then
      t_sim = [0.0]
    else
      t_sim = [0.0]
      do i = 1, size(nsweep)
        tmp = linspace(i-1.0, i+0.0, nsweep(i))
        t_sim = [t_sim, tmp(2:size(tmp))]
      end do
    end if
  end subroutine

  subroutine solve_steady_state()
    integer              :: si, sj, dir
    integer, allocatable :: sids(:)
    logical              :: log
    real,    allocatable :: V(:,:), t_inp(:), t(:)
    type(string)         :: name
    type(polygon_src)    :: input
    type(steady_state)   :: ss

    call runfile%get_sections("steady state", sids)
    call runfile%get_section("full newton params", sj)
    do si = 1, size(sids)
      print "(A)", "steady state"

      call runfile%get(sids(si), "name", name)

      ! input config
      call voltage_input_ss(sids(si), t_inp, V, t)
      call input%init(t_inp, V)

      gummel_restart = .true.
      gummel_once    = .false.
      gummel_enabled = .true.

      ! solve steady-state
      call runfile%get("full newton params", "log", log)
      call ss%init(dev%sys_full, log = log, msg = "Newton: ")
      call ss%set_params(runfile%sections%d(sj))
      print "(A,I0,A,I0)", "DEBUG: g%dim = ", dev%par%g%dim, ", g%idx_dim = ", dev%par%g%idx_dim

      ! Setup output variables
      call ss%init_output([string("pot"), string("ndens"), &
        & string("V_LEFT"), string("I_LEFT"), string("I_RIGHT")], name%s // ".fbs")
      call ss%run(input = input, t_input = t, gummel = gummel)

      ! Update electric field from converged potential
      do dir = 1, dev%par%g%dim
        call calc_ef(dir)%eval()
      end do

      ! Perform analytical comparison after steady-state
      call verify_electric_field()
    end do
  end subroutine

  subroutine verify_electric_field()
    !! Compare computed electric field with analytical solution

    integer :: i, n_points, ict
    real    :: V_left, V_right, L_device, E_analytical
    real    :: E_computed, E_error, max_error, avg_error
    real    :: p_vertex(dev%par%g%dim)
    real, allocatable :: efield_values(:), x_positions(:), pot_values(:)
    integer, allocatable :: idx_temp(:)
    character(len=16), allocatable :: vertex_type(:)

    print *, ""
    print *, "========================================="
    print *, "Electric Field Analytical Verification"
    print *, "========================================="

    ! Get applied voltages from contacts (use %x for scalar access)
    V_left  = dev%volt(1)%x
    V_right = dev%volt(2)%x
    print "(A,F10.4,A)", "V_LEFT  = ", denorm(V_left, "V"), " V"
    print "(A,F10.4,A)", "V_RIGHT = ", denorm(V_right, "V"), " V"

    ! Get device length from grid
    n_points = efield(1)%data%n
    allocate(idx_temp(dev%par%g%idx_dim))

    ! Find device length from first and last vertices
    idx_temp = [1]
    call dev%par%g%get_vertex(idx_temp, p_vertex)
    L_device = p_vertex(1)

    idx_temp = [n_points]
    call dev%par%g%get_vertex(idx_temp, p_vertex)
    L_device = p_vertex(1) - L_device

    print "(A,F10.4,A)", "L_device = ", denorm(L_device, "nm"), " nm"

    ! Analytical solution for uniform field: E = -(V_right - V_left) / L
    E_analytical = -(V_right - V_left) / L_device

    ! Get computed electric field and potential values
    allocate(efield_values(n_points), pot_values(n_points))
    allocate(x_positions(n_points), vertex_type(n_points))

    efield_values = efield(1)%get()
    pot_values = dev%pot%get()

    ! Get positions and classify vertices
    do i = 1, n_points
      idx_temp = [i]
      call dev%par%g%get_vertex(idx_temp, p_vertex)
      x_positions(i) = denorm(p_vertex(1), "nm")

      ! Get contact index for this vertex
      ict = dev%par%ict%get(idx_temp)
      if (ict > 0) then
        vertex_type(i) = dev%par%contacts(ict)%name
      else
        vertex_type(i) = "interior"
      end if
    end do

    ! Calculate error metrics
    max_error = 0.0
    avg_error = 0.0

    print *, ""
    print "(A)", "Idx  Type          X(nm)      Potential(V)     E_computed(V/cm)    E_analytical(V/cm)    Error(%)"
    print "(A)", "---  -----------   -------    ------------     ----------------    -----------------    --------"

    do i = 1, n_points
      E_computed = efield_values(i)

      if (abs(E_analytical) > 1e-20) then
        E_error = abs((E_computed - E_analytical) / E_analytical) * 100.0
      else
        E_error = 0.0
      end if

      print "(I3, 2X, A12, F10.2, 4ES20.6)", i, trim(vertex_type(i)), &
        & x_positions(i), denorm(pot_values(i), "V"), denorm(E_computed, "V/cm"), &
        & denorm(E_analytical, "V/cm"), E_error

      max_error = max(max_error, abs(E_computed - E_analytical))
      avg_error = avg_error + abs(E_computed - E_analytical)
    end do

    avg_error = avg_error / real(n_points)

    print *, ""
    print "(A,ES15.6,A)", "Analytical E-field:  ", denorm(E_analytical, "V/cm"), " V/cm"
    print "(A,ES15.6,A)", "Average computed E:  ", denorm(sum(efield_values)/real(n_points), "V/cm"), " V/cm"
    print "(A,ES15.6,A)", "Max absolute error:  ", denorm(max_error, "V/cm"), " V/cm"
    print "(A,ES15.6,A)", "Avg absolute error:  ", denorm(avg_error, "V/cm"), " V/cm"

    ! Check interior points uniformity (skip boundaries)
    if (n_points > 4) then
      print *, ""
      print *, "Interior field uniformity check (points 3 to n-2):"
      max_error = 0.0
      do i = 3, n_points-2
        max_error = max(max_error, abs(efield_values(i) - efield_values(3)))
      end do
      print "(A,ES15.6,A)", "Max variation in interior: ", denorm(max_error, "V/cm"), " V/cm"
    end if

    ! Summary
    print *, ""
    print *, "========================================="
    if (denorm(avg_error, "V/cm") < 1.0) then
      print *, "TEST PASSED: Average error < 1 V/cm"
    else
      print *, "TEST FAILED: Average error >= 1 V/cm"
    end if
    print *, "========================================="

    deallocate(efield_values, pot_values, x_positions, vertex_type, idx_temp)
  end subroutine

  subroutine gummel()
    !! gummel iteration

    integer            :: it, ci, dir, min_it, max_it, si
    logical            :: log
    real               :: error, err_pot, err_iref(2), atol
    real, allocatable  :: pot0(:), iref0(:,:)
    type(steady_state) :: ss_dd(2)

    if (.not. gummel_enabled) return
    if (gummel_once) gummel_enabled = .false.

    allocate (pot0(dev%pot%data%n), iref0(dev%iref(1)%data%n,2))

    if (gummel_restart) then
      ! approximate imrefs
      do ci = dev%par%ci0, dev%par%ci1
        call approx_imref(dev%par, dev%iref(ci), dev%volt)
      end do

      ! approximate potential
      call approx_potential(dev%par, dev%pot, dev%iref)
    end if

    call runfile%get_section("dd params", si)
    do ci = dev%par%ci0, dev%par%ci1
      call runfile%get("dd params", "log", log)
      call ss_dd(ci)%init(dev%sys_dd(ci), log = log, msg = CR_NAME(ci) // "DD: ")
      call ss_dd(ci)%set_params(runfile%sections%d(si))
    end do

    ! get gummel params
    call runfile%get("gummel params", "log", log)
    call runfile%get("gummel params", "atol", atol)
    call runfile%get("gummel params", "min_it", min_it)
    call runfile%get("gummel params", "max_it", max_it)

    ! gummel iteration
    it = 0
    error = huge(1.0)
    do while (((error > atol) .and. (it < max_it)) .or. (it < min_it))
      it = it + 1

      ! solve non-linear poisson equation
      pot0 = dev%pot%get()
      call solve_nlpe()
      err_pot = maxval(abs(dev%pot%get() - pot0))
      error = err_pot

      ! update electric field from new potential
      do dir = 1, dev%par%g%dim
        call calc_ef(dir)%eval()
      end do

      ! solve dd model for electrons and holes
      do ci = dev%par%ci0, dev%par%ci1
        iref0(:,ci) = dev%iref(ci)%get()

        call ss_dd(ci)%run()
        call ss_dd(ci)%select(1)
        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:,ci)))
        error = max(error, err_iref(ci))
      end do

      ! log
      if (log) then
        print "(A,I6,ES25.16E3)", "Gummel: ", it, denorm(error, "V")
      end if
    end do

    if ((it > max_it) .and. (error > atol)) then
      call program_error("Gummel iteration did not converge (maximum number of iterations reached)")
    end if
  end subroutine

  subroutine solve_nlpe()
    !! solve non-linear poisson equation

    integer                         :: it, nx, min_it, max_it, si
    logical                         :: status
    real                            :: err, res0, res1, damping, dx0, atol, dx_lim
    real,               allocatable :: x0(:), f(:), dx(:)
    type(block_real),   pointer     :: dfdx
    type(string)                    :: solver_name
    type(input_section)             :: solver_params, solver_params_tmp
    class(solver_real), allocatable :: solver

    ! memory
    nx = dev%sys_nlpe%n
    allocate (x0(nx), f(nx), dx(nx), source = 0.0)

    ! init iteration params
    it   = 0
    err  = huge(err)
    res0 = huge(res0)
    damping = 1.0

    ! initial approximation
    x0 = dev%sys_nlpe%get_x()

    ! get nlpe params
    call runfile%get_section("nlpe params", si)
    call runfile%get(si, "atol", atol)
    call runfile%get(si, "dx_lim", dx_lim)
    call runfile%get(si, "min_it", min_it)
    call runfile%get(si, "max_it", max_it)

    ! solver name
    call runfile%get(si, "solver", solver_name, status = status)
    if (.not. status) solver_name%s = "pardiso"

    ! solver parameters
    solver_params = default_solver_params(solver_name%s)
    call runfile%sections%d(si)%get_subsection(solver_name%s, solver_params_tmp)
    if (solver_params_tmp%params%n > 0) call solver_params%set_subsection(solver_params_tmp)

    ! init solver
    call init_solver_real(solver_name%s, solver_params, solver)

    ! newton iteration
    do while (((err > atol) .and. (it <= max_it)) .or. (it < min_it))
      it = it + 1

      ! evaluate system
      call dev%sys_nlpe%eval(f = f, df = dfdx)
      res1 = dot_product(f, f)
      write (*, "(A,I6,2ES25.16E3)", advance = "no") "NLPE: ", it, res1, damping

      ! repeat step with 0.5*dx if residual gets larger
      if (res1 > res0) then
        it = it - 1
        damping = damping * 0.5
        call dev%sys_nlpe%set_x(x0 + damping * dx)
        print *
        cycle
      end if

      ! accept step
      x0      = dev%sys_nlpe%get_x()
      res0    = res1
      damping = 1.0

      ! solve
      call solver%factorize(dfdx)
      call solver%solve(-f, dx)

      ! limit update
      dx0 = maxval(abs(dx) / dx_lim, dim=1)
      if (dx0 > 1) dx = dx / dx0

      ! absolute error
      err = maxval(abs(dx))

      ! update variables
      call dev%sys_nlpe%set_x(x0 + dx)

      write (*, "(ES25.16E3)") denorm(err, "V")
    end do

    call solver%destruct()
  end subroutine

end program
