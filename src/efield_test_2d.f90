program efield_test_2d

  use approx_m,           only: approx_imref, approx_potential
  use cl_options_m,       only: cl_option_descriptor, cl_option, get_cl_options
  use device_m,           only: dev
  use error_m,            only: program_error
  use grid_m,             only: IDX_VERTEX
  use harmonic_balance_m, only: harmonic_balance
  use input_m,            only: input_file
  use input_src_m,        only: polygon_src, harmonic_src
  use math_m,             only: linspace, logspace, PI
  use matrix_m,           only: block_real
  use newton_m,           only: newton_opt
  use normalization_m,    only: init_normconst, norm, denorm
  use semiconductor_m,    only: CR_NAME
  use small_signal_m,     only: small_signal
  use steady_state_m,     only: steady_state
  use storage_m,          only: storage, STORAGE_WRITE, DYNAMIC_EXT, DYNAMIC_APP
  use string_m,           only: string, new_string
  use transient_m,        only: transient, TRANS_TRBDF2
  use util_m,             only: get_hostname

  use current_density_m, only: time

  implicit none

  logical          :: gummel_restart, gummel_once, gummel_enabled
  real             :: temperature
  type(input_file) :: runfile

  print "(A)", "Start simulation on " // get_hostname()
  print *, "2D Electric Field Verification Test"

  ! parse command line arguments
  call command_line()

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
    integer              :: si
    integer, allocatable :: sids(:)
    logical              :: log
    real,    allocatable :: V(:,:), t_inp(:), t(:)
    type(string)         :: name
    type(polygon_src)    :: input
    type(steady_state)   :: ss

    call runfile%get_sections("steady state", sids)
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
      call ss%input_newton_params(runfile, "full newton params")
      call ss%input_var_params(runfile, "full newton params")
      print "(A,I0)", "DEBUG: Using solver = ", ss%solver
      print "(A,I0,A,I0)", "DEBUG: g%dim = ", dev%par%g%dim, ", g%idx_dim = ", dev%par%g%idx_dim
      ! Use LEFT and RIGHT contact variables for testing, include Ey for 2D
      call ss%init_output([new_string("pot"), new_string("ndens"), new_string("Ex"), new_string("Ey"), &
        & new_string("V_LEFT"), new_string("I_LEFT"), new_string("I_RIGHT")], name%s // ".fbs")
      call ss%run(input = input, t_input = t, gummel = gummel)

      ! Perform analytical comparison after steady-state
      call verify_electric_field()
    end do
  end subroutine

  subroutine verify_electric_field()
    !! Compare computed electric field with analytical solution

    integer :: i, n_points, ict, ix, iy, nx, ny
    real    :: V_applied, L_device, E_analytical
    real    :: E_computed, E_error, max_error, avg_error
    real    :: p_vertex(dev%par%g%dim), ct_surf_val, tr_vol_val
    real, allocatable :: efield_x(:), efield_y(:), x_positions(:), y_positions(:), pot_values(:), ct_surf_values(:), tr_vol_values(:)
    integer, allocatable :: idx_temp(:)
    character(len=10), allocatable :: vertex_type(:)

    print *, ""
    print *, "========================================="
    print *, "2D Electric Field Analytical Verification"
    print *, "========================================="

    ! Get device parameters (assuming uniform 1D grid from 0 to 100 nm)
    V_applied = norm(0.15, "V")  ! Normalize the applied voltage
    L_device = norm(100.0, "nm")  ! Normalize the device length

    ! Analytical solution for uniform field: E = -V/L (both in normalized units)
    E_analytical = -V_applied / L_device  ! This gives normalized E-field

    ! Get computed electric field values in grid order
    n_points = dev%efield(1)%data%n
    allocate(efield_x(n_points), efield_y(n_points), pot_values(n_points))
    efield_x = dev%efield(1)%get()  ! Ex component
    if (dev%par%g%dim >= 2) then
      efield_y = dev%efield(2)%get()  ! Ey component
    else
      efield_y = 0.0
    end if
    pot_values = dev%pot%get()

    ! Debug: Print potential values to check ordering
    print *, ""
    print *, "DEBUG: Potential values in data array order:"
    do i = 1, dev%pot%data%n
      print "(A,I3,A,ES15.8)", "  pot[", i, "] = ", denorm(pot_values(i), "V")
    end do

    ! Allocate arrays for position, type info, and contact surface
    allocate(x_positions(n_points), y_positions(n_points), vertex_type(n_points), ct_surf_values(n_points), tr_vol_values(n_points))
    allocate(idx_temp(dev%par%g%idx_dim))

    ! For 2D tensor grid, iterate through vertices systematically
    ! The data is stored in column-major order (Fortran style)
    nx = dev%par%g1D(1)%n  ! number of x points
    ny = dev%par%g1D(2)%n  ! number of y points

    i = 0
    do iy = 1, ny
      do ix = 1, nx
        i = i + 1
        if (i > n_points) exit
        idx_temp = [ix, iy]

        ! Get vertex position
        call dev%par%g%get_vertex(idx_temp, p_vertex)
        x_positions(i) = denorm(p_vertex(1), "nm")
        if (dev%par%g%dim >= 2) then
          y_positions(i) = denorm(p_vertex(2), "nm")
        else
          y_positions(i) = 0.0
        end if

        ! Get contact index for this vertex
        ict = dev%par%ict%get(idx_temp)

        ! Always try to get tr_vol if it exists (for both contact and interior vertices)
        if (dev%par%transport(IDX_VERTEX,0)%flags%get(idx_temp)) then
          tr_vol_values(i) = dev%par%tr_vol%get(idx_temp)
        else
          tr_vol_values(i) = 0.0
        end if

        ! Get values for comparison
        if (ict > 0) then
          ! Contact vertex: get contact boundary length [nm in 2D]
          ct_surf_val = dev%par%get_ct_surf(ict, idx_temp)
          ct_surf_values(i) = ct_surf_val
          vertex_type(i) = trim(dev%par%contacts(ict)%name)
        else
          ! Interior vertex: set to 0 since we display tr_vol separately
          ct_surf_values(i) = 0.0
          vertex_type(i) = "interior"
        end if
      end do
      if (i > n_points) exit
    end do

    ! Calculate error metrics
    max_error = 0.0
    avg_error = 0.0

    print *, ""
    print "(A)", "All grid points:"
    print "(A)", "Idx  Grid[x,y]  Type     X(nm)  Y(nm)   Pot(V)    Ex(V/cm)    Ey(V/cm)   ct_surf     tr_vol"
    print "(A)", "---  ---------  -------  -----  -----   -------   ---------   ---------  ----------  ----------"

    ! Re-iterate for printing (same ordering as above)
    i = 0
    do iy = 1, ny
      do ix = 1, nx
        i = i + 1
        if (i > n_points) exit
        idx_temp = [ix, iy]

        ! Print all points with both ct_surf and tr_vol
        print "(I3, 2X, '[',I2,',',I2,']', 2X, A8, 2F7.1, 3ES12.4, 2ES12.4)", i, idx_temp(1), idx_temp(2), vertex_type(i), &
          & x_positions(i), y_positions(i), denorm(pot_values(i), "V"), &
          & denorm(efield_x(i), "V/cm"), denorm(efield_y(i), "V/cm"), &
          & denorm(ct_surf_values(i), "nm"), denorm(tr_vol_values(i), "nm^2")
      end do
      if (i > n_points) exit
    end do

    ! Calculate statistics for Ex and Ey separately
    print *, ""
    print *, "Field Statistics:"
    print "(A,ES12.5,A)", "Analytical Ex: ", denorm(E_analytical, "V/cm"), " V/cm"
    print "(A,ES12.5,A)", "Average Ex:    ", denorm(sum(efield_x)/real(n_points), "V/cm"), " V/cm"
    print "(A,ES12.5,A)", "Average Ey:    ", denorm(sum(efield_y)/real(n_points), "V/cm"), " V/cm (should be ~0)"
    print "(A,ES12.5,A)", "Max |Ey|:      ", denorm(maxval(abs(efield_y)), "V/cm"), " V/cm"

    ! Check field uniformity in interior
    print *, ""
    print *, "Interior field uniformity check:"
    block
      real :: ex_min, ex_max, ey_min, ey_max
      integer :: n_interior
      n_interior = 0
      ex_min = huge(1.0)
      ex_max = -huge(1.0)
      ey_min = huge(1.0)
      ey_max = -huge(1.0)

      do i = 1, n_points
        if (vertex_type(i) == "interior") then
          n_interior = n_interior + 1
          ex_min = min(ex_min, efield_x(i))
          ex_max = max(ex_max, efield_x(i))
          ey_min = min(ey_min, efield_y(i))
          ey_max = max(ey_max, efield_y(i))
        end if
      end do

      if (n_interior > 0) then
        print "(A,I5)", "Number of interior points: ", n_interior
        print "(A,ES12.5,A,ES12.5,A)", "Ex range: [", denorm(ex_min, "V/cm"), ", ", denorm(ex_max, "V/cm"), "] V/cm"
        print "(A,ES12.5,A,ES12.5,A)", "Ey range: [", denorm(ey_min, "V/cm"), ", ", denorm(ey_max, "V/cm"), "] V/cm"
        print "(A,ES12.5,A)", "Ex variation: ", denorm(ex_max - ex_min, "V/cm"), " V/cm"
        print "(A,ES12.5,A)", "Ey variation: ", denorm(ey_max - ey_min, "V/cm"), " V/cm"
      end if
    end block

    deallocate(efield_x, efield_y, vertex_type, x_positions, y_positions, idx_temp, ct_surf_values, tr_vol_values)
  end subroutine

  subroutine gummel()
    !! gummel iteration

    integer            :: it, ci, dir, min_it, max_it
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

    do ci = dev%par%ci0, dev%par%ci1
      call runfile%get("dd params", "log", log)
      call ss_dd(ci)%init(dev%sys_dd(ci), log = log, msg = CR_NAME(ci) // "DD: ")
      call ss_dd(ci)%input_newton_params(runfile, "dd params")
      call ss_dd(ci)%input_var_params(runfile, "dd params")
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
        call dev%calc_efield(dir)%eval()
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

    integer                   :: it, nx, min_it, max_it
    real                      :: err, res0, res1, damping, dx0, atol, dx_lim
    real, allocatable         :: x0(:), f(:), dx(:)
    type(block_real), pointer :: dfdx

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
    call runfile%get("nlpe params", "atol", atol)
    call runfile%get("nlpe params", "dx_lim", dx_lim)
    call runfile%get("nlpe params", "min_it", min_it)
    call runfile%get("nlpe params", "max_it", max_it)

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
      call dfdx%factorize()
      call dfdx%solve_vec(-f, dx)
      call dfdx%reset(only_factorization = .true.)

      ! limit update
      dx0 = maxval(abs(dx) / dx_lim, dim=1)
      if (dx0 > 1) dx = dx / dx0

      ! absolute error
      err = maxval(abs(dx))

      ! update variables
      call dev%sys_nlpe%set_x(x0 + dx)

      write (*, "(ES25.16E3)") denorm(err, "V")
    end do
  end subroutine

end program efield_test_2D
