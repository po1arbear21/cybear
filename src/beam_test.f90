program beam_test
  !! EBIC Beam Generation Test Runner

  use device_m,        only: device
  use normalization_m, only: init_normconst, denorm, norm
  use semiconductor_m, only: CR_ELEC, CR_HOLE
  use input_src_m,     only: polygon_src
  use steady_state_m,  only: steady_state
  use string_m,        only: string
  use approx_m,        only: approx_imref, approx_potential
  use block_m,         only: block_real
  use solver_base_m,   only: solver_real
  use solver_m,        only: default_solver_params, init_solver_real
  use input_m,         only: input_section

  implicit none

  type(device) :: dev
  type(polygon_src) :: input
  type(steady_state) :: ss, ss_dd(2)
  class(solver_real), allocatable :: nlpe_solver
  integer :: ix, iy, nx, ny, ict, ninput, ict_n, ci, itest, n_sweep
  real :: G, beam_y, G_sum, G_tot, A, rel_err, I_N, y_start, y_end, dy_sweep
  real, allocatable :: t_inp(:), V(:,:)
  real, allocatable :: sweep_y(:), sweep_I(:)

  ! Initialize
  call init_normconst(300.0)
  call dev%init("beam_minimal.ini", 300.0)

  nx = size(dev%par%g1D(1)%x)
  ny = size(dev%par%g1D(2)%x)

  ! Test 1: Grid overview
  print "(A)", "Test 1: Grid structure"
  print "(A)", ""
  print "(A,I0,A,I0,A,I0)", "  Grid: ", nx, " x ", ny, " = ", nx*ny, " vertices"
  print "(A)", ""

  ! X positions
  print "(A)", "  X grid (thickness direction):"
  print "(A)", "  ┌───────┬────────────┐"
  print "(A)", "  │  ix   │   x [nm]   │"
  print "(A)", "  ├───────┼────────────┤"
  do ix = 1, nx
    print "(A,I4,A,F10.1,A)", "  │", ix, "   │", denorm(dev%par%g1D(1)%x(ix), 'nm'), "  │"
  end do
  print "(A)", "  └───────┴────────────┘"
  print "(A)", ""

  ! Y positions
  print "(A)", "  Y grid (length direction):"
  print "(A)", "  ┌───────┬────────────┐"
  print "(A)", "  │  iy   │   y [nm]   │"
  print "(A)", "  ├───────┼────────────┤"
  do iy = 1, ny
    print "(A,I4,A,F10.1,A)", "  │", iy, "   │", denorm(dev%par%g1D(2)%x(iy), 'nm'), "  │"
  end do
  print "(A)", "  └───────┴────────────┘"
  print "(A)", ""

  ! 2D grid map
  print "(A)", "  2D Grid Map:"
  print "(A)", ""

  ! Header row with y values
  write(*, "(A)", advance='no') "        x\\y │"
  do iy = 1, ny
    write(*, "(I6)", advance='no') nint(denorm(dev%par%g1D(2)%x(iy), 'nm'))
  end do
  print *

  ! Separator
  write(*, "(A)", advance='no') "      ──────┼"
  do iy = 1, ny
    write(*, "(A)", advance='no') "──────"
  end do
  print *

  ! Grid rows
  do ix = 1, nx
    write(*, "(F10.0,A)", advance='no') denorm(dev%par%g1D(1)%x(ix), 'nm'), " │"
    do iy = 1, ny
      write(*, "(A)", advance='no') "     ·"
    end do
    print *
  end do

  ! Test 2: Beam generation values
  print "(A)", ""
  print "(A)", "Test 2: bgen_n values"
  print "(A)", ""

  ! Set beam at y = 2500 nm (middle of device)
  beam_y = 1750.0
  dev%beam_pos%x = norm(beam_y, 'nm')
  print "(A,F8.1,A)", "  Beam position: y = ", beam_y, " nm"
  print "(A)", ""

  ! Evaluate beam generation
  call dev%calc_bgen(CR_ELEC)%eval()

  ! Header row
  write(*, "(A)", advance='no') "      x\\y │"
  do iy = 1, ny
    write(*, "(I10)", advance='no') nint(denorm(dev%par%g1D(2)%x(iy), 'nm'))
  end do
  print *

  ! Separator
  write(*, "(A)", advance='no') "     ──────┼"
  do iy = 1, ny
    write(*, "(A)", advance='no') "──────────"
  end do
  print *

  ! Grid rows with bgen values
  do ix = 1, nx
    write(*, "(F10.0,A)", advance='no') denorm(dev%par%g1D(1)%x(ix), 'nm'), " │"
    do iy = 1, ny
      G = denorm(dev%bgen(CR_ELEC)%x2(ix, iy), '1/cm^3/s')
      if (G > 0.0) then
        write(*, "(ES10.2)", advance='no') G
      else
        write(*, "(A)", advance='no') "         0"
      end if
    end do
    print *
  end do

  ! Test 3: Total generation rate validation
  print "(A)", ""
  print "(A)", "Test 3: Conservation check"
  print "(A)", ""

  ! tr_vol table
  print "(A)", "  tr_vol [cm^2]:"
  print "(A)", ""

  ! Header row
  write(*, "(A)", advance='no') "      x\\y │"
  do iy = 1, ny
    write(*, "(I10)", advance='no') nint(denorm(dev%par%g1D(2)%x(iy), 'nm'))
  end do
  print *

  ! Separator
  write(*, "(A)", advance='no') "     ──────┼"
  do iy = 1, ny
    write(*, "(A)", advance='no') "──────────"
  end do
  print *

  ! Grid rows with tr_vol values
  do ix = 1, nx
    write(*, "(F10.0,A)", advance='no') denorm(dev%par%g1D(1)%x(ix), 'nm'), " │"
    do iy = 1, ny
      A = denorm(dev%par%tr_vol%get([ix, iy]), 'cm^2')
      write(*, "(ES10.2)", advance='no') A
    end do
    print *
  end do
  print "(A)", ""

  G_tot = dev%calc_bgen(CR_ELEC)%G_tot

  G_sum = 0.0
  do ix = 1, nx
    do iy = 1, ny
      G = denorm(dev%bgen(CR_ELEC)%x2(ix, iy), '1/cm^3/s')
      A = denorm(dev%par%tr_vol%get([ix, iy]), 'cm^2')
      G_sum = G_sum + G * A
    end do
  end do

  rel_err = abs(G_sum - G_tot) / G_tot * 100.0

  print "(A,ES12.4,A)", "  G_tot      = ", G_tot, " pairs/(s·cm)"
  print "(A,ES12.4,A)", "  sum(G × A) = ", G_sum, " pairs/(s·cm)"
  print "(A,F8.4,A)",   "  Error      = ", rel_err, " %"
  if (rel_err < 1.0) then
    print "(A)", "  PASS"
  else
    print "(A)", "  FAIL"
  end if

  ! Test 4: Terminal current from DD solver
  print "(A)", ""
  print "(A)", "Test 4: Terminal current (DD solve)"
  print "(A)", ""

  ! Set up input: nct voltages (0V) + beam position
  ninput = dev%par%nct + 1
  allocate(t_inp(2), V(ninput, 2))

  t_inp = [0.0, 1.0]  ! two time points for polygon interpolation

  ! Set all voltages to 0V (short circuit)
  do ict = 1, dev%par%nct
    V(ict, :) = 0.0
  end do

  ! Set beam position (last row)
  V(ninput, :) = dev%beam_pos%x  ! already set from Test 2

  call input%init(t_inp, V)

  ! Set input values at t=0
  do ict = 1, dev%par%nct
    dev%volt(ict)%x = V(ict, 1)
  end do

  ! Step 1: Approximate initial conditions
  print "(A)", "  Step 1: Approximating initial conditions..."
  do ci = CR_ELEC, CR_HOLE
    call approx_imref(dev%par, dev%iref(ci), dev%volt)
  end do
  call approx_potential(dev%par, dev%pot, dev%iref)

  ! Step 2: Solve NLPE (non-linear Poisson equation)
  print "(A)", "  Step 2: Solving NLPE..."
  call solve_nlpe()

  ! Step 3: Gummel iteration
  print "(A)", "  Step 3: Running Gummel iteration..."
  call gummel()

  ! Step 4: Full Newton solver
  print "(A)", "  Step 4: Running full Newton solver..."
  call ss%init(dev%sys_full, log = .true., msg = "Newton: ")
  call ss%run(input = input, t_input = [0.0])

  ! Get N_CONTACT terminal current
  call dev%par%contact_map%get(string("N_CONTACT"), ict_n)
  I_N = denorm(dev%curr(ict_n)%x, 'A')

  print "(A)", ""
  print "(A,F8.1,A)", "  Beam position: y = ", beam_y, " nm"
  print "(A,ES12.4,A)", "  I_N_CONTACT (raw norm) = ", dev%curr(ict_n)%x, " (normalized)"
  print "(A,ES12.4,A)", "  I_N_CONTACT = ", I_N, " A"
  print "(A,ES12.4,A)", "  I_N_CONTACT = ", I_N * 1e9, " nA"

  ! Expected current: I = q * G_tot (ideal, 100% collection)
  print "(A)", ""
  print "(A,ES12.4,A)", "  I_expected  = ", 1.602e-19 * G_tot * 1e-4, " A (ideal)"

  ! Test 5: Fine beam sweep to detect kinks at grid nodes
  print "(A)", ""
  print "(A)", "========================================"
  print "(A)", "Test 5: Fine Beam Sweep (check for kinks at grid nodes)"
  print "(A)", "========================================"
  print "(A)", ""
  print "(A)", "Sweeping beam finely to detect discontinuities at grid nodes."
  print "(A)", "Grid nodes are at: 0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000 nm"
  print "(A)", ""

  ! Fine sweep around y=1500 nm (junction - grid node between p and n regions)
  y_start = 1300.0
  y_end = 1700.0
  dy_sweep = 25.0
  n_sweep = nint((y_end - y_start) / dy_sweep) + 1

  allocate(sweep_y(n_sweep), sweep_I(n_sweep))

  print "(A)", "  Sweeping from 1300 to 1700 nm in 25 nm steps (around junction at 1500 nm)..."
  print "(A)", ""

  do itest = 1, n_sweep
    beam_y = y_start + (itest - 1) * dy_sweep
    sweep_y(itest) = beam_y
    dev%beam_pos%x = norm(beam_y, 'nm')
    V(ninput, :) = dev%beam_pos%x

    call input%init(t_inp, V)

    do ict = 1, dev%par%nct
      dev%volt(ict)%x = V(ict, 1)
    end do

    do ci = CR_ELEC, CR_HOLE
      call approx_imref(dev%par, dev%iref(ci), dev%volt)
    end do
    call approx_potential(dev%par, dev%pot, dev%iref)

    call solve_nlpe_quiet()
    call gummel_quiet()
    call ss%init(dev%sys_full, log = .false., msg = "")
    call ss%run(input = input, t_input = [0.0])

    sweep_I(itest) = denorm(dev%curr(ict_n)%x, 'nA')
  end do

  ! Print results in table format for easy plotting
  print "(A)", "  Results (y [nm], I [nA], dI/dy [nA/nm]):"
  print "(A)", "  ─────────────────────────────────────────"
  do itest = 1, n_sweep
    if (itest == 1) then
      print "(A,F8.1,A,F12.6,A,A)", "  ", sweep_y(itest), "  ", sweep_I(itest), "      ", "---"
    else
      print "(A,F8.1,A,F12.6,A,F10.6)", "  ", sweep_y(itest), "  ", sweep_I(itest), "      ", &
        (sweep_I(itest) - sweep_I(itest-1)) / dy_sweep
    end if
  end do

  ! Check for kink at y=1500 (junction grid node)
  ! Find index closest to 1500
  print "(A)", ""
  print "(A)", "  Checking for kink at y=1500 nm (junction grid node):"
  do itest = 1, n_sweep
    if (abs(sweep_y(itest) - 1500.0) < 1.0) then
      if (itest > 1 .and. itest < n_sweep) then
        print "(A,F10.6,A)", "    dI/dy before: ", (sweep_I(itest) - sweep_I(itest-1)) / dy_sweep, " nA/nm"
        print "(A,F10.6,A)", "    dI/dy after:  ", (sweep_I(itest+1) - sweep_I(itest)) / dy_sweep, " nA/nm"
        print "(A,F10.6,A)", "    Jump in slope: ", &
          abs((sweep_I(itest+1) - sweep_I(itest)) - (sweep_I(itest) - sweep_I(itest-1))) / dy_sweep, " nA/nm"
      end if
      exit
    end if
  end do

  deallocate(sweep_y, sweep_I)

  deallocate(t_inp, V)

contains

  subroutine solve_nlpe()
    !! Solve non-linear Poisson equation using Newton iteration
    integer :: it, n_nlpe
    real :: err, res0, res1, damping, dx0
    real, allocatable :: x0(:), f(:), dx(:)
    type(block_real), pointer :: dfdx
    type(input_section) :: solver_params

    ! Parameters (hardcoded for test)
    real, parameter :: ATOL = 1e-10
    real, parameter :: DX_LIM = 1.0
    integer, parameter :: MAX_IT = 100

    n_nlpe = dev%sys_nlpe%n
    allocate(x0(n_nlpe), f(n_nlpe), dx(n_nlpe), source = 0.0)

    ! Initialize solver if not already done
    if (.not. allocated(nlpe_solver)) then
      solver_params = default_solver_params("pardiso")
      call init_solver_real("pardiso", solver_params, nlpe_solver)
    end if

    it = 0
    err = huge(err)
    res0 = huge(res0)
    damping = 1.0
    x0 = dev%sys_nlpe%get_x()

    do while ((err > ATOL) .and. (it <= MAX_IT))
      it = it + 1

      ! Evaluate system
      call dev%sys_nlpe%eval(f = f, df = dfdx)
      res1 = dot_product(f, f)

      ! Damping if residual increases
      if (res1 > res0) then
        it = it - 1
        damping = damping * 0.5
        call dev%sys_nlpe%set_x(x0 + damping * dx)
        if (damping < 1e-10) exit
        cycle
      end if

      ! Accept step
      x0 = dev%sys_nlpe%get_x()
      res0 = res1
      damping = 1.0

      ! Solve linear system
      call nlpe_solver%factorize(dfdx)
      call nlpe_solver%solve(-f, dx)

      ! Limit update
      dx0 = maxval(abs(dx) / DX_LIM)
      if (dx0 > 1) dx = dx / dx0

      err = maxval(abs(dx))

      ! Update
      call dev%sys_nlpe%set_x(x0 + dx)

      print "(A,I4,A,ES10.3,A,ES10.3)", "    NLPE it ", it, ": err = ", denorm(err, 'V'), ", res = ", res1
    end do

    deallocate(x0, f, dx)
  end subroutine

  subroutine gummel()
    !! Gummel iteration: alternate between NLPE and DD solves
    integer :: it, ci
    real :: err, err_pot, err_iref(2)
    real, allocatable :: pot0(:), iref0(:,:)

    ! Parameters (hardcoded for test)
    real, parameter :: ATOL = 1e-6
    integer, parameter :: MAX_IT = 50

    allocate(pot0(dev%pot%data%n), iref0(dev%iref(1)%data%n, 2))

    ! Initialize DD solvers
    do ci = CR_ELEC, CR_HOLE
      call ss_dd(ci)%init(dev%sys_dd(ci), log = .false., msg = "")
    end do

    it = 0
    err = huge(err)

    do while ((denorm(err, 'V') > ATOL) .and. (it < MAX_IT))
      it = it + 1

      ! Save old values
      pot0 = dev%pot%get()
      do ci = CR_ELEC, CR_HOLE
        iref0(:, ci) = dev%iref(ci)%get()
      end do

      ! Solve NLPE
      call solve_nlpe()
      err_pot = maxval(abs(dev%pot%get() - pot0))

      ! Solve DD for each carrier
      do ci = CR_ELEC, CR_HOLE
        call ss_dd(ci)%run()
        call ss_dd(ci)%select(1)
        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:, ci)))
      end do

      err = max(err_pot, err_iref(1), err_iref(2))
      print "(A,I3,A,ES10.3)", "    Gummel it ", it, ": err = ", denorm(err, 'V')
    end do

    deallocate(pot0, iref0)
  end subroutine

  subroutine solve_nlpe_quiet()
    !! Solve NLPE without printing (for batch tests)
    integer :: it, n_nlpe
    real :: err, res0, res1, damping, dx0
    real, allocatable :: x0(:), f(:), dx(:)
    type(block_real), pointer :: dfdx
    type(input_section) :: solver_params

    real, parameter :: ATOL = 1e-10
    real, parameter :: DX_LIM = 1.0
    integer, parameter :: MAX_IT = 100

    n_nlpe = dev%sys_nlpe%n
    allocate(x0(n_nlpe), f(n_nlpe), dx(n_nlpe), source = 0.0)

    if (.not. allocated(nlpe_solver)) then
      solver_params = default_solver_params("pardiso")
      call init_solver_real("pardiso", solver_params, nlpe_solver)
    end if

    it = 0
    err = huge(err)
    res0 = huge(res0)
    damping = 1.0
    x0 = dev%sys_nlpe%get_x()

    do while ((err > ATOL) .and. (it <= MAX_IT))
      it = it + 1
      call dev%sys_nlpe%eval(f = f, df = dfdx)
      res1 = dot_product(f, f)

      if (res1 > res0) then
        it = it - 1
        damping = damping * 0.5
        call dev%sys_nlpe%set_x(x0 + damping * dx)
        if (damping < 1e-10) exit
        cycle
      end if

      x0 = dev%sys_nlpe%get_x()
      res0 = res1
      damping = 1.0

      call nlpe_solver%factorize(dfdx)
      call nlpe_solver%solve(-f, dx)

      dx0 = maxval(abs(dx) / DX_LIM)
      if (dx0 > 1) dx = dx / dx0
      err = maxval(abs(dx))
      call dev%sys_nlpe%set_x(x0 + dx)
    end do

    deallocate(x0, f, dx)
  end subroutine

  subroutine gummel_quiet()
    !! Gummel iteration without printing (for batch tests)
    integer :: it, ci
    real :: err, err_pot, err_iref(2)
    real, allocatable :: pot0(:), iref0(:,:)

    real, parameter :: ATOL = 1e-6
    integer, parameter :: MAX_IT = 50

    allocate(pot0(dev%pot%data%n), iref0(dev%iref(1)%data%n, 2))

    do ci = CR_ELEC, CR_HOLE
      call ss_dd(ci)%init(dev%sys_dd(ci), log = .false., msg = "")
    end do

    it = 0
    err = huge(err)

    do while ((denorm(err, 'V') > ATOL) .and. (it < MAX_IT))
      it = it + 1

      pot0 = dev%pot%get()
      call solve_nlpe_quiet()
      err_pot = maxval(abs(dev%pot%get() - pot0))

      do ci = CR_ELEC, CR_HOLE
        iref0(:, ci) = dev%iref(ci)%get()
        call ss_dd(ci)%run()
        call ss_dd(ci)%select(1)
        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:, ci)))
      end do

      err = max(err_pot, err_iref(1), err_iref(2))
    end do

    deallocate(pot0, iref0)
  end subroutine

end program
