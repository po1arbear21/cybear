program schottky_test
  !! Schottky Contact Test Runner
  !!
  !! Test 1: Equilibrium (V=0) — solve and verify zero terminal currents

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
  integer :: ix, nx, ict, ninput, ci
  real :: phi_b_eV, n_contact, n_expected, rel_err
  real :: I_schottky, I_ohmic
  real, allocatable :: t_inp(:), V(:,:)

  ! ====================================================================
  ! Initialize
  ! ====================================================================
  call init_normconst(300.0)
  call dev%init("schottky_minimal.ini", 300.0)

  nx = size(dev%par%g1D(1)%x)

  ! ====================================================================
  ! Test 1: Grid overview
  ! ====================================================================
  print "(A)", ""
  print "(A)", "========================================"
  print "(A)", " Schottky Diode Test"
  print "(A)", "========================================"
  print "(A)", ""
  print "(A)", "Test 1: Grid and device structure"
  print "(A)", ""
  print "(A,I0,A)", "  Grid: ", nx, " vertices"
  print "(A)", ""

  ! X positions
  print "(A)", "  X grid:"
  print "(A)", "  ┌───────┬────────────┐"
  print "(A)", "  │  ix   │   x [nm]   │"
  print "(A)", "  ├───────┼────────────┤"
  do ix = 1, min(nx, 10)
    print "(A,I4,A,F10.1,A)", "  │", ix, "   │", denorm(dev%par%g1D(1)%x(ix), 'nm'), "  │"
  end do
  if (nx > 10) then
    print "(A)", "  │ ...   │    ...     │"
    print "(A,I4,A,F10.1,A)", "  │", nx, "   │", denorm(dev%par%g1D(1)%x(nx), 'nm'), "  │"
  end if
  print "(A)", "  └───────┴────────────┘"
  print "(A)", ""

  ! Contact info
  print "(A)", "  Contacts:"
  do ict = 1, dev%par%nct
    print "(A,I0,A,A,A,F8.4)", "    ", ict, ": ", dev%par%contacts(ict)%name, &
      "  phims = ", denorm(dev%par%contacts(ict)%phims, 'V')
  end do
  print "(A)", ""

  ! ====================================================================
  ! Test 2: Equilibrium solve (V = 0)
  ! ====================================================================
  print "(A)", "Test 2: Equilibrium solve (V = 0)"
  print "(A)", ""

  ! Set up input: all contacts at 0V
  ninput = dev%par%nct
  allocate(t_inp(2), V(ninput, 2))
  t_inp = [0.0, 1.0]
  V = 0.0
  call input%init(t_inp, V)

  ! Set voltages
  do ict = 1, dev%par%nct
    dev%volt(ict)%x = 0.0
  end do

  ! Step 1: Approximate initial conditions
  print "(A)", "  Step 1: Approximating initial conditions..."
  do ci = CR_ELEC, CR_HOLE
    call approx_imref(dev%par, dev%iref(ci), dev%volt)
  end do
  call approx_potential(dev%par, dev%pot, dev%iref)

  ! Step 2: Solve NLPE
  print "(A)", "  Step 2: Solving NLPE..."
  call solve_nlpe()

  ! Step 3: Gummel iteration
  print "(A)", "  Step 3: Running Gummel iteration..."
  call gummel()

  ! Step 4: Full Newton solver
  print "(A)", "  Step 4: Running full Newton solver..."
  call ss%init(dev%sys_full)
  ss%msg = "Newton: "
  call ss%init_output([string("pot"), string("ndens"), string("pdens"), &
    & string("V_SCHOTTKY"), string("V_OHMIC"), &
    & string("I_SCHOTTKY"), string("I_OHMIC")], "schottky_test.fbs")
  call ss%run(input = input, t_input = [0.0])

  ! ====================================================================
  ! Results
  ! ====================================================================
  print "(A)", ""
  print "(A)", "========================================"
  print "(A)", " Results"
  print "(A)", "========================================"
  print "(A)", ""

  ! Terminal currents
  I_schottky = denorm(dev%curr(1)%x, 'A')
  I_ohmic    = denorm(dev%curr(2)%x, 'A')
  print "(A,ES12.4,A)", "  I_SCHOTTKY = ", I_schottky, " A"
  print "(A,ES12.4,A)", "  I_OHMIC    = ", I_ohmic,    " A"
  print "(A)", ""

  ! Check: at equilibrium, currents should be ~0
  print "(A)", "  Equilibrium check (currents should be ~0):"
  if (abs(I_schottky) < 1e-15) then
    print "(A,ES10.2,A)", "    I_SCHOTTKY: PASS  (", I_schottky, " A)"
  else
    print "(A,ES10.2,A)", "    I_SCHOTTKY: FAIL  (", I_schottky, " A)"
  end if
  if (abs(I_ohmic) < 1e-15) then
    print "(A,ES10.2,A)", "    I_OHMIC:    PASS  (", I_ohmic, " A)"
  else
    print "(A,ES10.2,A)", "    I_OHMIC:    FAIL  (", I_ohmic, " A)"
  end if
  print "(A)", ""

  ! Built-in potential check
  print "(A,F10.6,A)", "  Built-in potential (pot at x=0):    ", &
    denorm(dev%pot%get([1]), 'V'), " V"
  print "(A,F10.6,A)", "  Built-in potential (pot at x=Lx):   ", &
    denorm(dev%pot%get([nx]), 'V'), " V"
  print "(A,F10.6,A)", "  V_bi = pot(Lx) - pot(0):            ", &
    denorm(dev%pot%get([nx]) - dev%pot%get([1]), 'V'), " V"
  print "(A)", ""

  ! Carrier density at Schottky contact
  print "(A)", "  Carrier density at Schottky contact (ix=1):"
  n_contact = denorm(dev%dens(CR_ELEC)%get([1]), '1/cm^3')
  phi_b_eV = denorm(dev%par%contacts(1)%phi_b, 'eV')
  n_expected = denorm(dev%par%smc%edos(CR_ELEC), '1/cm^3') * exp(-dev%par%contacts(1)%phi_b)
  print "(A,ES12.4,A)", "    n(0)      = ", n_contact,  " 1/cm^3"
  print "(A,ES12.4,A)", "    n_expected = ", n_expected, " 1/cm^3  (Nc * exp(-phi_b/kT))"
  rel_err = abs(n_contact - n_expected) / n_expected * 100.0
  print "(A,F8.2,A)",   "    Error     = ", rel_err, " %"
  if (rel_err < 5.0) then
    print "(A)", "    PASS"
  else
    print "(A)", "    FAIL"
  end if
  print "(A)", ""

  deallocate(t_inp, V)

contains

  subroutine solve_nlpe()
    !! Solve non-linear Poisson equation using Newton iteration
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

      print "(A,I4,A,ES10.3,A,ES10.3)", "    NLPE it ", it, ": err = ", denorm(err, 'V'), ", res = ", res1
    end do

    deallocate(x0, f, dx)
  end subroutine

  subroutine gummel()
    !! Gummel iteration: alternate between NLPE and DD solves
    integer :: it, ci
    real :: err, err_pot, err_iref(2)
    real, allocatable :: pot0(:), iref0(:,:)

    real, parameter :: ATOL = 1e-6
    integer, parameter :: MAX_IT = 50

    allocate(pot0(dev%pot%data%n), iref0(dev%iref(1)%data%n, 2))

    do ci = CR_ELEC, CR_HOLE
      call ss_dd(ci)%init(dev%sys_dd(ci))
    end do

    it = 0
    err = huge(err)

    do while ((denorm(err, 'V') > ATOL) .and. (it < MAX_IT))
      it = it + 1

      pot0 = dev%pot%get()
      call solve_nlpe()
      err_pot = maxval(abs(dev%pot%get() - pot0))

      do ci = CR_ELEC, CR_HOLE
        iref0(:, ci) = dev%iref(ci)%get()
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

end program
