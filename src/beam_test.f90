program beam_test
  !! EBIC Beam Generation Test Runner

  use device_m,        only: device
  use normalization_m, only: init_normconst, denorm, norm
  use semiconductor_m, only: CR_ELEC, CR_HOLE
  use input_src_m,     only: polygon_src
  use steady_state_m,  only: steady_state
  use string_m,        only: string
  use approx_m,        only: approx_imref, approx_potential
  use taylor_remainder_m, only: taylor_test, taylor_test_per_block
  use block_m,         only: block_real
  use solver_base_m,   only: solver_real
  use solver_m,        only: default_solver_params, init_solver_real
  use input_m,         only: input_section

  implicit none

  type(device) :: dev
  type(polygon_src) :: input
  type(steady_state) :: ss, ss_dd(2)
  class(solver_real), allocatable :: nlpe_solver
  integer :: ix, iy, nx, ny, ict, ninput, ict_n, ci
  real :: G, beam_y, G_sum, G_tot, A, rel_err, I_N
  real, allocatable :: t_inp(:), V(:,:)

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

  ! Test 1b: Surface vertices (charged surface)
  if (dev%par%smc%charged_surf) then
    print "(A)", ""
    print "(A)", "Test 1b: Surface vertices for charged surface"
    print "(A)", ""
    print "(A,I0)", "  Number of surface vertices: ", dev%calc_scharge%n_surf_vert
    print "(A)", ""

    ! Print surf_idx array
    print "(A)", "  surf_idx array:"
    print "(A)", "  ┌───────┬───────┬───────┐"
    print "(A)", "  │   j   │  ix   │  iy   │"
    print "(A)", "  ├───────┼───────┼───────┤"
    do iy = 1, dev%calc_scharge%n_surf_vert
      print "(A,I4,A,I4,A,I4,A)", "  │", iy, "   │", &
        dev%calc_scharge%surf_idx(1, iy), "   │", &
        dev%calc_scharge%surf_idx(2, iy), "   │"
    end do
    print "(A)", "  └───────┴───────┴───────┘"
    print "(A)", ""

    ! 2D grid map with surface vertices highlighted
    print "(A)", "  2D Grid Map (S = surface vertex):"
    print "(A)", ""

    ! Header row
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

    ! Grid rows with surface vertices marked
    ! S = surface (in surf_idx), C = contact, · = interior
    do ix = 1, nx
      write(*, "(F10.0,A)", advance='no') denorm(dev%par%g1D(1)%x(ix), 'nm'), " │"
      do iy = 1, ny
        ! Check if this is a contact vertex
        if (dev%par%ict%get([ix, iy]) /= 0) then
          write(*, "(A)", advance='no') "     C"
        ! Check if this should be a surface vertex (x-boundary, not contact)
        elseif (ix == 1 .or. ix == nx) then
          write(*, "(A)", advance='no') "     S"
        else
          write(*, "(A)", advance='no') "     ·"
        end if
      end do
      print *
    end do
    print "(A)", ""
    print "(A)", "  Legend: S = surface vertex, C = contact, · = interior"
  end if

  ! Test 2: Beam generation values
  print "(A)", ""
  print "(A)", "Test 2: bgen_n values"
  print "(A)", ""

  ! Set beam at y = 2500 nm (middle of device)
  beam_y = 2000.0
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

  ! Test 4: Taylor remainder (Jacobian verification)
  print "(A)", ""
  print "(A)", "Test 4: Taylor remainder (Jacobian verification)"
  print "(A,F8.1,A)", "  Beam position: y = ", beam_y, " nm"
  print "(A)", ""

  ! Initialize state (same as jacobian_test)
  do ict = 1, dev%par%nct
    dev%volt(ict)%x = 0.0
  end do
  do ci = dev%par%ci0, dev%par%ci1
    call approx_imref(dev%par, dev%iref(ci), dev%volt)
  end do
  call approx_potential(dev%par, dev%pot, dev%iref)
  do ci = dev%par%ci0, dev%par%ci1
    call dev%calc_bgen(ci)%eval()
  end do

  ! print "(A)", "========================================"
  ! print "(A)", " NLPE"
  ! print "(A)", "========================================"
  ! call taylor_test(dev%sys_nlpe)
  ! call taylor_test_per_block(dev%sys_nlpe)

  ! print "(A)", "========================================"
  ! print "(A)", " Full Newton"
  ! print "(A)", "========================================"
  ! call taylor_test(dev%sys_full)
  ! call taylor_test_per_block(dev%sys_full)

  ! do ci = CR_ELEC, CR_HOLE
  !   print "(A)", "========================================"
  !   print "(A,I0,A)", " DD (carrier ", ci, ")"
  !   print "(A)", "========================================"
  !   call taylor_test(dev%sys_dd(ci))
  !   call taylor_test_per_block(dev%sys_dd(ci))
  ! end do

  ! Test 5: Terminal current from DD solver
  print "(A)", ""
  print "(A)", "Test 5: Terminal current (DD solve)"
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

  ! Debug: initial state at contacts
  print "(A)", ""
  print "(A)", "  === Initial state after approximation ==="
  print "(A,ES12.4,A)", "  Nc (edos_n) = ", denorm(dev%par%smc%edos(CR_ELEC), '1/cm^3'), " 1/cm^3"
  print "(A,ES12.4,A)", "  Nv (edos_p) = ", denorm(dev%par%smc%edos(CR_HOLE), '1/cm^3'), " 1/cm^3"
  do ict = 1, dev%par%nct
    print "(A,I0,A,A)", "  Contact ", ict, ": ", dev%par%contacts(ict)%name
    print "(A,F10.6,A)", "    phims    = ", denorm(dev%par%contacts(ict)%phims, 'V'), " V"
    print "(A,F10.6,A)", "    phi_b    = ", denorm(dev%par%contacts(ict)%phi_b, 'eV'), " eV"
    print "(A,ES12.4,A)", "    n0b_elec = ", denorm(dev%par%smc%edos(CR_ELEC) * exp(-dev%par%contacts(ict)%phi_b), '1/cm^3'), " 1/cm^3"
    print "(A,ES12.4,A)", "    n0b_hole = ", denorm(dev%par%smc%edos(CR_HOLE) &
      & * exp(-(dev%par%smc%band_gap - dev%par%contacts(ict)%phi_b)), '1/cm^3'), " 1/cm^3"
  end do
  ! Density at contact vertices (corners: ix=1, ix=nx)
  print "(A)", "  Carrier density at P_CONTACT (iy=1):"
  print "(A,ES12.4,A)", "    n(1,1)   = ", denorm(dev%dens(CR_ELEC)%get([1,1]), '1/cm^3'), " 1/cm^3"
  print "(A,ES12.4,A)", "    p(1,1)   = ", denorm(dev%dens(CR_HOLE)%get([1,1]), '1/cm^3'), " 1/cm^3"
  print "(A)", "  Carrier density at N_CONTACT (iy=ny):"
  print "(A,ES12.4,A)", "    n(1,ny)  = ", denorm(dev%dens(CR_ELEC)%get([1,ny]), '1/cm^3'), " 1/cm^3"
  print "(A,ES12.4,A)", "    p(1,ny)  = ", denorm(dev%dens(CR_HOLE)%get([1,ny]), '1/cm^3'), " 1/cm^3"
  print "(A)", "  Carrier density at interior (ix=1, iy=ny/2):"
  print "(A,ES12.4,A)", "    n(1,mid) = ", denorm(dev%dens(CR_ELEC)%get([1,ny/2]), '1/cm^3'), " 1/cm^3"
  print "(A,ES12.4,A)", "    p(1,mid) = ", denorm(dev%dens(CR_HOLE)%get([1,ny/2]), '1/cm^3'), " 1/cm^3"
  print "(A)", ""

  ! Step 2: Solve NLPE (non-linear Poisson equation)
  print "(A)", "  Step 2: Solving NLPE..."
  call solve_nlpe()

  ! Debug: state after NLPE
  print "(A)", ""
  print "(A)", "  === State after NLPE ==="
  print "(A)", "  Carrier density at P_CONTACT (iy=1):"
  print "(A,ES12.4,A)", "    n(1,1)   = ", denorm(dev%dens(CR_ELEC)%get([1,1]), '1/cm^3'), " 1/cm^3"
  print "(A,ES12.4,A)", "    p(1,1)   = ", denorm(dev%dens(CR_HOLE)%get([1,1]), '1/cm^3'), " 1/cm^3"
  print "(A)", "  Carrier density at N_CONTACT (iy=ny):"
  print "(A,ES12.4,A)", "    n(1,ny)  = ", denorm(dev%dens(CR_ELEC)%get([1,ny]), '1/cm^3'), " 1/cm^3"
  print "(A,ES12.4,A)", "    p(1,ny)  = ", denorm(dev%dens(CR_HOLE)%get([1,ny]), '1/cm^3'), " 1/cm^3"
  print "(A)", "  Carrier density at interior (ix=1, iy=ny/2):"
  print "(A,ES12.4,A)", "    n(1,mid) = ", denorm(dev%dens(CR_ELEC)%get([1,ny/2]), '1/cm^3'), " 1/cm^3"
  print "(A,ES12.4,A)", "    p(1,mid) = ", denorm(dev%dens(CR_HOLE)%get([1,ny/2]), '1/cm^3'), " 1/cm^3"
  print "(A)", ""

  ! Step 3: Gummel iteration
  print "(A)", "  Step 3: Running Gummel iteration..."
  call gummel()

  ! Step 4: Full Newton solver
  print "(A)", "  Step 4: Running full Newton solver..."
  call ss%init(dev%sys_full)
  ss%log    = .true.
  ss%msg    = "Newton: "
  ss%max_it = 200
  ss%dx_lim = norm(5.0, 'V')   ! limit Newton step to 5V (~200 kT)
  call ss%init_output([string("pot"), string("ndens"), string("pdens"), string("Ex"), string("Ey"), &
    & string("bgen_n"), string("V_P_CONTACT"), string("V_N_CONTACT"), string("I_N_CONTACT"), string("I_P_CONTACT")], "device.fbs")
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
    ! Debug: track where the max error occurs
    integer :: dbg_i, dbg_ix, dbg_iy, dbg_nclamp, dbg_imax
    real :: dbg_diff, dbg_d, dbg_dmin, dens_floor
    real, allocatable :: dbg_delta(:), dens0(:)
    character(len=8) :: dbg_name(2) = [character(len=8) :: "electron", "hole    "]
    ! Parameters (hardcoded for test)
    real, parameter :: ATOL = 1e-6
    integer, parameter :: MAX_IT = 30

    allocate(pot0(dev%pot%data%n), iref0(dev%iref(1)%data%n, 2))

    ! Initialize DD solvers
    do ci = CR_ELEC, CR_HOLE
      call ss_dd(ci)%init(dev%sys_dd(ci))
      ss_dd(ci)%log    = .true.
      ss_dd(ci)%msg    = "    DD: "
      ss_dd(ci)%max_it = 50
      ss_dd(ci)%error_if_not_converged = .false.
    end do
    dens_floor = norm(1.0, '1/cm^3')

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
        ! Save density before DD solve for diagnostics
        dens0 = dev%dens(ci)%get()
        call ss_dd(ci)%run()
        call ss_dd(ci)%select(1)

        ! Diagnose: where did DD change density the most?
        dbg_delta = dev%dens(ci)%get() - dens0
        dbg_imax = maxloc(abs(dbg_delta), dim=1)
        dbg_ix = mod(dbg_imax - 1, nx) + 1
        dbg_iy = (dbg_imax - 1) / nx + 1
        print "(A,A,A)", "      [DD-DIAG] ", trim(dbg_name(ci)), " max density change:"
        print "(A,I0,A,I0,A,F7.1,A,F7.1,A)", &
          & "        at (", dbg_ix, ",", dbg_iy, ") = (", &
          & denorm(dev%par%g1D(1)%x(dbg_ix), 'nm'), ",", &
          & denorm(dev%par%g1D(2)%x(dbg_iy), 'nm'), " nm)"
        print "(A,ES12.4,A,ES12.4,A,ES12.4,A)", &
          & "        before: ", denorm(dens0(dbg_imax), '1/cm^3'), &
          & "  after: ", denorm(dens0(dbg_imax) + dbg_delta(dbg_imax), '1/cm^3'), &
          & "  delta: ", denorm(dbg_delta(dbg_imax), '1/cm^3'), " cm^-3"
        deallocate(dbg_delta, dens0)

        ! Debug: show raw density before clamping
        dbg_dmin = huge(dbg_dmin)
        dbg_nclamp = 0
        do dbg_iy = 1, ny
          do dbg_ix = 1, nx
            dbg_d = dev%dens(ci)%get([dbg_ix, dbg_iy])
            if (dbg_d < dbg_dmin) dbg_dmin = dbg_d
            if (dbg_d < dens_floor) then
              if (dbg_nclamp < 5) then
                print "(A,A,A,I0,A,I0,A,F7.1,A,F7.1,A,ES10.3,A)", &
                  & "      [CLAMP] ", trim(dbg_name(ci)), " dens(", dbg_ix, ",", dbg_iy, &
                  & ") = (", denorm(dev%par%g1D(1)%x(dbg_ix), 'nm'), ",", &
                  & denorm(dev%par%g1D(2)%x(dbg_iy), 'nm'), " nm): ", &
                  & denorm(dbg_d, '1/cm^3'), " cm^-3"
              end if
              dbg_nclamp = dbg_nclamp + 1
              call dev%dens(ci)%set([dbg_ix, dbg_iy], dens_floor)
            end if
          end do
        end do
        if (dbg_nclamp > 0) then
          print "(A,A,A,I0,A,I0,A,ES10.3,A)", &
            & "      [CLAMP] ", trim(dbg_name(ci)), ": ", dbg_nclamp, " of ", nx*ny, &
            & " points clamped (min was ", denorm(dbg_dmin, '1/cm^3'), " cm^-3)"
        end if

        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:, ci)))
      end do

      err = max(err_pot, err_iref(1), err_iref(2))
      print "(A,I3,A,ES10.3)", "    Gummel it ", it, ": err = ", denorm(err, 'V')

      ! Debug: identify which variable and grid node drives the error
      ! Potential
      dbg_delta = abs(dev%pot%get() - pot0)
      dbg_i = maxloc(dbg_delta, dim=1)
      dbg_ix = mod(dbg_i - 1, nx) + 1
      dbg_iy = (dbg_i - 1) / nx + 1
      print "(A,ES10.3,A,I0,A,I0,A,F7.1,A,F7.1,A)", &
        & "      pot:  d=", denorm(dbg_delta(dbg_i), 'V'), &
        & "  at (", dbg_ix, ",", dbg_iy, ") = (", &
        & denorm(dev%par%g1D(1)%x(dbg_ix), 'nm'), ",", &
        & denorm(dev%par%g1D(2)%x(dbg_iy), 'nm'), " nm)"
      deallocate(dbg_delta)
      ! Imref
      do ci = CR_ELEC, CR_HOLE
        dbg_delta = abs(dev%iref(ci)%get() - iref0(:, ci))
        dbg_i = maxloc(dbg_delta, dim=1)
        dbg_ix = mod(dbg_i - 1, nx) + 1
        dbg_iy = (dbg_i - 1) / nx + 1
        print "(A,A,A,ES10.3,A,I0,A,I0,A,F7.1,A,F7.1,A)", &
          & "      iref_", trim(dbg_name(ci)), ": d=", denorm(dbg_delta(dbg_i), 'V'), &
          & "  at (", dbg_ix, ",", dbg_iy, ") = (", &
          & denorm(dev%par%g1D(1)%x(dbg_ix), 'nm'), ",", &
          & denorm(dev%par%g1D(2)%x(dbg_iy), 'nm'), " nm)"
        deallocate(dbg_delta)
      end do
    end do

    deallocate(pot0, iref0)
  end subroutine

end program
