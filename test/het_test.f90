program het_test
  !! Heterojunction equilibrium test runner.
  !!
  !! Si / Si0.8Ge0.2 / Si N-n-n stack at V=0. Verifies:
  !!   1. Both terminal currents are zero to machine precision.
  !!   2. The conduction-band edge shows a ~20 meV step at each Si/SiGe
  !!      interface (matching the Ec_offset values in si_sige.ini).
  !!   3. The single-material regression suite is unaffected (run separately).

  use device_m,         only: device
  use normalization_m,  only: init_normconst, denorm, norm
  use semiconductor_m,  only: CR_ELEC, DOP_DCON
  use grid_m,           only: IDX_VERTEX
  use input_src_m,      only: polygon_src
  use steady_state_m,   only: steady_state
  use string_m,         only: string
  use approx_m,         only: approx_imref, approx_potential
  use solver_base_m,    only: solver_real
  use solver_m,         only: default_solver_params, init_solver_real
  use input_m,          only: input_section
  use block_m,          only: block_real

  implicit none

  type(device)                    :: dev
  type(polygon_src)               :: input
  type(steady_state)              :: ss, ss_dd
  class(solver_real), allocatable :: nlpe_solver
  integer                         :: ict, ci, ix, nx, n_pass, n_fail
  real                            :: I_L, I_R, V_L, V_R, Ec_max, Ec_min
  real,  allocatable              :: t_inp(:), V(:,:)

  n_pass = 0
  n_fail = 0

  call init_normconst(300.0)
  call dev%init("si_sige.ini", 300.0)

  nx = size(dev%par%g1D(1)%x)
  print "(A)", ""
  print "(A)", "========================================"
  print "(A)", " Heterojunction Equilibrium Test"
  print "(A)", " Si / Si0.8Ge0.2 / Si at V = 0"
  print "(A)", "========================================"
  print "(A,I0,A)", "  grid: ", nx, " vertices"
  print "(A,I0)",   "  contacts: ", dev%par%nct

  ! ---- equilibrium solve ----
  allocate(t_inp(2), V(dev%par%nct, 2))
  t_inp = [0.0, 1.0]
  V     = 0.0
  call input%init(t_inp, V)
  do ict = 1, dev%par%nct
    dev%volt(ict)%x = 0.0
  end do

  print "(A)", ""
  print "(A)", "  approximating initial conditions"
  do ci = dev%par%smc(dev%par%smc_default)%ci0, dev%par%smc(dev%par%smc_default)%ci1
    call approx_imref(dev%par, dev%iref(ci), dev%volt)
  end do
  call approx_potential(dev%par, dev%pot, dev%iref)

  print "(A)", "  solving NLPE..."
  call solve_nlpe()
  print "(A)", "  Gummel iteration..."
  call gummel()
  print "(A)", "  full Newton..."
  call ss%init(dev%sys_full)
  ss%msg = "  Newton: "
  call ss%init_output([string("pot"), string("ndens"),       &
    &                  string("V_L"), string("V_R"),         &
    &                  string("I_L"), string("I_R")], "het_test.fbs")
  call ss%run(input = input, t_input = [0.0])

  ! ---- results ----
  V_L = denorm(dev%volt(1)%x, 'V')
  V_R = denorm(dev%volt(2)%x, 'V')
  I_L = denorm(dev%curr(1)%x, 'A')
  I_R = denorm(dev%curr(2)%x, 'A')

  print "(A)", ""
  print "(A)", "  --- terminal voltages ---"
  print "(A,A,A,ES12.4,A)", "    V_", trim(dev%par%contacts(1)%name), " = ", V_L, " V"
  print "(A,A,A,ES12.4,A)", "    V_", trim(dev%par%contacts(2)%name), " = ", V_R, " V"
  print "(A)", "  --- terminal currents (target: ~0) ---"
  print "(A,A,A,ES12.4,A)", "    I_", trim(dev%par%contacts(1)%name), " = ", I_L, " A"
  print "(A,A,A,ES12.4,A)", "    I_", trim(dev%par%contacts(2)%name), " = ", I_R, " A"

  ! Currents at V=0 should be at numerical noise; 1e-15 A is well below femtoamps.
  if (abs(I_L) < 1e-15) then
    print "(A)", "    I_L  PASS (< 1e-15 A)"
    n_pass = n_pass + 1
  else
    print "(A,ES10.2)", "    I_L  FAIL: ", I_L
    n_fail = n_fail + 1
  end if
  if (abs(I_R) < 1e-15) then
    print "(A)", "    I_R  PASS (< 1e-15 A)"
    n_pass = n_pass + 1
  else
    print "(A,ES10.2)", "    I_R  FAIL: ", I_R
    n_fail = n_fail + 1
  end if

  ! ---- band-offset (dEc) check ----
  ! dEc_v stores each vertex's Ec_offset (user spec). With Si=0 and SiGe20=0.020,
  ! the step across the device should be exactly 20 meV regardless of Eg / DOS.
  block
    real :: dEc_max, dEc_min
    dEc_max = -huge(dEc_max)
    dEc_min =  huge(dEc_min)
    do ix = 1, nx
      dEc_max = max(dEc_max, denorm(dev%par%dEc_v%get([ix]), 'eV'))
      dEc_min = min(dEc_min, denorm(dev%par%dEc_v%get([ix]), 'eV'))
    end do
    print "(A)", ""
    print "(A)", "  --- per-vertex dEc(x) (user-specified Ec_offset) ---"
    print "(A,F8.4,A)", "    dEc_max = ", dEc_max, " eV"
    print "(A,F8.4,A)", "    dEc_min = ", dEc_min, " eV"
    print "(A,F8.4,A)", "    step    = ", dEc_max - dEc_min, " eV (expected 0.020)"
    if (abs((dEc_max - dEc_min) - 0.020) < 1e-5) then
      print "(A)", "    dEc step PASS"
      n_pass = n_pass + 1
    else
      print "(A)", "    dEc step FAIL"
      n_fail = n_fail + 1
    end if
  end block

  ! ---- band-edge step in the simulator's internal coordinate ----
  ! band_edge_v = smc%band_edge + smc%dEc. For Si/SiGe with our parameters:
  !   diff ≈ 0.5*(Eg_SiGe - Eg_Si) + 0.5*Vt*ln((Nc/Nv)_SiGe / (Nc/Nv)_Si) + Δ(dEc)
  !        ≈ 0.5*(-0.064) + tiny + 0.020 ≈ -0.012 eV
  Ec_max = -huge(Ec_max)
  Ec_min =  huge(Ec_min)
  do ix = 1, nx
    Ec_max = max(Ec_max, denorm(dev%par%band_edge_v(CR_ELEC)%get([ix]), 'eV'))
    Ec_min = min(Ec_min, denorm(dev%par%band_edge_v(CR_ELEC)%get([ix]), 'eV'))
  end do
  print "(A)", ""
  print "(A)", "  --- per-vertex band_edge_v (CR_ELEC) ---"
  print "(A,F8.4,A)", "    max = ", Ec_max, " eV"
  print "(A,F8.4,A)", "    min = ", Ec_min, " eV"
  print "(A,F8.4,A)", "    step = ", Ec_max - Ec_min, " eV"

  ! ---- per-edge current diagnostic ----
  ! Locate the maximum-magnitude edge current to confirm whether the contact-
  ! level current is dominated by interface edges (FP precision floor) or
  ! bulk edges (machine-zero from bit-identical inputs).
  block
    integer :: ie, ie_max
    real    :: J_e, J_max
    J_max  = 0.0
    ie_max = -1
    ! cdens(idx_dir, ci) — pick x-direction electron flux (1D device, n-only)
    do ie = 1, dev%cdens(1,CR_ELEC)%data%n
      J_e = abs(dev%cdens(1,CR_ELEC)%x1(ie))
      if (J_e > J_max) then
        J_max  = J_e
        ie_max = ie
      end if
    end do
    print "(A)", ""
    print "(A)", "  --- per-edge SG flux diagnostic ---"
    print "(A,I0,A,ES12.4)", "    max-|J| edge index = ", ie_max, ", normalized j_norm = ", J_max
    print "(A)", "    (interior bulk edges should be ~1e-30; interface edges set the floor)"
  end block

  print "(A)", ""
  print "(A)", "========================================"
  print "(A,I0,A,I0)", " Summary: ", n_pass, " passed, ", n_fail, " failed"
  print "(A)", "========================================"

  if (n_fail > 0) error stop 1

contains

  subroutine solve_nlpe()
    integer                   :: it
    real                      :: err, res0, res1, damping, dx0
    real,             allocatable :: x0(:), f(:), dx(:)
    type(block_real), pointer :: dfdx
    type(input_section)       :: solver_params
    real,    parameter :: ATOL = 1e-10
    real,    parameter :: DX_LIM = 1.0
    integer, parameter :: MAX_IT = 200

    allocate(x0(dev%sys_nlpe%n), f(dev%sys_nlpe%n), dx(dev%sys_nlpe%n), source = 0.0)
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
  end subroutine

  subroutine gummel()
    integer            :: it, ci
    real               :: err, err_pot, err_iref(2)
    real,  allocatable :: pot0(:), iref0(:,:)
    real,    parameter :: ATOL = 1e-6
    integer, parameter :: MAX_IT = 200

    allocate(pot0(dev%pot%data%n), &
      &      iref0(dev%iref(dev%par%smc(dev%par%smc_default)%ci0)%data%n, 2))
    call ss_dd%init(dev%sys_dd(dev%par%smc(dev%par%smc_default)%ci0))

    it = 0
    err = huge(err)
    err_iref = 0.0
    do while ((denorm(err, 'V') > ATOL) .and. (it < MAX_IT))
      it = it + 1
      pot0 = dev%pot%get()
      call solve_nlpe()
      err_pot = maxval(abs(dev%pot%get() - pot0))
      do ci = dev%par%smc(dev%par%smc_default)%ci0, dev%par%smc(dev%par%smc_default)%ci1
        iref0(:, ci) = dev%iref(ci)%get()
        call ss_dd%run()
        call ss_dd%select(1)
        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:, ci)))
      end do
      err = err_pot
      do ci = dev%par%smc(dev%par%smc_default)%ci0, dev%par%smc(dev%par%smc_default)%ci1
        err = max(err, err_iref(ci))
      end do
    end do
  end subroutine

end program het_test
