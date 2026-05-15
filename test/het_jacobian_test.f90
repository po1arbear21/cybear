program het_jacobian_test
  !! Validate the analytical sys_full Jacobian on the Si/SiGe heterojunction
  !! at converged V=0 equilibrium.
  !!
  !! The previous-session conclusion -- that the ~1e-21 A terminal current
  !! floor is structural FP cancellation in the (pot, dens) SG kernel -- rests
  !! on the assumption that Newton has actually converged to F(x) = 0 within
  !! machine epsilon. If the analytical Jacobian disagrees with the true
  !! Jacobian, Newton stalls early at a residual floor set by that
  !! inconsistency, which would look identical to an FP floor in the terminal
  !! current. This test independently checks which it is.
  !!
  !! Method: at the converged equilibrium state x*, run two threshold-free
  !! reference checks against the analytical Jacobian:
  !!
  !!   1. validate_taylor          ||F(x + h v) - F(x) - h J v|| in a log-log
  !!                               sweep. Slope ~ 2 iff J is correct, slope ~ 1
  !!                               otherwise. Threshold-free.
  !!   2. validate_jvp_richardson  Jv vs (4 g(h/2) - g(h))/3, accurate O(h^4).
  !!
  !! Two probe directions are exercised:
  !!   v_random  -- normalized random vector (catches widespread Jacobian bugs)
  !!   v_interf  -- nonzero only on the heterojunction-interface pot DoFs
  !!                (targets the band-edge term that this branch added)
  !!
  !! Build target: het_jacobian_test (see fargo.toml).

  use device_m,             only: device
  use normalization_m,      only: init_normconst, denorm
  use semiconductor_m,      only: CR_ELEC
  use grid_m,               only: IDX_VERTEX
  use input_src_m,          only: polygon_src
  use steady_state_m,       only: steady_state
  use string_m,             only: string
  use approx_m,             only: approx_imref, approx_potential
  use solver_base_m,        only: solver_real
  use solver_m,             only: default_solver_params, init_solver_real
  use input_m,              only: input_section
  use block_m,              only: block_real
  use jacobian_validator_m, only: validate_taylor, validate_jvp_richardson, log_log_slope
  use esystem_problem_m,    only: esystem_problem
  use dense_m,              only: dense_real

  implicit none

  type(device),       target      :: dev
  type(polygon_src)               :: input
  type(steady_state)              :: ss, ss_dd
  type(esystem_problem)           :: prob
  class(solver_real), allocatable :: nlpe_solver
  integer                         :: ict, ci, n, n_pass, n_fail, k_seed
  integer, allocatable            :: seed(:)
  real,    allocatable            :: t_inp(:), V(:,:), x_eq(:), v_rand(:), &
                                     v_pot_int(:), v_dens_int(:), v_bulk(:), f0(:)
  real                            :: I_L, I_R, res_norm
  integer, parameter              :: RAND_SEEDS(3) = [31337, 271828, 161803]
  character(len=32)               :: case_label

  ! h-sweep window. The slope-2 region is bounded below by the FP noise floor
  ! of ||F(x*)||, so hmin is chosen large enough that h*||Jv|| >> ||F(x*)||.
  real, parameter :: TAYLOR_HMIN  = 1e-5
  real, parameter :: TAYLOR_HMAX  = 1e-2
  integer, parameter :: TAYLOR_NH = 10
  real, parameter :: SLOPE_PASS   = 1.8   ! slope >= 1.8 counts as quadratic
  real, parameter :: JVP_TOL      = 1e-6  ! Richardson floor in normalized DoFs

  n_pass = 0
  n_fail = 0

  call init_normconst(300.0)
  call dev%init("si_sige.ini", 300.0)

  print "(A)", ""
  print "(A)", "================================================================"
  print "(A)", " het_jacobian_test  --  sys_full Jacobian validation at V=0 eq."
  print "(A)", "================================================================"

  ! ----- bring device to equilibrium (same path as het_test) -----
  allocate(t_inp(2), V(dev%par%nct, 2))
  t_inp = [0.0, 1.0]
  V     = 0.0
  call input%init(t_inp, V)
  do ict = 1, dev%par%nct
    dev%volt(ict)%x = 0.0
  end do

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
  call ss%init_output([string("pot"), string("ndens"), string("V_L"), string("V_R"), &
    &                  string("I_L"), string("I_R")], "het_jacobian_test.fbs")
  call ss%run(input = input, t_input = [0.0])

  I_L = denorm(dev%curr(1)%x, 'A')
  I_R = denorm(dev%curr(2)%x, 'A')
  print "(A)", ""
  print "(A,ES12.4,A)", "  terminal current I_L = ", I_L, " A"
  print "(A,ES12.4,A)", "  terminal current I_R = ", I_R, " A"

  ! ----- extract converged state and residual -----
  n = dev%sys_full%n
  allocate(x_eq(n), f0(n), v_rand(n), v_pot_int(n), v_dens_int(n), v_bulk(n))
  x_eq = dev%sys_full%get_x()
  call dev%sys_full%eval(f = f0)
  res_norm = norm2(f0)
  print "(A,I0)",       "  sys_full dof count n = ", n
  print "(A,ES12.4)",   "  ||F(x*)||_2          = ", res_norm
  print "(A)",          "  (a Jacobian bug would manifest as Newton stalling"
  print "(A)",          "   above the FP floor; slope-1 in the Taylor sweep confirms.)"

  prob%sys => dev%sys_full

  ! ----- structured probes: interface (pot, dens) and bulk-pot -----
  call build_one_var_dir(v_pot_int,  pot_active = .true.,  at_interface = .true. )
  call build_one_var_dir(v_dens_int, pot_active = .false., at_interface = .true. )
  call build_one_var_dir(v_bulk,     pot_active = .true.,  at_interface = .false.)

  if (norm2(v_pot_int) > 0.0) v_pot_int  = v_pot_int  / norm2(v_pot_int )
  if (norm2(v_dens_int)> 0.0) v_dens_int = v_dens_int / norm2(v_dens_int)
  if (norm2(v_bulk)    > 0.0) v_bulk     = v_bulk     / norm2(v_bulk    )

  call run_probe("v_pot_int  (interface, pot)", v_pot_int)
  call run_probe("v_dens_int (interface, ndens)", v_dens_int)
  call run_probe("v_bulk     (Si bulk, pot)", v_bulk)

  ! ----- random directions across three different seeds -----
  do k_seed = 1, size(RAND_SEEDS)
    call init_seed(seed_value = RAND_SEEDS(k_seed))
    call random_number(v_rand)
    v_rand = v_rand - 0.5
    v_rand = v_rand / norm2(v_rand)
    write (case_label, "(A,I0)") "v_rand  s=", RAND_SEEDS(k_seed)
    call run_probe(trim(case_label), v_rand)
    deallocate(seed)
  end do

  ! restore the converged state in case anything reads it after this program
  call dev%sys_full%set_x(x_eq)

  print "(A)", ""
  print "(A)", "================================================================"
  print "(A,I0,A,I0)", " Summary: ", n_pass, " passed, ", n_fail, " failed"
  print "(A)", "================================================================"

  if (n_fail > 0) error stop 1

contains

  subroutine solve_nlpe()
    integer                       :: it
    real                          :: err, res0, res1, damping, dx0
    real,             allocatable :: x0(:), f(:), dx(:)
    type(block_real), pointer     :: dfdx
    type(input_section)           :: solver_params
    real,    parameter :: ATOL   = 1e-10
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
    real,    parameter :: ATOL   = 1e-6
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

  subroutine init_seed(seed_value)
    integer, intent(in) :: seed_value
    integer :: m
    call random_seed(size = m)
    allocate(seed(m))
    seed = seed_value
    call random_seed(put = seed)
  end subroutine

  subroutine build_one_var_dir(vd, pot_active, at_interface)
    !! Construct a sys_full-sized perturbation that is nonzero on a single
    !! variable family (pot or ndens) at a specific subset of vertices.
    !! at_interface = true  -> vertices straddling a band_edge_v jump;
    !! at_interface = false -> a small patch of bulk vertices far from any
    !!                         jump (here: midpoints of each transport region).
    !! Returns the zero vector if no matching vertex set was found.
    real,    intent(out) :: vd(:)
    logical, intent(in)  :: pot_active
    logical, intent(in)  :: at_interface

    integer :: nx, ix
    real    :: be_left, be_right
    real, parameter :: EPS_JUMP = 1e-10
    logical, allocatable :: vert_mask(:)
    real, allocatable    :: x_save(:)

    vd = 0.0
    nx = size(dev%par%g1D(1)%x)
    allocate(vert_mask(nx), source = .false.)

    if (at_interface) then
      do ix = 1, nx - 1
        be_left  = dev%par%band_edge_v(CR_ELEC)%get([ix])
        be_right = dev%par%band_edge_v(CR_ELEC)%get([ix + 1])
        if (abs(be_right - be_left) > EPS_JUMP) then
          vert_mask(ix)     = .true.
          vert_mask(ix + 1) = .true.
        end if
      end do
    else
      ! pick a single bulk vertex near 25% of the device length; the test
      ! device's left Si region is contiguous so this is interior-to-Si
      ix = max(2, nx / 4)
      vert_mask(ix) = .true.
    end if

    ! Inject the perturbation through cybear's grid_data accessors so that the
    ! sys_full flat layout is queried, not assumed. This is robust to changes
    ! in the (pot, ndens, currents, V_*) ordering across configurations.
    allocate(x_save(dev%sys_full%n))
    x_save = dev%sys_full%get_x()
    do ix = 1, nx
      if (.not. vert_mask(ix)) cycle
      if (pot_active) then
        call dev%pot%set([ix], dev%pot%get([ix]) + 1.0)
      else
        call dev%dens(CR_ELEC)%set([ix], dev%dens(CR_ELEC)%get([ix]) + 1.0)
      end if
    end do
    vd = dev%sys_full%get_x() - x_save
    call dev%sys_full%set_x(x_save)
  end subroutine

  subroutine run_probe(label, v)
    character(*), intent(in) :: label
    real,         intent(in) :: v(:)

    real    :: slope_loc, slope_base_loc, jvp_err_loc, r1_max
    logical :: ok_taylor, ok_floor, ok_jvp

    if (norm2(v) == 0.0) then
      print "(A)", ""
      print "(A,A,A)", "  Skipping probe (", trim(label), ") — direction is zero"
      return
    end if

    print "(A)", ""
    print "(A,A)", "  Probe: ", trim(label)
    print "(A)",   "  ---------------------------------------------------------------"
    call dev%sys_full%set_x(x_eq)
    call validate_taylor(prob, x_eq, v, ok_taylor, slope_loc, slope_base_loc, &
      &                  hmin = TAYLOR_HMIN, hmax = TAYLOR_HMAX, nh = TAYLOR_NH, &
      &                  slope_tol = SLOPE_PASS)
    call dev%sys_full%set_x(x_eq)
    call print_taylor_table(label, v, r1_max)

    ! When R1 stays at the FP noise floor (||F(x*)||~1e-14) throughout the
    ! sweep, hJv matches F(x+hv)-F(x) to machine precision — the strongest
    ! possible "Jacobian correct" outcome. The fitted slope is meaningless
    ! noise in that regime, so accept "max R1 below floor" as a separate PASS.
    ok_floor = r1_max < 1e-10
    call report_taylor(label, ok_taylor .or. ok_floor, slope_loc, slope_base_loc, &
      &                floor_pass = ok_floor)

    call dev%sys_full%set_x(x_eq)
    call validate_jvp_richardson(prob, x_eq, v, ok_jvp, jvp_err_loc, tol = JVP_TOL)
    call report_jvp(label, ok_jvp, jvp_err_loc)
  end subroutine

  subroutine print_taylor_table(label, v, r1_max_out)
    !! Print R0(h) = ||F(x+hv) - F(x)|| and R1(h) = ||...-hJv|| at each h.
    !! Slope-by-eye is more robust than a fitted slope when the sweep window
    !! brushes against the FP noise floor of ||F(x*)||.
    character(*), intent(in)  :: label
    real,         intent(in)  :: v(:)
    real,         intent(out) :: r1_max_out

    type(dense_real)  :: D
    real, allocatable :: J(:,:), Jv(:), fh(:)
    real              :: h, R0_loc, R1_loc, sR0, sR1
    integer           :: k, ndof

    ndof = size(v)
    allocate(J(ndof, ndof), Jv(ndof), fh(ndof))

    call dev%sys_full%set_x(x_eq)
    call dev%sys_full%eval()
    call dev%sys_full%get_df(D)
    J = D%d
    Jv = matmul(J, v)

    print "(A,A,A)",         "    Taylor R(h) sweep (", trim(label), "):"
    print "(A)",             "      h            R0(h)=||F(x+hv)-F(x)||   R1(h)=||...-hJv||"

    r1_max_out = 0.0
    do k = 1, TAYLOR_NH
      h = TAYLOR_HMAX * (TAYLOR_HMIN / TAYLOR_HMAX) ** &
        & (real(k - 1) / real(TAYLOR_NH - 1))
      call dev%sys_full%set_x(x_eq + h * v)
      call dev%sys_full%eval(f = fh)
      R0_loc = norm2(fh - f0)
      R1_loc = norm2(fh - f0 - h * Jv)
      r1_max_out = max(r1_max_out, R1_loc)
      print "(A,ES10.3,A,ES15.6,A,ES15.6)", "      ", h, "  ", R0_loc, "  ", R1_loc
    end do

    call dev%sys_full%set_x(x_eq)

    ! Report local slopes on the two halves of the sweep (small-h vs large-h):
    ! a Jacobian bug shows as slope ~ 1 in the small-h half; an FP-noise floor
    ! shows as slope ~ 0 in the small-h half. Slope ~ 2 throughout means the
    ! Jacobian is correct in this direction.
    block
      real :: hs(TAYLOR_NH), R1s(TAYLOR_NH), R0s(TAYLOR_NH)
      integer :: half
      do k = 1, TAYLOR_NH
        hs(k) = TAYLOR_HMAX * (TAYLOR_HMIN / TAYLOR_HMAX) ** &
          &     (real(k - 1) / real(TAYLOR_NH - 1))
        call dev%sys_full%set_x(x_eq + hs(k) * v)
        call dev%sys_full%eval(f = fh)
        R0s(k) = norm2(fh - f0)
        R1s(k) = norm2(fh - f0 - hs(k) * Jv)
      end do
      half = TAYLOR_NH / 2
      sR0 = log_log_slope(hs(half:),     R0s(half:))
      sR1 = log_log_slope(hs(half:),     R1s(half:))
      print "(A,F5.2,A,F5.2,A)", &
        "      small-h slope (R0, R1) = (", sR0, ", ", sR1, ")"
      sR0 = log_log_slope(hs(:half),     R0s(:half))
      sR1 = log_log_slope(hs(:half),     R1s(:half))
      print "(A,F5.2,A,F5.2,A)", &
        "      large-h slope (R0, R1) = (", sR0, ", ", sR1, ")"
    end block

    call dev%sys_full%set_x(x_eq)
  end subroutine

  subroutine report_taylor(label, ok, slope, base, floor_pass)
    character(*), intent(in) :: label
    logical,      intent(in) :: ok
    real,         intent(in) :: slope, base
    logical,      intent(in), optional :: floor_pass
    character(:), allocatable :: verdict, note
    if (ok) then
      n_pass = n_pass + 1
      verdict = "PASS"
    else
      n_fail = n_fail + 1
      verdict = "FAIL"
    end if
    note = ""
    if (present(floor_pass)) then
      if (floor_pass) note = " (R1 at FP floor — Jacobian agrees to machine precision)"
    end if
    print "(A,A,A,F5.2,A,F5.2,A,A,A,A)", &
      "    taylor   (", label, ") slope=", slope, "  baseline=", base, &
      "   [", verdict, "]", note
  end subroutine

  subroutine report_jvp(label, ok, err)
    character(*), intent(in) :: label
    logical,      intent(in) :: ok
    real,         intent(in) :: err
    character(:), allocatable :: verdict
    if (ok) then
      n_pass = n_pass + 1
      verdict = "PASS"
    else
      n_fail = n_fail + 1
      verdict = "FAIL"
    end if
    print "(A,A,A,ES12.5,A,A,A)", &
      "    jvp+rich (", label, ") err=", err, "             [", verdict, "]"
  end subroutine

end program het_jacobian_test
