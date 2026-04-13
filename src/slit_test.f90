program slit_test
  !! Quick single-point EBIC solver for FIB slit geometry.
  !! Edit BEAM_X/BEAM_Y below, then: fargo run slit_test

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
  use taylor_remainder_m, only: taylor_test, taylor_test_per_block
  use grid_m,              only: IDX_VERTEX

  implicit none

  ! ═══════════════════════════════════════════════
  ! CONFIGURE: beam position (um)
  real, parameter :: BEAM_X = 13.75
  real, parameter :: BEAM_Y = 8.25
  ! ═══════════════════════════════════════════════

  type(device) :: dev
  type(polygon_src) :: input
  type(steady_state) :: ss, ss_dd(2)
  class(solver_real), allocatable :: nlpe_solver
  integer :: nx, ny, ict, ninput, ict_n, ict_p, ci
  real :: I_N, I_P, G_tot
  real, allocatable :: t_inp(:), V(:,:)

  ! Initialize
  call init_normconst(300.0)
  call dev%init("stem_ebic_slit.ini", 300.0)

  nx = size(dev%par%g1D(1)%x)
  ny = size(dev%par%g1D(2)%x)

  print "(A)",                       "═══ STEM-EBIC Slit Test ═══"
  print "(A,I0,A,I0,A,I0)", "Grid: ", nx, " x ", ny, " = ", nx*ny
  print "(A,F6.2,A,F6.2,A)", "Beam: (x, y) = (", BEAM_X, ", ", BEAM_Y, ") um"
  print "(A)", ""

  ! Set beam position (x = fixed coordinate, y = sweep coordinate)
  dev%par%reg_beam(1)%beam_x = norm(BEAM_X, 'um')
  dev%beam_pos%x = norm(BEAM_Y, 'um')

  ! Set up input: V=0 on all contacts, beam at BEAM_Y
  ninput = dev%par%nct + 1
  allocate(t_inp(2), V(ninput, 2))
  t_inp = [0.0, 1.0]
  do ict = 1, dev%par%nct
    V(ict, :) = 0.0
  end do
  V(ninput, :) = dev%beam_pos%x
  call input%init(t_inp, V)
  do ict = 1, dev%par%nct
    dev%volt(ict)%x = V(ict, 1)
  end do

  ! Initialize state before Taylor tests
  do ci = CR_ELEC, CR_HOLE
    call approx_imref(dev%par, dev%iref(ci), dev%volt)
  end do
  call approx_potential(dev%par, dev%pot, dev%iref)

  ! Evaluate beam generation for ALL carriers
  do ci = dev%par%ci0, dev%par%ci1
    call dev%calc_bgen(ci)%eval()
  end do
  G_tot = dev%calc_bgen(CR_ELEC)%G_tot
  print "(A,ES12.4,A)", "G_tot = ", G_tot, " pairs/(s·cm)"

  ! ── Taylor tests (at initialized state) ──

  ! Warm-up: single eval of sys_nlpe before DD Taylor test
  ! (sys_nlpe does NOT contain calc_srh — testing if it initializes shared state)
  call dev%sys_nlpe%eval()

  ! do ci = CR_ELEC, CR_HOLE
  !   print "(A)", "========================================"
  !   print "(A,I0,A)", " DD (carrier ", ci, ")"
  !   print "(A)", "========================================"

  !   ! Check if set_x leaks into the other carrier's density
  !   block
  !     real :: p_before, p_after
  !     real, allocatable :: x0(:), h(:)
  !     integer :: isz
  !     x0 = dev%sys_dd(ci)%get_x()
  !     allocate(h(size(x0)))
  !     call random_seed(size = isz)
  !     block
  !       integer, allocatable :: s(:)
  !       allocate(s(isz)); s = 99; call random_seed(put = s)
  !     end block
  !     call random_number(h); h = (h - 0.5) * 0.1
  !     if (ci == CR_ELEC) then
  !       p_before = dev%dens(CR_HOLE)%get([1,1])
  !     else
  !       p_before = dev%dens(CR_ELEC)%get([1,1])
  !     end if
  !     call dev%sys_dd(ci)%set_x(x0 + h)
  !     if (ci == CR_ELEC) then
  !       p_after = dev%dens(CR_HOLE)%get([1,1])
  !     else
  !       p_after = dev%dens(CR_ELEC)%get([1,1])
  !     end if
  !     print "(A,ES12.4,A,ES12.4,A,ES12.4)", "  PROVIDED LEAK CHECK: before=", p_before, &
  !       & " after=", p_after, " diff=", p_after - p_before
  !     call dev%sys_dd(ci)%set_x(x0)
  !     deallocate(x0, h)
  !   end block

  !   ! call taylor_test(dev%sys_dd(ci))
  !   ! call taylor_test_per_block(dev%sys_dd(ci))

  !   ! Diagnose SRH chain rule: compare V*dR/dn with chain rule result
  !   block
  !     integer :: n_dd, itest, istep, iloc
  !     integer, allocatable :: idx(:)
  !     real, allocatable :: x0_dd(:), fp_dd(:), fm_dd(:), Jcol_dd(:), ei(:)
  !     type(block_real), pointer :: df_dd
  !     real :: delta_fd, Ja_dd, Jn_dd, V_cell, dRdn_val, srh_expected
  !     real :: n_val, p_val, np_val, numer, denom, ni_sq, ni_val, tau_n, tau_p

  !     n_dd = dev%sys_dd(ci)%n
  !     delta_fd = 1e-7
  !     allocate(x0_dd(n_dd), fp_dd(n_dd), fm_dd(n_dd), Jcol_dd(n_dd), ei(n_dd))
  !     allocate(idx(dev%par%g%idx_dim))

  !     ! Get SRH parameters
  !     ni_sq = dev%calc_srh%ni_sq
  !     ni_val = dev%calc_srh%ni
  !     tau_n = dev%calc_srh%tau_n
  !     tau_p = dev%calc_srh%tau_p

  !     x0_dd = dev%sys_dd(ci)%get_x()
  !     call dev%sys_dd(ci)%eval(df = df_dd)

  !     print "(A)", ""
  !     print "(A)", "  SRH chain rule diagnostic:"
  !     print "(A)", "  DOF   J_chain   J_numeric  V*dR/dn    chain_err  srh_err    (ix,iy)"

  !     do istep = 0, 9
  !       iloc = 1 + istep * ((dev%par%transport(IDX_VERTEX, 0)%n - 1) / 9)
  !       itest = dev%sys_dd(ci)%i0(1) + iloc - 1
  !       idx = dev%par%transport(IDX_VERTEX, 0)%get_idx(iloc)

  !       ! Analytical from chain rule: (J * e_i)_i
  !       ei = 0.0; ei(itest) = 1.0
  !       Jcol_dd = 0.0
  !       call df_dd%mul_vec(ei, Jcol_dd)
  !       Ja_dd = Jcol_dd(itest)

  !       ! Numerical: centered FD
  !       call dev%sys_dd(ci)%set_x(x0_dd + delta_fd * ei)
  !       call dev%sys_dd(ci)%eval(f = fp_dd)
  !       call dev%sys_dd(ci)%set_x(x0_dd - delta_fd * ei)
  !       call dev%sys_dd(ci)%eval(f = fm_dd)
  !       Jn_dd = (fp_dd(itest) - fm_dd(itest)) / (2.0 * delta_fd)
  !       call dev%sys_dd(ci)%set_x(x0_dd)

  !       ! Directly compute V * dR/dn at this vertex
  !       V_cell = dev%par%tr_vol%get(idx)
  !       n_val = dev%dens(CR_ELEC)%get(idx)
  !       p_val = dev%dens(CR_HOLE)%get(idx)
  !       np_val = n_val * p_val
  !       numer = np_val - ni_sq
  !       denom = tau_p * (n_val + ni_val) + tau_n * (p_val + ni_val)
  !       if (ci == CR_ELEC) then
  !         dRdn_val = (p_val * denom - numer * tau_p) / (denom * denom)
  !       else
  !         dRdn_val = (n_val * denom - numer * tau_n) / (denom * denom)
  !       end if
  !       srh_expected = V_cell * dRdn_val

  !       print "(A,I6,5ES11.3,A,I4,A,I4,A)", "  ", itest, Ja_dd, Jn_dd, &
  !         & srh_expected, abs(Ja_dd - Jn_dd), abs(Jn_dd - (Ja_dd - srh_expected)) , &
  !         & " (", idx(1), ",", idx(2), ")"
  !     end do

  !     deallocate(x0_dd, fp_dd, fm_dd, Jcol_dd, ei, idx)
  !   end block
  ! end do

  ! (state already initialized above)

  ! 2. NLPE
  print "(A)", "NLPE..."
  call solve_nlpe()

  ! 3. Gummel — SKIPPED (unstable with Schottky spanning junction)
  ! call gummel()


  ! 4. Full Newton + save output
  print "(A)", "Newton..."
  call ss%init(dev%sys_full)
  ss%log    = .true.
  ss%msg    = "Newton: "
  ss%max_it = 200
  ss%dx_lim = norm(0.1, 'V')
  ss%atol   = norm(1e-9, 'V')
  call ss%init_output([ &
    & string("pot"), string("ndens"), string("pdens"), &
    ! & string("Ex"), string("Ey"), string("bgen_n"), &
    ! & string("V_P_CONTACT"), string("V_N_CONTACT"), &
    ! & string("I_N_CONTACT"), string("I_P_CONTACT")], "device.fbs")
     & string("Ex"), string("Ey"), string("bgen_n")],  "device.fbs")
  call ss%run(input = input, t_input = [0.0])

  ! ── Results ──

  print "(A)", ""
  print "(A)",          "═══ Results ═══"
  print "(A,F6.2,A,F6.2,A)", "Beam (x,y)  = (", BEAM_X, ", ", BEAM_Y, ") um"
  do ict = 1, dev%par%nct
    print "(A,A,A,ES12.4,A)", "I_", trim(dev%par%contacts(ict)%name), " = ", &
      & denorm(dev%curr(ict)%x, 'A') * 1e9, " nA"
  end do
  print "(A,ES12.4,A)", "I_expected  = ", 1.602e-19 * G_tot * 1e-4 * 1e9, " nA (ideal)"

  deallocate(t_inp, V)

contains

  subroutine solve_nlpe()
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

      print "(A,I4,A,ES10.3,A,ES10.3)", "  NLPE it ", it, &
        & ": err = ", denorm(err, 'V'), ", res = ", res1
    end do

    deallocate(x0, f, dx)
  end subroutine

  subroutine gummel()
    integer :: it, ci
    real :: err, err_pot, err_iref(2)
    real, allocatable :: pot0(:), iref0(:,:)

    real, parameter :: ATOL = 1e-6
    integer, parameter :: MAX_IT = 50

    allocate(pot0(dev%pot%data%n), iref0(dev%iref(1)%data%n, 2))

    do ci = CR_ELEC, CR_HOLE
      call ss_dd(ci)%init(dev%sys_dd(ci))
      ss_dd(ci)%log    = .true.
      ss_dd(ci)%msg    = "    DD: "
      ss_dd(ci)%max_it = 200
      ss_dd(ci)%dx_lim = 1e20
      ss_dd(ci)%error_if_not_converged = .false.
    end do

    it = 0
    err = huge(err)

    do while ((denorm(err, 'V') > ATOL) .and. (it < MAX_IT))
      it = it + 1

      pot0 = dev%pot%get()
      do ci = CR_ELEC, CR_HOLE
        iref0(:, ci) = dev%iref(ci)%get()
      end do

      call solve_nlpe()
      err_pot = maxval(abs(dev%pot%get() - pot0))

      do ci = CR_ELEC, CR_HOLE
        call ss_dd(ci)%run()
        call ss_dd(ci)%select(1)
        call dev%calc_iref(ci)%eval()
        err_iref(ci) = maxval(abs(dev%iref(ci)%get() - iref0(:, ci)))
      end do

      err = max(err_pot, err_iref(1), err_iref(2))
      print "(A,I3,A,ES10.3)", "  Gummel it ", it, ": err = ", denorm(err, 'V')
    end do

    deallocate(pot0, iref0)
  end subroutine


end program
