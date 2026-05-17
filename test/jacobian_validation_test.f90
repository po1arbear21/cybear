module bratu_problem_m
  !! Test fixture for jacobian_validation_test.
  !!
  !! 1D Bratu residual on a uniform grid with Dirichlet BCs:
  !!
  !!   F_i(x) = x_{i-1} - 2 x_i + x_{i+1} + lam * exp(x_i)   (interior)
  !!   F_1    = x_1,    F_n = x_n                            (Dirichlet rows)
  !!
  !! The exact analytical Jacobian is tridiagonal; eval_jac returns it after
  !! optionally injecting one of several controlled defects so the test can
  !! exercise each validator against known-wrong Jacobians.
  !!
  !! The "exp_term" flag makes the lam*exp(x_i) contribution toggleable; this
  !! plays the role of a physics term in cybear DD and lets the test suite
  !! exercise the paired-flag harness and toggle-diff workflow on a fixture
  !! whose analytical answer is known by hand.

  use jacobian_validator_m, only: jacobian_problem_cmplx

  implicit none

  private
  public :: bratu_problem
  public :: BREAK_NONE, BREAK_SCALED, BREAK_DROPPED_DIAG, BREAK_TRANSPOSED, BREAK_NOISE
  public :: BREAK_TERM_JAC_STALE

  integer, parameter :: BREAK_NONE           = 0
  integer, parameter :: BREAK_SCALED         = 1   ! J' = 1.5 J
  integer, parameter :: BREAK_DROPPED_DIAG   = 2   ! zero one diagonal entry
  integer, parameter :: BREAK_TRANSPOSED     = 3   ! swap one off-diagonal pair, asymmetrically
  integer, parameter :: BREAK_NOISE          = 4   ! small additive noise
  integer, parameter :: BREAK_TERM_JAC_STALE = 5   ! exp_term=off in F but on in J  (paired-flag bug)

  type, extends(jacobian_problem_cmplx) :: bratu_problem
    integer :: n        = 0
    real    :: lam      = 1.0
    integer :: break    = BREAK_NONE
    logical :: exp_term = .true.       ! toggles the nonlinear exp contribution
  contains
    procedure :: eval_real  => bratu_eval_real
    procedure :: eval_jac   => bratu_eval_jac
    procedure :: eval_cmplx => bratu_eval_cmplx
  end type

contains

  subroutine bratu_eval_real(this, x, f)
    class(bratu_problem), intent(in)  :: this
    real,                 intent(in)  :: x(:)
    real,                 intent(out) :: f(:)
    integer :: i
    f(1)      = x(1)
    f(this%n) = x(this%n)
    do i = 2, this%n - 1
      f(i) = x(i-1) - 2.0 * x(i) + x(i+1)
      if (this%exp_term) f(i) = f(i) + this%lam * exp(x(i))
    end do
  end subroutine

  subroutine bratu_eval_cmplx(this, z, fz)
    class(bratu_problem), intent(in)  :: this
    complex,              intent(in)  :: z(:)
    complex,              intent(out) :: fz(:)
    integer :: i
    fz(1)      = z(1)
    fz(this%n) = z(this%n)
    do i = 2, this%n - 1
      fz(i) = z(i-1) - 2.0 * z(i) + z(i+1)
      if (this%exp_term) fz(i) = fz(i) + this%lam * exp(z(i))
    end do
  end subroutine

  subroutine bratu_eval_jac(this, x, J)
    !! Return the (possibly intentionally broken) claimed Jacobian.
    class(bratu_problem), intent(in)  :: this
    real,                 intent(in)  :: x(:)
    real,                 intent(out) :: J(:,:)
    integer :: i, k, kt
    real    :: tmp
    real, allocatable :: noise(:,:)
    integer, allocatable :: seed(:)
    integer :: m
    logical :: include_exp_in_jac

    ! BREAK_TERM_JAC_STALE simulates a paired-flag bug: F drops the exp term
    ! but the Jacobian assembly forgot to drop the matching contribution.
    ! Any well-built physics toggle in cybear must avoid this failure mode.
    if (this%break == BREAK_TERM_JAC_STALE) then
      include_exp_in_jac = .true.
    else
      include_exp_in_jac = this%exp_term
    end if

    J = 0.0
    J(1, 1)           = 1.0
    J(this%n, this%n) = 1.0
    do i = 2, this%n - 1
      J(i, i-1) = 1.0
      J(i, i  ) = -2.0
      if (include_exp_in_jac) J(i, i) = J(i, i) + this%lam * exp(x(i))
      J(i, i+1) = 1.0
    end do

    select case (this%break)
    case (BREAK_NONE)
      ! nothing to inject
    case (BREAK_SCALED)
      J = 1.5 * J
    case (BREAK_DROPPED_DIAG)
      k = this%n / 2
      J(k, k) = 0.0
    case (BREAK_TRANSPOSED)
      ! Swap J(k,k+1) <-> J(k+1,k) and add 0.7 so the swap is asymmetric
      ! and cannot self-cancel against the symmetric structure.
      k         = this%n / 3
      kt        = k + 1
      tmp       = J(k, kt)
      J(k, kt)  = J(kt, k) + 0.7
      J(kt, k)  = tmp + 0.7
    case (BREAK_NOISE)
      allocate (noise(size(J,1), size(J,2)))
      call random_seed(size = m)
      allocate (seed(m))
      seed = 7
      call random_seed(put = seed)
      call random_number(noise)
      J = J + 1.0e-9 * (noise - 0.5)
      deallocate (noise, seed)
    case (BREAK_TERM_JAC_STALE)
      ! Nothing to do here -- the bug is encoded by include_exp_in_jac above.
    end select
  end subroutine

end module


program jacobian_validation_test
  !! Proof-test for the four validators in jacobian_validator_m.
  !!
  !! For each of {correct, scaled, dropped_diag, transposed, noise} the
  !! claimed Jacobian is altered and all four validators are run. A test
  !! "passes" iff each method's verdict (ok / not ok) matches the expectation
  !! (the correct J must be accepted; broken J's must be rejected -- except
  !! that noise at 1e-9 is below the FD cancellation floor and is expected to
  !! slip past FD, demonstrating CSD's stronger discrimination).

  use jacobian_validator_m, only: validate_finite_diff,   &
    &                              validate_complex_step, &
    &                              validate_taylor,       &
    &                              validate_jvp_richardson, &
    &                              print_top_defects

  use jacval_snapshot_m,    only: jacval_snapshot_dump,       &
    &                              index_map_t, physics_flags_t

  use bratu_problem_m,      only: bratu_problem,                              &
    &                              BREAK_NONE, BREAK_SCALED, BREAK_DROPPED_DIAG, &
    &                              BREAK_TRANSPOSED, BREAK_NOISE,               &
    &                              BREAK_TERM_JAC_STALE

  implicit none

  type(bratu_problem) :: prob
  real, allocatable   :: x(:), v(:)
  integer             :: n, i, m, n_pass, n_fail
  integer, allocatable :: seed(:)

  n = 12
  allocate (x(n), v(n))

  do i = 1, n
    x(i) = 0.3 * sin(3.0 * real(i) / real(n))
  end do

  call random_seed(size = m)
  allocate (seed(m))
  seed = 42
  call random_seed(put = seed)
  call random_number(v)
  v = v - 0.5
  v = v / norm2(v)

  prob%n   = n
  prob%lam = 1.0

  n_pass = 0
  n_fail = 0

  print "(A)", "======================================================================"
  print "(A)", " Jacobian validator proof-test  (1D Bratu, n = 12)"
  print "(A)", "======================================================================"
  print "(A)", ""
  print "(A)", "  Convention: PASS = method's verdict matches expectation"
  print "(A)", "              (correct J -> ok; broken J -> not ok)"
  print "(A)", ""

  call run_case("J correct (control)",                BREAK_NONE,         &
    expect_fd = .true.,  expect_cs = .true.,  expect_ta = .true.,  expect_jvp = .true.)

  call run_case("J scaled (J' = 1.5 J)",              BREAK_SCALED,       &
    expect_fd = .false., expect_cs = .false., expect_ta = .false., expect_jvp = .false.)

  call run_case("J dropped diagonal at i=n/2",        BREAK_DROPPED_DIAG, &
    expect_fd = .false., expect_cs = .false., expect_ta = .false., expect_jvp = .false.)

  call run_case("J transposed off-diag pair at i=n/3", BREAK_TRANSPOSED,  &
    expect_fd = .false., expect_cs = .false., expect_ta = .false., expect_jvp = .false.)

  ! Noise at 1e-9 is intentionally just below FD's cancellation floor
  ! (~sqrt(eps) * ||F'|| ~ 1e-8) but well above CSD's machine-precision floor
  ! and just above JVP+Richardson's O(h^4) floor. The Taylor slope test is
  ! asymptotic and remains at 2.0 because R1(h) at the chosen h-window is
  ! dominated by the residual curvature, not the injected noise.
  ! Net effect: only FD is too coarse to reliably flag this defect.
  call run_case("J + noise (amplitude 1e-9)",          BREAK_NOISE,        &
    expect_fd = .true.,  expect_cs = .false., expect_ta = .true.,  expect_jvp = .false.)

  print "(A)", ""
  print "(A)", "======================================================================"
  print "(A)", " Paired-flag harness: exp_term toggle must vanish in BOTH F and J"
  print "(A)", "======================================================================"
  call run_paired_flag_test(prob, x)

  print "(A)", ""
  print "(A)", "======================================================================"
  print "(A)", " Snapshot dump: producing JSON+bin pairs for the Python jacval tool"
  print "(A)", "======================================================================"
  call dump_demo_snapshots(prob, x)

  print "(A)", "======================================================================"
  print "(A,I0,A,I0,A)", " Summary: ", n_pass, " passed, ", n_fail, " failed"
  print "(A)", "======================================================================"

  if (n_fail > 0) error stop 1

contains

  subroutine run_paired_flag_test(prob, x)
    !! Acceptance gate for §3.3 of jacval.md:
    !!   1. norm(F_on - F_off) > eps_term
    !!   2. norm(J_on - J_off) > eps_term
    !!   3. J_off validates against CSD with exp_term=false
    !!   4. J_on  validates against CSD with exp_term=true
    !! and one negative case:
    !!   5. BREAK_TERM_JAC_STALE = exp_term off in F but on in J  ==>  must FAIL
    type(bratu_problem), intent(inout) :: prob
    real,                intent(in)    :: x(:)

    real, allocatable :: F_on(:), F_off(:), J_on(:,:), J_off(:,:)
    real :: norm_F_diff, norm_J_diff, max_err
    logical :: ok_off, ok_on, ok_stale
    integer :: n, original_break

    n = size(x)
    allocate (F_on(n), F_off(n), J_on(n,n), J_off(n,n))

    original_break = prob%break
    prob%break = BREAK_NONE

    prob%exp_term = .true.
    call prob%eval_real(x, F_on)
    call prob%eval_jac (x, J_on)
    call validate_complex_step(prob, x, ok_on,  max_err)

    prob%exp_term = .false.
    call prob%eval_real(x, F_off)
    call prob%eval_jac (x, J_off)
    call validate_complex_step(prob, x, ok_off, max_err)

    norm_F_diff = norm2(F_on - F_off)
    norm_J_diff = sqrt(sum((J_on - J_off)**2))

    call report("F contributes      ", norm_F_diff > 1.0e-6, .true., norm_F_diff)
    call report("J contributes      ", norm_J_diff > 1.0e-6, .true., norm_J_diff)
    call report("CSD ok @ off       ", ok_off,               .true., max_err)
    call report("CSD ok @ on        ", ok_on,                .true., max_err)

    ! Negative case: stale Jacobian must be caught.
    prob%exp_term = .false.
    prob%break    = BREAK_TERM_JAC_STALE
    call validate_complex_step(prob, x, ok_stale, max_err)
    call report("stale-J caught     ", ok_stale,             .false., max_err)

    prob%break    = original_break
    prob%exp_term = .true.

    deallocate (F_on, F_off, J_on, J_off)
  end subroutine

  subroutine dump_demo_snapshots(prob, x)
    !! Generate three snapshots that the Python jacval CLI can act on:
    !!   bratu_correct.{json,bin}        :  exp_term=on,  no defect
    !!   bratu_dropped_diag.{json,bin}   :  exp_term=on,  BREAK_DROPPED_DIAG
    !!   bratu_exp_off.{json,bin}        :  exp_term=off, no defect
    !!
    !! With these in hand, toggle_diff(bratu_correct, bratu_exp_off) localizes
    !! the exp-term contribution, and inspect(bratu_dropped_diag) shows the
    !! dropped-diagonal bug ranked at the top.
    type(bratu_problem), intent(inout) :: prob
    real,                intent(in)    :: x(:)

    type(index_map_t)     :: im
    type(physics_flags_t) :: flags
    integer               :: i, n, original_break
    logical               :: original_exp

    n = size(x)
    allocate (im%equation(n), im%cell(n), im%labels(1))
    im%labels(1) = "x"
    do i = 1, n
      im%equation(i) = 1
      im%cell(i)     = i
    end do

    original_break  = prob%break
    original_exp    = prob%exp_term

    flags = physics_flags_t()
    flags%thermionic    = .false.   ! Bratu has no thermionic etc.; only one
    flags%tunneling     = .false.   ! toggled term here, modeled as "ifbl".
    flags%drift_contact = .false.
    flags%srh           = .false.
    flags%auger         = .false.

    prob%exp_term = .true.
    prob%break    = BREAK_NONE
    flags%ifbl = .true.
    call jacval_snapshot_dump("bratu_correct", prob, x, im, flags,           &
                              bias_V=0.0, build_hash="bratu_fixture")

    prob%exp_term = .true.
    prob%break    = BREAK_DROPPED_DIAG
    flags%ifbl = .true.
    call jacval_snapshot_dump("bratu_dropped_diag", prob, x, im, flags,      &
                              bias_V=0.0, build_hash="bratu_fixture")

    prob%exp_term = .false.
    prob%break    = BREAK_NONE
    flags%ifbl = .false.
    call jacval_snapshot_dump("bratu_exp_off", prob, x, im, flags,           &
                              bias_V=0.0, build_hash="bratu_fixture")

    prob%break    = original_break
    prob%exp_term = original_exp
  end subroutine

  subroutine run_case(label, break_mode, expect_fd, expect_cs, expect_ta, expect_jvp)
    character(*), intent(in) :: label
    integer,      intent(in) :: break_mode
    logical,      intent(in) :: expect_fd, expect_cs, expect_ta, expect_jvp

    logical :: ok_fd, ok_cs, ok_ta, ok_jvp
    real    :: err_fd, err_cs, err_jvp, slope_ta, slope_baseline

    real, allocatable :: J_claim(:,:), J_ref(:,:)
    integer           :: nn

    prob%break = break_mode

    print "(A)", ""
    print "(A,A)", "  Case: ", label
    print "(A)", "  ----------------------------------------------------------------"

    nn = size(x)
    allocate (J_claim(nn,nn), J_ref(nn,nn))

    call validate_finite_diff    (prob, x, ok_fd,  err_fd)
    call validate_complex_step   (prob, x, ok_cs,  err_cs, &
      &                            J_claim_out = J_claim, J_ref_out = J_ref)
    call validate_taylor         (prob, x, v, ok_ta,  slope_ta, slope_baseline)
    call validate_jvp_richardson (prob, x, v, ok_jvp, err_jvp)

    call report      ("finite_diff   ", ok_fd,  expect_fd,  err_fd)
    call report      ("complex_step  ", ok_cs,  expect_cs,  err_cs)
    call print_top_defects(J_claim, J_ref, k = 5, label = "CSD entrywise")
    call report_slope("taylor_slope  ", ok_ta,  expect_ta,  slope_ta, slope_baseline)
    call report      ("jvp_richardson", ok_jvp, expect_jvp, err_jvp)

    deallocate (J_claim, J_ref)
  end subroutine

  subroutine report(name, ok, expect, value)
    character(*), intent(in) :: name
    logical,      intent(in) :: ok, expect
    real,         intent(in) :: value
    character(:), allocatable :: verdict
    if (ok .eqv. expect) then
      n_pass = n_pass + 1
      verdict = "PASS"
    else
      n_fail = n_fail + 1
      verdict = "FAIL"
    end if
    print "(A,A,A,L1,A,L1,A,ES12.5,A,A,A)", &
      "    ", name, "  ok=", ok, "  expect=", expect, "  diag=", value, "   [", verdict, "]"
  end subroutine

  subroutine report_slope(name, ok, expect, slope, base)
    character(*), intent(in) :: name
    logical,      intent(in) :: ok, expect
    real,         intent(in) :: slope, base
    character(:), allocatable :: verdict
    if (ok .eqv. expect) then
      n_pass = n_pass + 1
      verdict = "PASS"
    else
      n_fail = n_fail + 1
      verdict = "FAIL"
    end if
    print "(A,A,A,L1,A,L1,A,F5.2,A,F5.2,A,A,A)", &
      "    ", name, "  ok=", ok, "  expect=", expect, "  slope=", slope, &
      "  baseline=", base, "   [", verdict, "]"
  end subroutine

end program
