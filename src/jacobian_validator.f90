module jacobian_validator_m
  !! Validate an analytical Jacobian J(x) of a nonlinear vector residual F(x).
  !!
  !! Four complementary tests are provided:
  !!
  !!   1. validate_finite_diff       Baseline central FD, O(h^2) truncation.
  !!                                 Cancellation floor ~ sqrt(eps).
  !!   2. validate_complex_step      J_ij = Im(F_i(x + i h e_j))/h. No subtractive
  !!                                 cancellation; verifies J to machine precision
  !!                                 when F admits a complex-arithmetic evaluator.
  !!   3. validate_taylor            Fits the log-log slope of the Taylor remainder
  !!                                 ||F(x+hv) - F(x) - h J v|| over a sweep of h.
  !!                                 Slope 2 iff J is correct; slope 1 if not.
  !!                                 Dimensionless and threshold-free.
  !!   4. validate_jvp_richardson    Matrix-free directional check. Compares J v
  !!                                 with a Richardson-extrapolated central FD of
  !!                                 F along v, accurate to O(h^4).
  !!
  !! Each method returns ok plus a numeric diagnostic. The Taylor test also
  !! reports the fitted slope so callers can inspect convergence order.

  use, intrinsic :: iso_fortran_env, only: real64

  implicit none

  private
  public :: jacobian_problem
  public :: jacobian_problem_cmplx
  public :: validate_finite_diff
  public :: validate_complex_step
  public :: validate_taylor
  public :: validate_jvp_richardson
  public :: log_log_slope
  public :: print_top_defects

  type, abstract :: jacobian_problem
    !! Residual + analytical Jacobian provider.
  contains
    procedure(eval_real_i), deferred :: eval_real
    procedure(eval_jac_i),  deferred :: eval_jac
  end type

  type, abstract, extends(jacobian_problem) :: jacobian_problem_cmplx
    !! Adds a complex-arithmetic residual evaluator (required by CSD).
  contains
    procedure(eval_cmplx_i), deferred :: eval_cmplx
  end type

  abstract interface
    subroutine eval_real_i(this, x, f)
      import jacobian_problem
      class(jacobian_problem), intent(in)  :: this
      real,                    intent(in)  :: x(:)
      real,                    intent(out) :: f(:)
    end subroutine

    subroutine eval_jac_i(this, x, J)
      import jacobian_problem
      class(jacobian_problem), intent(in)  :: this
      real,                    intent(in)  :: x(:)
      real,                    intent(out) :: J(:,:)
    end subroutine

    subroutine eval_cmplx_i(this, z, fz)
      import jacobian_problem_cmplx
      class(jacobian_problem_cmplx), intent(in)  :: this
      complex,                       intent(in)  :: z(:)
      complex,                       intent(out) :: fz(:)
    end subroutine
  end interface

contains

  subroutine validate_finite_diff(prob, x, ok, max_err, tol, h, h_vec, positive_mask, J_claim_out, J_ref_out)
    !! Central finite differences as numerical reference.
    !!
    !! Optional out-args expose the claimed and FD-reference Jacobians so the
    !! caller can run element-wise diagnostics (see print_top_defects). The
    !! FD reference carries a per-entry roundoff floor of ~sqrt(eps), which
    !! is far coarser than CSD -- so for ratio-based localization on FD,
    !! choose ratio_threshold accordingly (typically 1e-3..1e-4, not 1e-6).
    !!
    !! Step-size resolution (precedence): h_vec(j) > h > default.
    !!   - h_vec(:): per-column step; required for mixed-scale systems where a
    !!     single global step is too large for tiny variables and too small for
    !!     large ones (e.g. cybear DD: dens at depleted Schottky O(1e-14) vs
    !!     pot O(10)).
    !!   - h scalar override: same h for every column.
    !!   - default: sqrt(eps) * (1 + ||x||_inf); fine for uniformly-scaled
    !!     problems, broken for mixed-scale ones.
    !!
    !! positive_mask(:) optional: when .true. for column j, switch to the
    !! second-order forward stencil  ( -3 F(x) + 4 F(x+h) - F(x+2h) ) / (2h)
    !! whenever x(j) - h(j) would cross zero. For positive-only variables
    !! (carrier densities, concentrations) this avoids the catastrophic
    !! branch on x <= 0 in the residual (e.g. the schottky_tunneling
    !! early-return path) while keeping O(h^2) accuracy.
    !!
    !! On return, the routine restores the underlying eval state to x (calls
    !! eval_real(x, ...) once after the FD loop). Callers can rely on the
    !! problem being at baseline x after this routine returns.
    class(jacobian_problem), intent(in)            :: prob
    real,                    intent(in)            :: x(:)
    logical,                 intent(out)           :: ok
    real,                    intent(out)           :: max_err
    real,                    intent(in),  optional :: tol
    real,                    intent(in),  optional :: h
    real,                    intent(in),  optional :: h_vec(:)
    logical,                 intent(in),  optional :: positive_mask(:)
    real,                    intent(out), optional :: J_claim_out(:,:)
    real,                    intent(out), optional :: J_ref_out(:,:)

    real, allocatable :: J_claim(:,:), J_fd(:,:), xp(:), xm(:), x2(:), fp(:), fm(:), f2(:)
    real, allocatable :: hcol(:), f0(:)
    logical, allocatable :: posmask(:)
    real    :: hloc, tolloc
    integer :: n, j
    logical :: use_forward

    n = size(x)
    allocate (hcol(n), posmask(n))
    if (present(h_vec)) then
      hcol = h_vec
    else
      if (present(h)) then
        hloc = h
      else
        hloc = sqrt(epsilon(1.0)) * (1.0 + maxval(abs(x)))
      end if
      hcol = hloc
    end if
    if (present(positive_mask)) then
      posmask = positive_mask
    else
      posmask = .false.
    end if
    if (present(tol)) then
      tolloc = tol
    else
      tolloc = 1.0e-6
    end if

    allocate (J_claim(n,n), J_fd(n,n))
    allocate (xp(n), xm(n), x2(n), fp(n), fm(n), f2(n), f0(n))

    call prob%eval_jac(x, J_claim)
    call prob%eval_real(x, f0)

    do j = 1, n
      ! Switch to second-order forward stencil for positive-only variables
      ! whose backward step would cross zero. The 3-point forward formula
      ! (-3f0 + 4f+ - f++)/(2h) preserves O(h^2) accuracy.
      use_forward = posmask(j) .and. (x(j) - hcol(j) <= 0.0)
      if (use_forward) then
        xp = x; xp(j) = xp(j) +       hcol(j)
        x2 = x; x2(j) = x2(j) + 2.0 * hcol(j)
        call prob%eval_real(xp, fp)
        call prob%eval_real(x2, f2)
        J_fd(:,j) = (-3.0 * f0 + 4.0 * fp - f2) / (2.0 * hcol(j))
      else
        xp = x; xp(j) = xp(j) + hcol(j)
        xm = x; xm(j) = xm(j) - hcol(j)
        call prob%eval_real(xp, fp)
        call prob%eval_real(xm, fm)
        J_fd(:,j) = (fp - fm) / (2.0 * hcol(j))
      end if
    end do

    ! Restore state to x so downstream code sees a clean baseline.
    call prob%eval_real(x, fp)

    max_err = maxval(abs(J_claim - J_fd))
    ok      = max_err <= tolloc

    if (present(J_claim_out)) J_claim_out = J_claim
    if (present(J_ref_out))   J_ref_out   = J_fd
  end subroutine

  subroutine validate_complex_step(prob, x, ok, max_err, tol, h, J_claim_out, J_ref_out)
    !! Complex-step derivative reference: J_ij = Im(F_i(x + i h e_j))/h.
    !!
    !! Optional out-args expose the claimed and reference Jacobians so the
    !! caller can run element-wise diagnostics (see print_top_defects). CSD
    !! has no subtractive cancellation, so J_ref is accurate to machine
    !! precision and an entrywise diff localizes the bug.

    class(jacobian_problem_cmplx), intent(in)            :: prob
    real,                          intent(in)            :: x(:)
    logical,                       intent(out)           :: ok
    real,                          intent(out)           :: max_err
    real,                          intent(in),  optional :: tol
    real,                          intent(in),  optional :: h
    real,                          intent(out), optional :: J_claim_out(:,:)
    real,                          intent(out), optional :: J_ref_out(:,:)

    real,    allocatable :: J_claim(:,:), J_cs(:,:)
    complex, allocatable :: z(:), fz(:)
    real    :: hloc, tolloc
    integer :: n, j

    n = size(x)
    if (present(h)) then
      hloc = h
    else
      ! h chosen far below sqrt(eps) so the h^2 truncation term in CSD is
      ! invisible in floating point; no subtractive cancellation either way.
      ! Exponent kept in single-precision range to avoid a spurious lexer
      ! underflow warning on the sp-parsed literal.
      hloc = 1.0e-30_real64
    end if
    if (present(tol)) then
      tolloc = tol
    else
      tolloc = 1.0e-10
    end if

    allocate (J_claim(n,n), J_cs(n,n))
    allocate (z(n), fz(n))

    call prob%eval_jac(x, J_claim)

    do j = 1, n
      z      = cmplx(x, 0.0)
      z(j)   = z(j) + cmplx(0.0, hloc)
      call prob%eval_cmplx(z, fz)
      J_cs(:,j) = aimag(fz) / hloc
    end do

    max_err = maxval(abs(J_claim - J_cs))
    ok      = max_err <= tolloc

    if (present(J_claim_out)) J_claim_out = J_claim
    if (present(J_ref_out))   J_ref_out   = J_cs
  end subroutine

  subroutine validate_taylor(prob, x, v, ok, slope, slope_baseline, hmin, hmax, nh, slope_tol)
    !! Taylor-remainder log-log slope test.
    !!
    !! R0(h) := ||F(x + h v) - F(x)||           expected slope ~ 1 (Lipschitz)
    !! R1(h) := ||F(x + h v) - F(x) - h J v||   expected slope ~ 2 iff J correct,
    !!                                          slope ~ 1 if J is wrong.
    !!
    !! Fits slope by ordinary least squares in log-log space; passes iff
    !! slope >= slope_tol AND R0's slope >= 0.9 (sanity check that F is not
    !! locally constant along v).
    class(jacobian_problem), intent(in)            :: prob
    real,                    intent(in)            :: x(:)
    real,                    intent(in)            :: v(:)
    logical,                 intent(out)           :: ok
    real,                    intent(out)           :: slope
    real,                    intent(out)           :: slope_baseline
    real,                    intent(in),  optional :: hmin, hmax
    integer,                 intent(in),  optional :: nh
    real,                    intent(in),  optional :: slope_tol

    real, allocatable :: J(:,:), f0(:), fh(:), Jv(:), hs(:), R0(:), R1(:)
    real    :: hmn, hmx, stol
    integer :: n, k, nh_

    n = size(x)
    if (present(hmin))      then; hmn = hmin;      else; hmn  = 1.0e-6; end if
    if (present(hmax))      then; hmx = hmax;      else; hmx  = 1.0e-2; end if
    if (present(nh))        then; nh_ = nh;        else; nh_  = 8;      end if
    if (present(slope_tol)) then; stol = slope_tol; else; stol = 1.8;   end if

    allocate (J(n,n), f0(n), fh(n), Jv(n))
    allocate (hs(nh_), R0(nh_), R1(nh_))

    call prob%eval_jac (x, J)
    call prob%eval_real(x, f0)
    Jv = matmul(J, v)

    do k = 1, nh_
      hs(k) = hmx * (hmn / hmx)**(real(k - 1) / real(nh_ - 1))
      call prob%eval_real(x + hs(k) * v, fh)
      R0(k) = norm2(fh - f0)
      R1(k) = norm2(fh - f0 - hs(k) * Jv)
    end do

    slope          = log_log_slope(hs, R1)
    slope_baseline = log_log_slope(hs, R0)

    ok = (slope >= stol) .and. (slope_baseline >= 0.9)

    ! Restore state to x so downstream probes see a clean baseline.
    call prob%eval_real(x, fh)
  end subroutine

  subroutine validate_jvp_richardson(prob, x, v, ok, err, tol, h)
    !! Matrix-free check: compare J v with Richardson-extrapolated central FD.
    !!
    !!   g(h)   = (F(x + h v) - F(x - h v)) / (2 h)            = J v + O(h^2)
    !!   g_R    = (4 g(h/2) - g(h)) / 3                        = J v + O(h^4)
    !!
    !! Four residual evaluations regardless of system size.
    class(jacobian_problem), intent(in)            :: prob
    real,                    intent(in)            :: x(:)
    real,                    intent(in)            :: v(:)
    logical,                 intent(out)           :: ok
    real,                    intent(out)           :: err
    real,                    intent(in),  optional :: tol
    real,                    intent(in),  optional :: h

    real, allocatable :: J(:,:), Jv(:), Jv_ref(:)
    real, allocatable :: fp1(:), fm1(:), fp2(:), fm2(:), g1(:), g2(:)
    real    :: hloc, tolloc
    integer :: n

    n = size(x)
    if (present(h)) then
      hloc = h
    else
      hloc = 1.0e-4
    end if
    if (present(tol)) then
      tolloc = tol
    else
      tolloc = 1.0e-9
    end if

    allocate (J(n,n), Jv(n), Jv_ref(n))
    allocate (fp1(n), fm1(n), fp2(n), fm2(n), g1(n), g2(n))

    call prob%eval_jac(x, J)
    Jv = matmul(J, v)

    call prob%eval_real(x + hloc       * v, fp1)
    call prob%eval_real(x - hloc       * v, fm1)
    g1 = (fp1 - fm1) / (2.0 * hloc)

    call prob%eval_real(x + 0.5 * hloc * v, fp2)
    call prob%eval_real(x - 0.5 * hloc * v, fm2)
    g2 = (fp2 - fm2) / (hloc)

    Jv_ref = (4.0 * g2 - g1) / 3.0

    err = norm2(Jv - Jv_ref)
    ok  = err <= tolloc

    ! Restore state to x so downstream probes see a clean baseline.
    call prob%eval_real(x, fp1)
  end subroutine

  subroutine print_top_defects(J_claim, J_ref, unit, k, ratio_threshold, floor_frac, label)
    !! Print the largest entrywise disagreements between a claimed Jacobian
    !! and a numerical reference (typically the CSD result).
    !!
    !! Score:
    !!   denom_ij = max(|J_ref(i,j)|, floor_frac * ||J_ref||_inf)
    !!   ratio_ij = |J_claim(i,j) - J_ref(i,j)| / denom_ij
    !! Entries with ratio_ij > ratio_threshold are flagged; up to k are
    !! printed, sorted by ratio descending.
    !!
    !! The norm-tied floor avoids two failure modes:
    !!   1. An eps-floored per-entry ratio false-positives all over the
    !!      off-stencil zeros, where CSD's machine-precision noise blows
    !!      past any reasonable ratio threshold.
    !!   2. A single global-norm denominator masks small-magnitude block
    !!      defects under whichever block dominates ||J||_inf -- in DD
    !!      that means Poisson hides Schottky / tunneling / recombination.
    real,         intent(in)           :: J_claim(:,:)
    real,         intent(in)           :: J_ref(:,:)
    integer,      intent(in), optional :: unit
    integer,      intent(in), optional :: k
    real,         intent(in), optional :: ratio_threshold
    real,         intent(in), optional :: floor_frac
    character(*), intent(in), optional :: label

    integer :: n, m, i, j, kk, kmax, u, n_flag, idx_max
    real    :: norm_inf, floor_val, denom, ratio, rt, ff

    integer, allocatable :: ii(:), jj(:)
    real,    allocatable :: rr(:)

    n = size(J_claim, 1)
    m = size(J_claim, 2)
    if (present(unit))            then; u    = unit;            else; u    = 6;       end if
    if (present(k))               then; kmax = k;               else; kmax = 5;       end if
    if (present(ratio_threshold)) then; rt   = ratio_threshold; else; rt   = 1.0e-6;  end if
    if (present(floor_frac))      then; ff   = floor_frac;      else; ff   = 1.0e-12; end if

    norm_inf  = maxval(abs(J_ref))
    floor_val = ff * norm_inf

    allocate (ii(n*m), jj(n*m), rr(n*m))
    n_flag = 0
    do j = 1, m
      do i = 1, n
        denom = max(abs(J_ref(i,j)), floor_val)
        ratio = abs(J_claim(i,j) - J_ref(i,j)) / denom
        if (ratio > rt) then
          n_flag     = n_flag + 1
          ii(n_flag) = i
          jj(n_flag) = j
          rr(n_flag) = ratio
        end if
      end do
    end do

    if (present(label)) then
      write (u, "(A,A,A)") "      defects [", trim(label), "]:"
    else
      write (u, "(A)")     "      defects:"
    end if
    write (u, "(A,ES10.3,A,ES10.3,A,ES10.3,A)") &
      "        ||J_ref||_inf = ", norm_inf,                                       &
      "   floor = ", floor_val, " (", ff, " * ||J_ref||_inf)"

    if (n_flag == 0) then
      write (u, "(A,ES10.3,A)") &
        "        no entries above ratio threshold ", rt, "   -- J matches J_ref entrywise."
      return
    end if

    write (u, "(A,I0,A,ES10.3,A)") &
      "        ", n_flag, " entries above ratio threshold ", rt, ":"

    do kk = 1, min(kmax, n_flag)
      idx_max = maxloc(rr(1:n_flag), dim = 1)
      write (u, "(A,I3,A,I3,A,ES12.4,A,ES12.4,A,ES10.3)") &
        "          (i=", ii(idx_max), ", j=", jj(idx_max),       &
        ")  J_claim=", J_claim(ii(idx_max), jj(idx_max)),        &
        "  J_ref=",    J_ref  (ii(idx_max), jj(idx_max)),        &
        "  ratio=",    rr(idx_max)
      rr(idx_max) = -1.0   ! consumed; next maxloc picks the next-largest
    end do

    if (n_flag > kmax) then
      write (u, "(A,I0,A)") "          ... (", n_flag - kmax, " more entries omitted)"
    end if
  end subroutine

  pure function log_log_slope(x, y) result(m)
    !! OLS slope of (log x, log y); used by validate_taylor.
    real, intent(in) :: x(:)
    real, intent(in) :: y(:)
    real             :: m

    real    :: lx(size(x)), ly(size(y)), mx, my, num, den
    integer :: n

    n  = size(x)
    lx = log(x)
    ly = log(max(y, tiny(1.0)))   ! guard against 0 from exact cancellation
    mx = sum(lx) / n
    my = sum(ly) / n
    num = sum((lx - mx) * (ly - my))
    den = sum((lx - mx)**2)
    m   = num / den
  end function

end module
