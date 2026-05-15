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

  subroutine validate_finite_diff(prob, x, ok, max_err, tol, h)
    !! Central finite differences as numerical reference.
    class(jacobian_problem), intent(in)            :: prob
    real,                    intent(in)            :: x(:)
    logical,                 intent(out)           :: ok
    real,                    intent(out)           :: max_err
    real,                    intent(in),  optional :: tol
    real,                    intent(in),  optional :: h

    real, allocatable :: J_claim(:,:), J_fd(:,:), xp(:), xm(:), fp(:), fm(:)
    real    :: hloc, tolloc
    integer :: n, j

    n = size(x)
    if (present(h)) then
      hloc = h
    else
      ! optimum for central FD on real64: trades O(h^2) truncation against eps/h roundoff
      hloc = sqrt(epsilon(1.0)) * (1.0 + maxval(abs(x)))
    end if
    if (present(tol)) then
      tolloc = tol
    else
      tolloc = 1.0e-6
    end if

    allocate (J_claim(n,n), J_fd(n,n))
    allocate (xp(n), xm(n), fp(n), fm(n))

    call prob%eval_jac(x, J_claim)

    do j = 1, n
      xp = x; xp(j) = xp(j) + hloc
      xm = x; xm(j) = xm(j) - hloc
      call prob%eval_real(xp, fp)
      call prob%eval_real(xm, fm)
      J_fd(:,j) = (fp - fm) / (2.0 * hloc)
    end do

    max_err = maxval(abs(J_claim - J_fd))
    ok      = max_err <= tolloc
  end subroutine

  subroutine validate_complex_step(prob, x, ok, max_err, tol, h)
    !! Complex-step derivative reference: J_ij = Im(F_i(x + i h e_j))/h.
    class(jacobian_problem_cmplx), intent(in)            :: prob
    real,                          intent(in)            :: x(:)
    logical,                       intent(out)           :: ok
    real,                          intent(out)           :: max_err
    real,                          intent(in),  optional :: tol
    real,                          intent(in),  optional :: h

    real,    allocatable :: J_claim(:,:), J_cs(:,:)
    complex, allocatable :: z(:), fz(:)
    real    :: hloc, tolloc
    integer :: n, j

    n = size(x)
    if (present(h)) then
      hloc = h
    else
      hloc = 1.0e-200_real64  ! safe: no cancellation; only floor is overflow in F
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
