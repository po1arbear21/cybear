program max_test

  use error_m,         only: program_error
  use ieee_arithmetic, only: ieee_is_finite
  use gauss_m,         only: gauss_legendre
  use math_m,          only: ber, linspace, logspace
  use quad_m,          only: quad
  use util_m,          only: int2str

  implicit none

  integer, parameter :: NMAX = 256

  integer           :: N_ws, N_ts, N
  real              :: eta1, eta2, dpot, j, err_ws, err_ts
  real, allocatable :: x(:,:), w(:,:)

  allocate (x(NMAX,NMAX), w(NMAX,NMAX))
  !$omp parallel do schedule(dynamic) default(none) private(N) shared(x,w)
  do N = 1, NMAX
    print "(A,I0)", "Gauss-Legendre: ", N
    call gauss_legendre(x(1:N,N), w(1:N,N))
  end do
  !$omp end parallel do
  print *

  eta1 = -10.0
  eta2 = -5.0
  dpot = 3.0
  call scharfetter_gummel(eta1, eta2, dpot, j)
  print *, "j1 = ", j
  call output_integrand(eta1, eta2, dpot, j, "u1.csv")
  call test_weierstrass(eta1, eta2, dpot, j, 1, NMAX, "err_ws1.csv", "2", "1")
  call test_tanh_sinh(eta1, eta2, dpot, j, 1e-1, 1e-16, NMAX, "err_ts1.csv", "4", "2")

  eta1 = -10.0
  eta2 = -5.0
  dpot = -7.0
  call scharfetter_gummel(eta1, eta2, dpot, j)
  print *, "j2 = ", j
  call output_integrand(eta1, eta2, dpot, j, "u2.csv")
  call test_weierstrass(eta1, eta2, dpot, j, 1, NMAX, "err_ws2.csv", "2", "3")
  call test_tanh_sinh(eta1, eta2, dpot, j, 1e-1, 1e-16, NMAX, "err_ts2.csv", "4", "4")

contains

  subroutine dist(eta, F)
    real, intent(in)  :: eta
    real, intent(out) :: F

    F = exp(eta)
  end subroutine

  subroutine scharfetter_gummel(eta1, eta2, dpot, j)
    real, intent(in)  :: eta1
    real, intent(in)  :: eta2
    real, intent(in)  :: dpot
    real, intent(out) :: j

    real :: F1, F2, B1, B2

    call dist(eta1, F1)
    call dist(eta2, F2)
    B1 = ber( dpot)
    B2 = ber(-dpot)
    j  = B2 * F1 - B1 * F2
  end subroutine

  subroutine weierstrass(eta1, eta2, dpot, j, x, w, f, dfdj)
    real,    intent(in)  :: eta1
    real,    intent(in)  :: eta2
    real,    intent(in)  :: dpot
    real,    intent(in)  :: j
    real,    intent(in)  :: x(:)
    real,    intent(in)  :: w(:)
    real,    intent(out) :: f
    real,    intent(out) :: dfdj

    integer :: k
    real    :: eta, u, dudeta, dudp(2)

    f    = 0
    dfdj = 0
    do k = 1, size(x)
      eta = 0.5 * (eta1 + eta2) + 0.5 * (eta2 - eta1) * x(k)
      call my_integrand(eta, [dpot, j], u, dudeta, dudp)
      f = f + 0.5 * (eta2 - eta1) * w(k) * u
      dfdj = dfdj + 0.5 * (eta2 - eta1) * w(k) * dudp(2)
    end do
    f = f - 1
  end subroutine

  subroutine solve_weierstrass(eta1, eta2, dpot, j, N)
    real,    intent(in)    :: eta1
    real,    intent(in)    :: eta2
    real,    intent(in)    :: dpot
    real,    intent(inout) :: j
    integer, intent(in)    :: N

    real,    parameter :: ATOL = 1e-16, RTOL = 1e-14
    integer, parameter :: MIN_IT = 5, MAX_IT = 50

    integer           :: it
    real              :: err, f, dfdj, dj

    ! Newton iteration
    print "(A,I0)", "weierstrass: ", N
    err = huge(1.0)
    it  = 0
    do while ((it < MIN_IT) .or. (err > RTOL * abs(j) + ATOL))
      it = it + 1
      if (it > MAX_IT) then
        j = 0
        j = j / j
        return
      end if

      ! evaluate residual and get Newton update
      call weierstrass(eta1, eta2, dpot, j, x(1:N,N), w(1:N,N), f, dfdj)

      ! Newton update
      dj  = - f / dfdj
      err = abs(dj)

      ! update solution
      j = j + dj

      print "(I6,A,ES25.16E3,A,ES25.16E3,A,ES25.16E3)", it, ": j = ", j, "  +/-", err, "  ,", err / j
    end do
    print *
  end subroutine

  subroutine test_weierstrass(eta1, eta2, dpot, j, N1, N2, fname, mt, mc)
    real,         intent(in) :: eta1
    real,         intent(in) :: eta2
    real,         intent(in) :: dpot
    real,         intent(in) :: j
    integer,      intent(in) :: N1
    integer,      intent(in) :: N2
    character(*), intent(in) :: fname
    character(*), intent(in) :: mt
    character(*), intent(in) :: mc

    integer :: N, funit
    real    :: err, j1

    open (newunit = funit, file = fname, status = "replace", action = "write")
    write (funit, "(A)") "% lt=0 mt="//mt//" mc="//mc//" ms=4"
    do N = N1, N2
      j1 = j
      call solve_weierstrass(eta1, eta2, dpot, j1, N)
      err = abs(j - j1) / abs(j)
      if (.not. ieee_is_finite(err)) err = 1e3
      write (funit, "(I0, ES25.16E3)") N, max(err, 1e-16)
    end do
    close (funit)
  end subroutine

  subroutine tanh_sinh(eta1, eta2, dpot, j, eps, N, f, dfdj)
    real,    intent(in)  :: eta1
    real,    intent(in)  :: eta2
    real,    intent(in)  :: dpot
    real,    intent(in)  :: j
    real,    intent(in)  :: eps
    integer, intent(out) :: N
    real,    intent(out) :: f
    real,    intent(out) :: dfdj

    real :: dfdeta1, dfdeta2, dfdp(2)

    call quad(my_integrand, eta1, eta2, [dpot, j], f, dfdeta1, dfdeta2, dfdp, eps = eps, ncalls = N)
    f    = f - 1
    dfdj = dfdp(2)
  end subroutine

  subroutine solve_tanh_sinh(eta1, eta2, dpot, j, eps, N)
    real,    intent(in)    :: eta1
    real,    intent(in)    :: eta2
    real,    intent(in)    :: dpot
    real,    intent(inout) :: j
    real,    intent(in)    :: eps
    integer, intent(out)   :: N

    real,    parameter :: ATOL = 1e-16, RTOL = 1e-14
    integer, parameter :: MIN_IT = 5, MAX_IT = 20

    integer :: it
    real    :: err, f, dfdj, dj

    ! Newton iteration
    print "(A,ES25.16E3)", "tanh_sinh: ", eps
    err = huge(1.0)
    it  = 0
    do while ((it < MIN_IT) .or. (err > RTOL * abs(j) + ATOL))
      it = it + 1
      if (it > MAX_IT) call program_error("No convergence after " // int2str(MAX_IT) // " iterations")

      ! evaluate residual and get Newton update
      call tanh_sinh(eta1, eta2, dpot, j, eps, N, f, dfdj)

      ! Newton update
      dj  = - f / dfdj
      err = abs(dj)

      ! update solution
      j = j + dj

      print "(I6,A,ES25.16E3,A,ES25.16E3,A,ES25.16E3)", it, ": j = ", j, "  +/-", err, "  ,", err / j
    end do
  end subroutine

  subroutine test_tanh_sinh(eta1, eta2, dpot, j, eps1, eps2, Neps, fname, mt, mc)
    real,         intent(in) :: eta1
    real,         intent(in) :: eta2
    real,         intent(in) :: dpot
    real,         intent(in) :: j
    real,         intent(in) :: eps1
    real,         intent(in) :: eps2
    integer,      intent(in) :: Neps
    character(*), intent(in) :: fname
    character(*), intent(in) :: mt
    character(*), intent(in) :: mc

    integer           :: i, N, funit
    real              :: j1, err
    real, allocatable :: eps(:)

    eps = logspace(eps1, eps2, Neps)

    open (newunit = funit, file = fname, status = "replace", action = "write")
    write (funit, "(A)") "% lt=0 mt="//mt//" mc="//mc//" ms=4"
    do i = 1, Neps
      j1 = j
      call solve_tanh_sinh(eta1, eta2, dpot, j1, eps(i), N)
      err = abs(j - j1) / abs(j)
      if (.not. ieee_is_finite(err)) err = 1e3
      write (funit, "(I0, ES25.16E3)") N, max(err, 1e-16)
    end do
    close (funit)
  end subroutine

  subroutine my_integrand(eta, p, u, dudeta, dudp)
    real, intent(in)  :: eta
    real, intent(in)  :: p(:)
    real, intent(out) :: u
    real, intent(out) :: dudeta
    real, intent(out) :: dudp(:)

    real :: dpot, j, F

    dpot = p(1)
    j    = p(2)
    call dist(eta, F)

    u       = 1.0 / (dpot - j / F)
    dudeta  = 0.0
    dudp(1) = 0.0
    dudp(2) = 1.0 / (F * (dpot - j / F)**2)
  end subroutine

  subroutine output_integrand(eta1, eta2, dpot, j, fname)
    real,         intent(in) :: eta1
    real,         intent(in) :: eta2
    real,         intent(in) :: dpot
    real,         intent(in) :: j
    character(*), intent(in) :: fname

    integer           :: funit, i
    real, allocatable :: eta(:)
    real              :: u, dudeta, dudp(2)

    eta = linspace(eta1, eta2, 1001)
    open (newunit = funit, file = fname, status = "replace", action = "write")
    do i = 1, size(eta)
      call my_integrand(eta(i), [dpot, j], u, dudeta, dudp)
      write (funit, "(2ES25.16E3)") eta(i), u
    end do
    close (funit)
  end subroutine

end program
