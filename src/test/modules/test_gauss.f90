m4_include(../../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_mpfr},0,-1))

module test_gauss_m

  use gauss_m,     only: gauss
  use math_m,      only: PI
  use mpfr_m,      only: div, mpfr
  use test_case_m, only: test_case

  implicit none

  private
  public test_gauss

contains

  subroutine test_gauss()
    integer, parameter :: N = 8
    real,    parameter :: rtol = 1e-15, atol = 1e-16

    type(test_case) :: tc
    integer         :: i, j
    real            :: x(N), w(N), f, FF

    call tc%init("gauss")

    ! test Gauss-Legendre
    call gauss("legendre", x, w)
    FF = 0
    do i = 1, N
      f = 1.0
      do j = 5, 1, -1
        f = j + x(i) * f
      end do
      FF = FF + w(i) * f
    end do
    call tc%assert_eq(6.0, FF, rtol, atol, "Gauss-Legendre")

    ! test Gauss-Laguerre
    call gauss("laguerre", x, w)
    FF = 0
    do i = 1, N
      f = 1.0
      do j = 5, 1, -1
        f = j + x(i) * f
      end do
      FF = FF + w(i) * f
    end do
    call tc%assert_eq(273.0, FF, rtol, atol, "Gauss-Laguerre")

    ! test Gauss-Hermite
    call gauss("hermite", x, w)
    FF = 0
    do i = 1, N
      f = 1.0
      do j = 5, 1, -1
        f = j + x(i) * f
      end do
      FF = FF + w(i) * f
    end do
    call tc%assert_eq(25 * sqrt(PI) / 4, FF, rtol, atol, "Gauss-Hermite")

    ! test Gauss quadrature with custom weight function (sqrt(x) in interval 0 to 1)
    call gauss("custom", x, w, momfun = moments)
    FF = 0
    do i = 1, N
      f = 1.0
      do j = 5, 1, -1
        f = j + x(i) * f
      end do
      FF = FF + w(i) * f
    end do
    call tc%assert_eq(192596.0 / 45045.0, FF, rtol, atol, "Gauss sqrt(x)")

    call tc%finish()
  end subroutine

  subroutine moments(s)
    !! moments for weight function sqrt(x) in interval 0 to 1
    type(mpfr), intent(inout) :: s(0:)

    integer :: i

    ! s_n = integral_0^1 sqrt(x) * x**n dx = 2 / (2 * n + 3)
    do i = 0, ubound(s, 1)
      call s(i)%set(2 * i + 3)
      call div(s(i), 2.0, s(i))
    end do
  end subroutine

end module

m4_divert(0)
