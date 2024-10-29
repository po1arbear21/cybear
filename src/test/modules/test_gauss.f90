m4_include(../../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_mpfr},0,-1))

module test_gauss_m

  use gauss_m,     only: gauss, gauss_legendre, gauss_laguerre, gauss_hermite
  use math_m,      only: PI
  use mpfr_m,      only: add, div, mpfr, mul, sqr
  use test_case_m, only: test_case

  implicit none

  private
  public test_gauss

contains

  subroutine test_gauss()
    integer, parameter :: N = 8
    real,    parameter :: rtol = 1e-15, atol = 1e-16

    integer         :: i, j
    real            :: x(N), w(N), s, dxdp(N,1), dwdp(N,1), dsdp(1), f, dfdp, FF, dFFdp
    type(gauss)     :: gs
    type(test_case) :: tc

    ! test polynomial: f(x) = 1 + x*(2 + x*(3 + x*(4 + x*(5 + x))))

    call tc%init("gauss")

    ! test Gauss-Legendre
    call gauss_legendre(x, w)
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
    call gauss_laguerre(x, w)
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
    call gauss_hermite(x, w)
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
    call gs%init(N, NP=1)
    call gs%generate(moments, [1.0], x, w, s, dxdp, dwdp, dsdp)
    FF    = 0
    dFFdp = 0
    do i = 1, N
      f    = 1.0
      dfdp = 0.0
      do j = 5, 1, -1
        dfdp = dxdp(i,1) * f + x(i) * dfdp
        f    = j + x(i) * f
      end do
      FF    = FF + w(i) * s * f
      dFFdp = dFFdp + dwdp(i,1) * s * f + w(i) * dsdp(1) * f + w(i) * s * dfdp
    end do
    call tc%assert_eq(11.0 / 6.0, FF, rtol, atol, "Gauss custom 1")
    call tc%assert_eq(83.0 / 72.0, dFFdp, rtol, atol, "Gauss custom 1 d/dp")

    call gs%generate(moments, [1.5], x, w, s, dxdp, dwdp, dsdp)
    FF    = 0
    dFFdp = 0
    do i = 1, N
      f    = 1.0
      dfdp = 0.0
      do j = 5, 1, -1
        dfdp = dxdp(i,1) * f + x(i) * dfdp
        f    = j + x(i) * f
      end do
      FF    = FF + w(i) * s * f
      dFFdp = dFFdp + dwdp(i,1) * s * f + w(i) * dsdp(1) * f + w(i) * s * dfdp
    end do
    call tc%assert_eq(454.0 / 195.0, FF, rtol, atol, "Gauss custom 2")
    call tc%assert_eq(97304.0 / 114075.0, dFFdp, rtol, atol, "Gauss custom 2 d/dp")

    call gs%destruct()

    call tc%finish()
  end subroutine

  subroutine moments(p, s, dsdp)
    !! moments for weight function 1 - abs(x)**p in interval -1 to 1
    real,       intent(in)    :: p(:)
      !! parameters
    type(mpfr), intent(inout) :: s(0:)
      !! moments, already initialized
    type(mpfr), intent(inout) :: dsdp(0:,:)
      !! derivatives of s wrt p

    integer :: k, n

    n = (ubound(s, 1) + 1) / 2

    ! s_{2k}   = integral_{-1}^1 (1 - abs(x)**p) * x**{2k} dx = 2*p/((2k+1)*(p+2k+1))
    ! s_{2k+1} = 0
    do k = 0, n - 1
      ! s    = 2*p/((2*k+1)*(2*k+1 + p))
      ! dsdp = 2 / (2*k+1 + p)**2
      call s(2*k)%set(p(1))
      call add(s(2*k), s(2*k), 2*k+1)
      call dsdp(2*k,1)%set(s(2*k))
      call mul(s(2*k), s(2*k), 2*k+1)
      call div(s(2*k), p(1), s(2*k))
      call add(s(2*k), s(2*k), s(2*k))
      call sqr(dsdp(2*k,1), dsdp(2*k,1))
      call div(dsdp(2*k,1), 2, dsdp(2*k,1))

      call s(2*k+1)%set(0)

      call dsdp(2*k+1,1)%set(0)
    end do
  end subroutine

end module

m4_divert(0)
