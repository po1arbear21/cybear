m4_include(../../util/macro.f90.inc)

module test_quad_m

  use ieee_arithmetic, only: ieee_value, ieee_positive_inf
  use math_m,          only: PI
  use quad_m,          only: quad
  use test_case_m,     only: test_case

  implicit none

  private
  public test_quad

contains

  subroutine test_quad()
    real            :: I, dIda, dIdb, dIdp(2), inf, dum(0), dum2(0)
    type(test_case) :: tc

    call tc%init("quad")

    ! 3 * x**2
    call quad(integrand1, 1.0, 2.0, [3.0], I, dIda, dIdb, dIdp(1:1))
    call tc%assert_eq(7.0, I, 1e-14, 1e-16, "I1")
    call tc%assert_eq(-3.0, dIda, 1e-14, 1e-16, "dI1da")
    call tc%assert_eq(12.0, dIdb, 1e-14, 1e-16, "dI1db")
    call tc%assert_eq(7.0 / 3.0, dIdp(1), 1e-14, 1e-16, "dI1dp")

    ! singularity close to 1
    call quad(integrand2, -1.0, 1.0, [2.0, 1.015625], I, dIda, dIdb, dIdp, max_levels = 16)
    call tc%assert_eq(-9.7196248087233439,     I, 1e-14, 1e-16, "I2")
    call tc%assert_eq( 0.99224806201550386, dIda, 1e-14, 1e-16, "dI2da")
    call tc%assert_eq(-128.0,               dIdb, 1e-14, 1e-16, "dI2db")
    call tc%assert_eq([-4.8598124043616721, 127.007751937984496], dIdp, 1e-14, 1e-16, "dI2dp")

    ! singularity at 0
    call quad(integrand3, 0.00390625, 10.0, dum, I, dIda, dIdb, dum2)
    call tc%assert_eq(5.5470845327363952, I, 1e-14, 1e-16, "I3")
    call tc%assert_eq(-255.50032552075055, dIda, 1e-14, 1e-16, "dI3da")
    call tc%assert_eq(0.000045401991009687768, dIdb, 1e-14, 1e-16, "dI3db")

    ! ! lower bound > upper bound
    call quad(integrand3, 10.0, 0.00390625, dum, I, dIda, dIdb, dum2)
    call tc%assert_eq(-5.5470845327363952, I, 1e-14, 1e-16, "I3a")
    call tc%assert_eq(-0.000045401991009687768, dIda, 1e-14, 1e-16, "dI3da")
    call tc%assert_eq(255.50032552075055, dIdb, 1e-14, 1e-16, "dI3db")

    ! int_{0.5}^{inf} 1.5 * exp(-x/4) dx
    inf = ieee_value(0.0, ieee_positive_inf)
    call quad(integrand4, 0.5, inf, [1.5, 0.25], I, dIda, dIdb, dIdp)
    call tc%assert_eq(6*exp(-0.125), I, 1e-14, 1e-16, "I4")
    call tc%assert_eq(-1.5*exp(-0.125), dIda, 1e-14, 1e-16, "dI4da")
    call tc%assert_eq([4*exp(-0.125), -27*exp(-0.125)], dIdp, 1e-14, 1e-16, "dI4dp")

    ! int_{-inf}^{inf} exp(-2*x^2) dx
    call quad(integrand5, -inf, inf, [2.0], I, dIda, dIdb, dIdp(1:1))
    call tc%assert_eq(sqrt(PI/2), I, 1e-14, 1e-16, "I5")
    call tc%assert_eq(-0.25*sqrt(PI/2), dIdp(1), 1e-14, 1e-16, "dI5dp")

    call tc%finish()
  end subroutine

  subroutine integrand1(x, p, f, dfdx, dfdp)
    real, intent(in)  :: x
    real, intent(in)  :: p(:)
    real, intent(out) :: f
    real, intent(out) :: dfdx
    real, intent(out) :: dfdp(:)

    f       = p(1) * x * x
    dfdx    = 2 * p(1) * x
    dfdp(1) = x * x
  end subroutine

  subroutine integrand2(x, p, f, dfdx, dfdp)
    real, intent(in)  :: x
    real, intent(in)  :: p(:)
    real, intent(out) :: f
    real, intent(out) :: dfdx
    real, intent(out) :: dfdp(:)

    f       = p(1) / (x - p(2))
    dfdx    = - p(1) / (x - p(2))**2
    dfdp(1) = 1 / (x - p(2))
    dfdp(2) = p(1) / (x - p(2))**2
  end subroutine

  subroutine integrand3(x, p, f, dfdx, dfdp)
    use math_m, only: expm1

    real, intent(in)  :: x
    real, intent(in)  :: p(:)
    real, intent(out) :: f
    real, intent(out) :: dfdx
    real, intent(out) :: dfdp(:)

    m4_ignore(p)
    m4_ignore(dfdp)

    f = 1 / expm1(x)
    dfdx = - exp(x) * f**2
  end subroutine

  subroutine integrand4(x, p, f, dfdx, dfdp)
    real, intent(in)  :: x
    real, intent(in)  :: p(:)
    real, intent(out) :: f
    real, intent(out) :: dfdx
    real, intent(out) :: dfdp(:)

    f       =   p(1)        * exp(-p(2) * x)
    dfdx    = - p(1) * p(2) * exp(-p(2) * x)
    dfdp(1) =                 exp(-p(2) * x)
    dfdp(2) = - p(1) *    x * exp(-p(2) * x)
  end subroutine

  subroutine integrand5(x, p, f, dfdx, dfdp)
    real, intent(in)  :: x
    real, intent(in)  :: p(:)
    real, intent(out) :: f
    real, intent(out) :: dfdx
    real, intent(out) :: dfdp(:)

    f       =                  exp(-p(1) * x**2)
    dfdx    = - 2 * p(1) * x * exp(-p(1) * x**2)
    dfdp(1) = - x**2         * exp(-p(1) * x**2)
  end subroutine

end module
