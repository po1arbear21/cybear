m4_include(../../util/macro.f90.inc)

module test_quad_m

  use ieee_arithmetic, only: ieee_value, ieee_positive_inf
  m4_ifdef({m4_mpfr},{
  use mpfr_m, only: div, exp, expm1, mpfr, mpfr_startup, mpfr_cleanup, neg, sqr, sub
  })
  use math_m,      only: PI
  use quad_m,      only: quad m4_ifdef({m4_mpfr},{, quad_mpfr})
  use string_m,    only: string
  use test_case_m, only: test_case

  implicit none

  private
  public test_quad

  integer :: count = 0
  integer :: count2 = 0

contains

  subroutine test_quad()
    real            :: I, dIda, dIdb, dIdp(2), inf, dum(0)
    type(test_case) :: tc

    call tc%init("quad")

    ! 3 * x**2
    call quad(integrand1, 1.0, 2.0, [3.0], I, dIda, dIdb, dIdp(1:1))
    call tc%assert_eq(7.0, I, 1e-14, 1e-16, "I1")
    call tc%assert_eq(-3.0, dIda, 1e-14, 1e-16, "dI1da")
    call tc%assert_eq(12.0, dIdb, 1e-14, 1e-16, "dI1db")
    call tc%assert_eq(7.0 / 3.0, dIdp(1), 1e-14, 1e-16, "dI1dp")

    ! singularity close to 1
    call quad(integrand2, -1.0, 1.0, [2.0, 1.015625], I, dIda, dIdb, dIdp)
    call tc%assert_eq(-9.7196248087233442,    I, 1e-14, 1e-16, "I2")
    call tc%assert_eq( 0.9922480620155039, dIda, 1e-14, 1e-16, "dI2da")
    call tc%assert_eq(-128.0,              dIdb, 1e-14, 1e-16, "dI2db")
    call tc%assert_eq([-4.8598124043616721, 127.007751937984496], dIdp, 1e-14, 1e-16, "dI2dp")

    ! singularity at 0
    call quad(integrand3, 0.00390625, 10.0, dum, I, dIda, dIdb, dum)
    call tc%assert_eq(5.5470845327363952, I, 1e-14, 1e-16, "I3")
    call tc%assert_eq(-255.50032552075055, dIda, 1e-14, 1e-16, "dI3da")
    call tc%assert_eq(0.000045401991009687768, dIdb, 1e-14, 1e-16, "dI3db")

    ! ! lower bound > upper bound
    call quad(integrand3, 10.0, 0.00390625, dum, I, dIda, dIdb, dum)
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

    m4_divert(m4_ifdef({m4_mpfr},0,-1))
    call test_quad_mpfr(tc)
    m4_divert(0)

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

  m4_divert(m4_ifdef({m4_mpfr},0,-1))

  subroutine test_quad_mpfr(tc)
    type(test_case), intent(inout) :: tc

    real         :: inf
    type(mpfr)   :: I
    type(string) :: e, s

    call mpfr_startup(prec = 256)
    call I%init()

    call quad_mpfr(integrand1_mpfr, 1.0, 2.0, I)
    call tc%assert_eq(7.0 / 3.0, I%to_real(), 1e-15, 1e-16, "I1 MPFR")

    call quad_mpfr(integrand2_mpfr, -1.0, 1.0, I)
    call tc%assert_eq(-4.8598124043616721, I%to_real(), 1e-15, 1e-16, "I2 MPFR")

    call quad_mpfr(integrand3_mpfr, 0.00390625, 10.0, I)
    s%s = I%to_string(n = 31)
    e%s = "5.547084532736395225915596138204E+0"
    call tc%assert_eq(e, s, "I3 MPFR")

    call quad_mpfr(integrand3_mpfr, 10.0, 0.00390625, I)
    s%s = I%to_string(n = 31)
    e%s = "-5.547084532736395225915596138204E+0"
    call tc%assert_eq(e, s, "I3a MPFR")

    ! int_{0}^{inf} exp(-x) dx
    inf = ieee_value(0.0, ieee_positive_inf)
    call quad_mpfr(integrand4_mpfr, 0.0, inf, I)
    s%s = I%to_string(n = 31)
    e%s = "1.000000000000000000000000000000E+0"
    call tc%assert_eq(e, s, "I4 MPFR")

    ! int_{-inf}^{inf} exp(-x^2) dx
    call quad_mpfr(integrand5_mpfr, -inf, inf, I)
    s%s = I%to_string(n = 31)
    e%s = "1.772453850905516027298167483341E+0"
    call tc%assert_eq(e, s, "I5 MPFR")
  end subroutine

  subroutine integrand1_mpfr(x, f)
    !! f = x**2
    type(mpfr), intent(in)    :: x
    type(mpfr), intent(inout) :: f

    call sqr(f, x)
  end subroutine

  subroutine integrand2_mpfr(x, f)
    !! f = 1 / (x - 1.015625)
    type(mpfr), intent(in)    :: x
    type(mpfr), intent(inout) :: f

    call sub(f, x, 1.015625) ! f = x - 1.015625
    call div(f, 1, f)        ! f = 1 / (x - 1.015625)
  end subroutine

  subroutine integrand3_mpfr(x, f)
    !! f = 1 / expm1(x)
    type(mpfr), intent(in)    :: x
    type(mpfr), intent(inout) :: f

    call expm1(f, x)  ! f = expm1(x)
    call div(f, 1, f) ! f = 1 / expm1(x)
  end subroutine

  subroutine integrand4_mpfr(x, f)
    !! f = exp(-x)
    type(mpfr), intent(in)    :: x
    type(mpfr), intent(inout) :: f

    call neg(f, x) ! f = -x
    call exp(f, f) ! f = exp(-x)
  end subroutine

  subroutine integrand5_mpfr(x, f)
    !! f = exp(-x**2)
    type(mpfr), intent(in)    :: x
    type(mpfr), intent(inout) :: f

    call sqr(f, x) ! f = x**2
    call neg(f, f) ! f = -x**2
    call exp(f, f) ! f = exp(-x**2)
  end subroutine

  m4_divert(0)

end module
