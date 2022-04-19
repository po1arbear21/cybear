m4_include(../../util/macro.f90.inc)

m4_divert(m4_ifdef({m4_quadpack},0,-1))

module test_quadpack_m

  use, intrinsic :: ieee_arithmetic
  use quadpack_m,  only: quadpack_int
  use test_case_m, only: test_case

  implicit none

  private
  public test_quadpack

contains

  subroutine test_quadpack()
    type(test_case) :: tc
    real            :: alpha, xmin, xmax, res, res_exp

    call tc%init("quadpack")

    alpha   = 0.8
    xmin    = 1.5
    xmax    = 3.0
    res_exp = alpha * (xmax**3 - xmin**3) / 3.0
    call quadpack_int(fun1, xmin, xmax, 1e-12, 1e-12, res)
    call tc%assert_eq(res_exp, res, 1e-10, "integral 1")

    alpha = 2.0
    xmin = 0.3
    xmax = ieee_value(xmax, IEEE_POSITIVE_INF)
    res_exp = exp(- 0.3 * alpha) / alpha
    call quadpack_int(fun2, xmin, xmax, 1e-12, 1e-12, res)
    call tc%assert_eq(res_exp, res, 1e-10, "integral 2")

    call tc%finish()

  contains

    function fun1(x) result(f)
      real, intent(in) :: x
      real             :: f

      f = alpha * x**2
    end function

    function fun2(x) result(f)
      real, intent(in) :: x
      real             :: f

      f = exp(-alpha * x)
    end function

  end subroutine

end module

m4_divert(0)
