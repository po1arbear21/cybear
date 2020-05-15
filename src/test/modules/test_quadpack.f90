module test_quadpack_m
  use test_case_m
  use quadpack_m
  implicit none

contains

  subroutine test_quadpack()
    type(test_case) :: tc
    real            :: alpha, xmin, xmax, res, res_exp

    print "(1A)", "test_quadpack"
    call tc%init("quadpack")

    alpha   = 0.8
    xmin    = 1.5
    xmax    = 3.0
    res_exp = alpha * (xmax**3 - xmin**3) / 3.0
    call qags(fun, xmin, xmax, 1e-12, 1e-12, res)

    call tc%assert_eq(res_exp, res, 1e-10, "integral")

    call tc%finish()

  contains

    function fun(x) result(f)
      real, intent(in) :: x
      real             :: f

      f = alpha * x**2
    end function

  end subroutine

end module
