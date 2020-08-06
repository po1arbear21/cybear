module test_high_precision_m
  use test_case_m
  use high_precision_m
  implicit none

contains

  subroutine test_high_precision()
    type(test_case) :: tc

    print "(A)", "test_high_precision"
    call tc%init("high_precision")

    ! hp addition
    block
      type(hp_real) :: h1, h2, h3, h3_exp

      h1 = TwoSum(1.234567890123453,6.234578902142152e-12)
      h2 = TwoSum(4.539873945644359,3.335434545702941e-12)
      h3 = h1 + h2

      h3_exp%x = 5.7744418357773819
      h3_exp%y = 3.3506478609411802e-16

      call tc%assert_eq(h3_exp%x, h3%x, abs(h3_exp%x)*epsilon(h3_exp%x), "hp addition principal")
      call tc%assert_eq(h3_exp%y, h3%y, abs(h3_exp%y)*epsilon(h3_exp%y), "hp addition correction")
    end block

    ! hp subtraction
    block
      type(hp_real) :: h1, h2, h3, h3_exp

      h1 = TwoSum(1.234567890123453,6.234578902142152e-12)
      h2 = TwoSum(4.539873945644359,3.335434545702941e-12)
      h3 = h1 - h2

      h3_exp%x = -3.3053060555180074
      h3_exp%y =  1.2999453800233868e-16

      call tc%assert_eq(h3_exp%x, h3%x, abs(h3_exp%x)*epsilon(h3_exp%x), "hp subtraction principal")
      call tc%assert_eq(h3_exp%y, h3%y, abs(h3_exp%y)*epsilon(h3_exp%y), "hp subtraction correction")
    end block

    ! hp multiplication
    block
      type(hp_real) :: h1, h2, h3, h3_exp

      h1 = TwoSum(1.234567890123453,6.234578902142152e-12)
      h2 = TwoSum(4.539873945644359,3.335434545702941e-12)
      h3 = h1 * h2

      h3_exp%x = 5.6047825985330135
      h3_exp%y = 5.3491707868658234e-16

      call tc%assert_eq(h3_exp%x, h3%x, abs(h3_exp%x)*epsilon(h3_exp%x), "hp multiplication principal")
      call tc%assert_eq(h3_exp%y, h3%y, abs(h3_exp%y)*epsilon(h3_exp%y), "hp multiplication correction")
    end block

    ! hp division
    block
      type(hp_real) :: h1, h2, h3, h3_exp

      h1 = TwoSum(1.234567890123453,6.234578902142152e-12)
      h2 = TwoSum(4.539873945644359,3.335434545702941e-12)
      h3 = h1 / h2

      h3_exp%x = 2.7193880378842855e-01
      h3_exp%y = 1.9298233693481757e-17

      call tc%assert_eq(h3_exp%x, h3%x, abs(h3_exp%x)*epsilon(h3_exp%x), "hp division principal")
      call tc%assert_eq(h3_exp%y, h3%y, abs(h3_exp%y)*epsilon(h3_exp%y), "hp division correction")
    end block

    ! hp square root
    block
      type(hp_real) :: h1, h2, h2_exp

      h1 = TwoSum(1.234567890123453,6.234578902142152e-12)
      h2 = sqrt(h1)

      h2_exp%x =  1.1111111061139149
      h2_exp%y = -3.1190733475658156e-17

      call tc%assert_eq(h2_exp%x, h2%x, abs(h2_exp%x)*epsilon(h2_exp%x), "hp square root principal")
      call tc%assert_eq(h2_exp%y, h2%y, abs(h2_exp%y)*epsilon(h2_exp%y), "hp square root correction")
    end block

    ! hp exponential
    block
      type(hp_real) :: h1, h2, h2_exp

      h1 = TwoSum(1.234567890123453,6.234578902142152e-12)
      h2 = exp(h1)

      h2_exp%x =  3.4368930843674224
      h2_exp%y = -8.3963109305985143e-17

      call tc%assert_eq(h2_exp%x, h2%x, abs(h2_exp%x)*epsilon(h2_exp%x), "hp exponential principal")
      call tc%assert_eq(h2_exp%y, h2%y, abs(h2_exp%y)*epsilon(h2_exp%y), "hp exponential correction")
    end block

    ! hp exponential minus 1
    block
      type(hp_real) :: h1, h2, h2_exp

      h1 = TwoSum(1.234567890123453e-6,6.234578902142152e-18)
      h2 = expm1(h1)

      h2_exp%x =  1.2345686522089389e-06
      h2_exp%y = -1.0310008648997058e-22

      call tc%assert_eq(h2_exp%x, h2%x, abs(h2_exp%x)*epsilon(h2_exp%x), "hp exponential minus 1 principal")
      call tc%assert_eq(h2_exp%y, h2%y, abs(h2_exp%y)*epsilon(h2_exp%y), "hp exponential minus 1 correction")
    end block

    ! hp_sum
    block
      real, allocatable :: p(:)
      real :: res, res_exp, tol

      tol     = epsilon(1.0)

      p       = [1e100, 1.0, -1e100]
      res_exp = 1.0
      res     = hp_sum(p)

      call tc%assert_eq(res_exp, res, tol, "sum1")

      p       = [1e200, 1e100, 1.0, -1e100, -1e200]
      res_exp = 1.0
      res     = hp_sum(p, K=3)

      call tc%assert_eq(res_exp, res, tol, "sum2")
    end block

    ! hp_dot
    block
      real, allocatable :: x(:), y(:)
      real :: res, res_exp, tol

      tol     = epsilon(1.0)

      x       = [1e100, 1.0, -1e100]
      y       = [1.0, 1.0, 1.0]
      res_exp = 1.0
      res     = hp_dot(x, y)

      call tc%assert_eq(res_exp, res, tol, "dot1")

      x       = [1.0, 1.0, -1e100]
      y       = [1e100, 1.0, 1.0]
      res_exp = 1.0
      res     = hp_dot(x, y)

      call tc%assert_eq(res_exp, res, tol, "dot2")

      x       = [1e200, 1.0, 1.0, -1e100, 1.0]
      y       = [1.0, 1e100, 1.0, 1.0, -1e200]
      res_exp = 1.0
      res     = hp_dot(x, y, K=3)

      call tc%assert_eq(res_exp, res, tol, "dot3")
    end block

    call tc%finish
  end subroutine

end module
