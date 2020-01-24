module test_high_precision_m
  use test_case_m
  use high_precision_m

  implicit none

contains

  subroutine test_high_precision()
    type(test_case) :: tc

    print "(1A)", "test_high_precision"
    call tc%init("high_precision")

    ! testing sum
    block
      real, allocatable :: p(:)
      integer :: K
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

    ! testing dot
    block
      real, allocatable :: x(:), y(:)
      integer :: K
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

    call tc%finish()
  end subroutine

end module
