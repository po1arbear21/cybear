module test_dual_m
  use test_case_m
  use dual_m
  implicit none

contains

  subroutine test_dual()
    type(test_case) :: tc
    real, parameter :: tol = 1e-14

    print "(A)", "test_dual"
    call tc%init("dual")

    ! test addition
    block
      type(dual) :: x, y

      call x%init(1, 1.0, 1)

      y = 3.0 + ((x + x) + 2.0)

      call tc%assert_eq(7.0, y%x    , tol, "addition value")
      call tc%assert_eq(2.0, y%dx(1), tol, "addition derivative")
    end block

    ! test subtraction
    block
      type(dual) :: x, y

      call x%init(1, 1.0, 1)

      y = 5.0 - (((-x) - x) - 1.0)

      call tc%assert_eq(8.0, y%x    , tol, "subtraction value")
      call tc%assert_eq(2.0, y%dx(1), tol, "subtraction derivative")
    end block

    ! test multiplication
    block
      type(dual) :: x, y

      call x%init(1, 2.0, 1)

      y = 4.0 * ((x * x) * 0.5)

      call tc%assert_eq(8.0, y%x    , tol, "multiplication value")
      call tc%assert_eq(8.0, y%dx(1), tol, "multiplication derivative")
    end block

    ! test division
    block
      type(dual) :: x, y, z

      call x%init(2, 4.0, 1)
      call y%init(2, 2.0, 2)

      z = 1.0 / ((x / y) / 0.5)

      call tc%assert_eq( 0.2500, z%x    , tol, "division value")
      call tc%assert_eq(-0.0625, z%dx(1), tol, "division derivative 1")
      call tc%assert_eq( 0.1250, z%dx(2), tol, "division derivative 2")
    end block

    ! test power
    block
      type(dual) :: x, y, z

      call x%init(2, 1.5, 1)
      call y%init(2, 2.0, 2)


      z = (((2.0 ** x) ** y) ** 0.25) ** 4

      call tc%assert_eq( 8.0           , z%x    , tol, "power value")
      call tc%assert_eq(16.0 * log(2.0), z%dx(1), tol, "power derivative 1")
      call tc%assert_eq(12.0 * log(2.0), z%dx(2), tol, "power derivative 2")
    end block

    ! test abs
    block
      type(dual) :: x, y, z

      call x%init(1, 5.0, 1)
      call y%init(1, -5.0, 1)

      z = abs(x)
      call tc%assert_eq(5.0, z%x    , tol, "abs value 1")
      call tc%assert_eq(1.0, z%dx(1), tol, "abs derivative 1")

      z = abs(y)
      call tc%assert_eq( 5.0, z%x    , tol, "abs value 2")
      call tc%assert_eq(-1.0, z%dx(1), tol, "abs derivative 2")
    end block

    ! test sqrt
    block
      type(dual) :: x, y

      call x%init(1, 9.0, 1)

      y = sqrt(x)

      call tc%assert_eq(3.0    , y%x    , tol, "sqrt value")
      call tc%assert_eq(1.0/6.0, y%dx(1), tol, "sqrt derivative")
    end block

    ! test exp
    block
      type(dual) :: x, y

      call x%init(1, 2.0, 1)

      y = exp(x)

      call tc%assert_eq(7.38905609893065, y%x    , tol, "exp value")
      call tc%assert_eq(7.38905609893065, y%dx(1), tol, "exp derivative")
    end block

    ! test log
    block
      type(dual) :: x, y

      call x%init(1, 2.0, 1)

      y = log(x)

      call tc%assert_eq(0.6931471805599453, y%x    , tol, "log value")
      call tc%assert_eq(0.5               , y%dx(1), tol, "log derivative")
    end block

    ! test sin
    block
      type(dual) :: x, y

      call x%init(1, 3.0, 1)

      y = sin(x)

      call tc%assert_eq( 0.1411200080598672, y%x    , tol, "sin value")
      call tc%assert_eq(-0.9899924966004455, y%dx(1), tol, "sin derivative")
    end block

    ! test cos
    block
      type(dual) :: x, y

      call x%init(1, 3.0, 1)

      y = cos(x)

      call tc%assert_eq(-0.9899924966004455, y%x    , tol, "cos value")
      call tc%assert_eq(-0.1411200080598672, y%dx(1), tol, "cos derivative")
    end block

    ! test tan
    block
      type(dual) :: x, y

      call x%init(1, 3.0, 1)

      y = tan(x)

      call tc%assert_eq(-0.1425465430742778, y%x    , tol, "tan value")
      call tc%assert_eq( 1.0203195169424269, y%dx(1), tol, "tan derivative")
    end block

    call tc%finish()
  end subroutine

end module
