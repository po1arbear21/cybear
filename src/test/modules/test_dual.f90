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
      integer      :: i
      type(dual_1) :: x, y(0:2)

      call x%init(1.0, 1)

      y(0)  =  3.0       + (( x     +  x    ) +  2.0      )        ! test for single addition
      y(1:) = [3.0, 3.0] + (([x, x] + [x, x]) + [2.0, 2.0])        ! test for elemental addition

      do i = 0, 2
        call tc%assert_eq(7.0, y(i)%x    , tol, "addition value"     )
        call tc%assert_eq(2.0, y(i)%dx(1), tol, "addition derivative")
      end do
    end block

    ! test subtraction
    block
      integer      :: i
      type(dual_1) :: x, y(0:2)

      call x%init(1.0, 1)

      y(0)  =  5.0       - (((- x    ) -  x    ) -  1.0      )
      y(1:) = [5.0, 5.0] - (((-[x, x]) - [x, x]) - [1.0, 1.0])

      do i = 0, 2
        call tc%assert_eq(8.0, y(i)%x    , tol, "subtraction value"     )
        call tc%assert_eq(2.0, y(i)%dx(1), tol, "subtraction derivative")
      end do
    end block

    ! test multiplication
    block
      integer      :: i
      type(dual_1) :: x, y(0:2)

      call x%init(2.0, 1)

      y(0)  =  4.0       * ((x      *  x    ) *  0.5      )
      y(1:) = [4.0, 4.0] * (([x, x] * [x, x]) * [0.5, 0.5])

      do i = 0, 2
        call tc%assert_eq(8.0, y(i)%x    , tol, "multiplication value"     )
        call tc%assert_eq(8.0, y(i)%dx(1), tol, "multiplication derivative")
      end do
    end block

    ! test division
    block
      integer      :: i
      type(dual_2) :: x, y, z(0:2)

      call x%init(4.0, 1)
      call y%init(2.0, 2)

      z(0)  =  1.0       / (( x     /  y    ) /  0.5      )
      z(1:) = [1.0, 1.0] / (([x, x] / [y, y]) / [0.5, 0.5])

      do i = 0, 2
        call tc%assert_eq( 0.2500, z(i)%x    , tol, "division value"       )
        call tc%assert_eq(-0.0625, z(i)%dx(1), tol, "division derivative 1")
        call tc%assert_eq( 0.1250, z(i)%dx(2), tol, "division derivative 2")
      end do
    end block

    ! test power
    block
      integer      :: i
      type(dual_2) :: x, y, z(0:2)

      call x%init(1.5, 1)
      call y%init(2.0, 2)

      z(0)  = ((( 2.0       **  x    ) **  y    ) **  0.25       ) **  4
      z(1:) = ((([2.0, 2.0] ** [x, x]) ** [y, y]) ** [0.25, 0.25]) ** [4, 4]

      do i = 0, 2
        call tc%assert_eq( 8.0           , z(i)%x    , tol, "power value"       )
        call tc%assert_eq(16.0 * log(2.0), z(i)%dx(1), tol, "power derivative 1")
        call tc%assert_eq(12.0 * log(2.0), z(i)%dx(2), tol, "power derivative 2")
      end do
    end block

    ! test abs
    block
      integer      :: i
      type(dual_1) :: x, y, z(0:2)

      call x%init( 5.0, 1)
      call y%init(-5.0, 1)

      z(0)  = abs(x     )
      z(1:) = abs([x, x])
      do i = 0, 2
        call tc%assert_eq(5.0, z(i)%x    , tol, "abs value 1"     )
        call tc%assert_eq(1.0, z(i)%dx(1), tol, "abs derivative 1")
      end do

      z(0)  = abs(y)
      z(1:) = abs([y, y])
      do i = 0, 2
        call tc%assert_eq( 5.0, z(i)%x    , tol, "abs value 2"     )
        call tc%assert_eq(-1.0, z(i)%dx(1), tol, "abs derivative 2")
      end do
    end block

    ! test sqrt
    block
      integer      :: i
      type(dual_1) :: x, y(0:2)

      call x%init(9.0, 1)

      y(0)  = sqrt(x)
      y(1:) = sqrt([x, x])

      do i = 0, 2
        call tc%assert_eq(3.0    , y(i)%x    , tol, "sqrt value"     )
        call tc%assert_eq(1.0/6.0, y(i)%dx(1), tol, "sqrt derivative")
      end do
    end block

    ! test exp
    block
      integer      :: i
      type(dual_1) :: x, y(0:2)

      call x%init(2.0, 1)

      y(0)  = exp(x)
      y(1:) = exp([x, x])

      do i = 0, 2
        call tc%assert_eq(7.38905609893065, y(i)%x    , tol, "exp value"     )
        call tc%assert_eq(7.38905609893065, y(i)%dx(1), tol, "exp derivative")
      end do
    end block

    ! test log
    block
      integer      :: i
      type(dual_1) :: x, y(0:2)

      call x%init(2.0, 1)

      y(0)  = log(x)
      y(1:) = log([x, x])

      do i = 0, 2
        call tc%assert_eq(0.6931471805599453, y(i)%x    , tol, "log value"     )
        call tc%assert_eq(0.5               , y(i)%dx(1), tol, "log derivative")
      end do
    end block

    ! test sin
    block
      integer      :: i
      type(dual_1) :: x, y(0:2)

      call x%init(3.0, 1)

      y(0)  = sin(x)
      y(1:) = sin([x, x])

      do i = 0, 2
        call tc%assert_eq( 0.1411200080598672, y(i)%x    , tol, "sin value"     )
        call tc%assert_eq(-0.9899924966004455, y(i)%dx(1), tol, "sin derivative")
      end do
    end block

    ! test cos
    block
      integer      :: i
      type(dual_1) :: x, y(0:2)

      call x%init(3.0, 1)

      y(0)  = cos(x)
      y(1:) = cos([x, x])

      do i = 0, 2
        call tc%assert_eq(-0.9899924966004455, y(i)%x    , tol, "cos value"     )
        call tc%assert_eq(-0.1411200080598672, y(i)%dx(1), tol, "cos derivative")
      end do
    end block

    ! test tan
    block
      integer      :: i
      type(dual_1) :: x, y(0:2)

      call x%init(3.0, 1)

      y(0)  = tan(x)
      y(1:) = tan([x, x])

      do i = 0, 2
        call tc%assert_eq(-0.1425465430742778, y(i)%x    , tol, "tan value"     )
        call tc%assert_eq( 1.0203195169424269, y(i)%dx(1), tol, "tan derivative")
      end do
    end block

    ! test dot_product
    block
      integer                   :: i
      integer,      parameter   :: n=5
      type(dual_5), allocatable :: x(:), y(:)
      type(dual_5)              :: z, z_exp

      ! simple test for simple 1-element sum
      allocate (x(1))
      call x(1)%init(1.0, i=1)

      z = x .dot. x

      call tc%assert_eq(1.0, z%x,     tol, "dot product (1 element): value"     )
      call tc%assert_eq(2.0, z%dx(1), tol, "dot product (1 element): derivative")
      deallocate (x)

      ! more complex test. n-elements
      allocate (x(n), y(n))
      do i = 1, n
        call x(i)%init(  real(i), i)
        call y(i)%init(2*real(i), i)
      end do

      z = x .dot. y

      call z_exp%init(0.0)
      do i = 1, n
        z_exp = z_exp + x(i) * y(i)
      end do

      call tc%assert_eq(z_exp%x,  z%x,  tol, "dot product (n elements): value"     )
      call tc%assert_eq(z_exp%dx, z%dx, tol, "dot product (n elements): derivative")
    end block

    call tc%finish
  end subroutine

end module
