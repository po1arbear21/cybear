module test_polynomial_m
  use test_case_m
  use math_m
  use polynomial_m
  implicit none

contains

  subroutine test_polynomial()
    type(test_case) :: tc

    print "(1A)", "test_polynomial"
    call tc%init("polynomial")

    ! 1D quadratic
    block
      integer, parameter :: N = 20
      integer            :: i
      real               :: x1(1,3), f1(3), x2(1,N), f2(N), f(N)
      real,    parameter :: tol = 1e-13
      type(polynomial)   :: poly

      ! f(x) = x^2 + 2x + 3
      x1(1,1) = -2.0
      x1(1,2) =  1.0
      x1(1,3) =  1.7
      f1 = x1(1,:)**2 + 2.0*x1(1,:) + 3.0

      x2(1,:) = linspace(-5.0, 5.0, N)
      f2 = x2(1,:)**2 + 2.0*x2(1,:) + 3.0

      call poly%init(1, [2])
      call poly%interpolate(x1, f1)

      do i = 1, N
        call poly%eval(x2(:,i), f(i))
      end do
      call tc%assert_eq(f2, f, tol, "1D quadratic")
    end block

    ! 2D bilinear
    block
      integer, parameter :: N = 20
      integer            :: i
      real               :: r1(2,4), f1(4), r2(2,N*N), f2(N*N), f(N*N)
      real,    parameter :: tol = 1e-13
      type(polynomial)   :: poly

      ! unit square
      r1(:,1) = [ 0.0, 0.0 ]
      r1(:,2) = [ 1.0, 0.0 ]
      r1(:,3) = [ 0.0, 1.0 ]
      r1(:,4) = [ 1.0, 1.0 ]
      f1(1)   = 0.7
      f1(2)   = 1.5
      f1(3)   = 9.2
      f1(4)   = 3.8

      do i = 1, N
        r2(1,1+(i-1)*N:i*N) = linspace(r1(1,1), r1(1,2), N)
        r2(2,i:i+N-1:N) = linspace(r1(2,1), r1(2,3), N)
      end do

      do i = 1, N*N
        f2(i) = f1(1) * (1.0 - r2(1,i))*(1.0 - r2(2,i)) &
              + f1(2) *        r2(1,i) *(1.0 - r2(2,i)) &
              + f1(3) * (1.0 - r2(1,i))*       r2(2,i)  &
              + f1(4) *        r2(1,i) *       r2(2,i)
      end do

      call poly%init(2, [1, 1])
      call poly%interpolate(r1, f1)

      do i = 1, N*N
        call poly%eval(r2(:,i), f(i))
      end do
      call tc%assert_eq(f2, f, tol, "2D bilinear")

    end block

    call tc%finish()
  end subroutine

end module
