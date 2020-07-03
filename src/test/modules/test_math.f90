module test_math_m
  use test_case_m
  use math_m
  implicit none

contains

  subroutine test_math
    type(test_case) :: tc

    print "(A)", "test_math"
    call tc%init("math")

    ! cross product
    block
      integer :: i, j
      real    :: a(3), b(3), v(3), v_exp(3)

      do i = 1, 3
        a    = 0.0
        a(i) = real(i)

        do j = 1, 3
          b    = 0.0
          b(j) = real(j)

          v = cross_product(a, b)

          v_exp = 0.0
          if (j == modulo(i, 3)+1) then
            v_exp(modulo(i-2, 3)+1) = real(i*j)
          else if (j == modulo(i-2, 3)+1) then
            v_exp(modulo(i, 3)+1) = -real(i*j)
          end if

          call tc%assert_eq(v_exp, v, 1e-15, "cross product")
        end do
      end do
    end block

    ! cross product 2d
    block
      integer :: i, j
      real    :: a(3), b(3), v, v3d(3), v_exp

      do i = 1, 3
        a    = 0.0
        a(i) = real(i)

        do j = 1, 3
          b    = 0.0
          b(j) = real(j)

          v = cross_product_2d(a(1:2), b(1:2))
          v3d = cross_product(a, b)
          v_exp = v3d(3)

          call tc%assert_eq(v_exp, v, 1e-15, "cross product 2d")
        end do
      end do
    end block

    call tc%finish
  end subroutine

end module
