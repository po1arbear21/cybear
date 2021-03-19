module test_math_m
  use test_case_m
  use math_m
  implicit none

  private
  public test_math

contains

  subroutine test_math()
    type(test_case) :: tc

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

    ! check linear dependence
    block
      real    :: M(3,8)
      logical :: l

      M(1,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
      M(2,:) = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
      M(3,:) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
      l = check_lin_dep(M)
      call tc%assert(.not. l, "check linear dependence 1")

      M(1,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
      M(2,:) = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
      M(3,:) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
      l = check_lin_dep(M)
      call tc%assert(.not. l, "check linear dependence 2")

      M(1,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
      M(2,:) = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
      M(3,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
      l = check_lin_dep(M)
      call tc%assert(l, "check linear dependence 3")

      M(1,:) = [0.3371, 0.6020, 0.0838, 0.9961, 0.1622, 0.2630, 0.2290, 0.0782]
      M(2,:) = [0.7943, 0.6541, 0.9133, 0.4427, 0.3112, 0.6892, 0.1524, 0.1067]
      M(3,:) = [2.1427, 3.0621, 1.2485, 4.4271, 0.9600, 1.7412, 1.0684, 0.4195]
      l = check_lin_dep(M)
      call tc%assert(l, "check linear dependence 4")

      M(2,8) = M(2,8) + 0.1
      l = check_lin_dep(M)
      call tc%assert(.not. l, "check linear dependence 5")
    end block

    call tc%finish
  end subroutine

end module
