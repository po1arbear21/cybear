submodule(test_matrix_m) test_triang_m
  use matrix_m
  implicit none

contains

  module subroutine test_triang
    type(test_case)   :: tc
    type(triang_real) :: d_lower, d_upper

    call tc%init("triang")

    call d_upper%init(3, .true.)
    d_upper%d = reshape([1,2,3,0,4,5,0,0,6], [3, 3], order=[2, 1])

    call d_lower%init(3, .false.)
    d_lower%d = reshape([1,0,0,2,3,0,4,5,6], [3, 3], order=[2, 1])

    ! mul_vec
    block
      real :: x(3), y(3), y_exp(3)

      ! test for upper triang
      x = [-1,2,5]
      y = [2,5,7]
      y_exp = [25.000, 50.500, 54.500]
      call d_upper%mul_vec(x, y, fact_y=3.5)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec: upper")

      ! test for lower triang
      x = [-1,2,5]
      y = [2,5,7]
      y_exp = [6.000, 21.500, 60.500]
      call d_lower%mul_vec(x, y, fact_y=3.5)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec: lower")
    end block

    ! mul_mat
    block
      real :: x(3,3), y(3,3), y_exp(3,3)

      ! test mul_mat: y <- this * x + fact_y * y
      x = reshape([8,3,4,1,5,9,6,7,2], [3,3])                     ! magic(3)
      y = reshape([8,3,4,1,5,9,6,7,2], [3,3], order=[2,1])        ! magic(3)'

      call d_upper%mul_mat(x, y, fact_y=3.0)
      y_exp = reshape([50,35,42,47,80,75,38,65,18], [3,3])        ! reshape(M_up *magic(3) + 3*magic(3)',[1,9])
      call tc%assert_eq(y_exp, y, 1e-12, "mul_mat: add")

      ! test mul_mat2: y <- fact * op(x) * y
      y = reshape([8,3,4,1,5,9,6,7,2], [3,3], order=[2,1])        ! magic(3)'
      call d_lower%mul_mat2(y, fact=5.0)
      y_exp = reshape([40,95,365,15,105,395,20,175,365], [3,3])   ! reshape(5 * M_low * magic(3)', [1,9])
      call tc%assert_eq(y_exp, y, 1e-12, "mul_mat: sides")
    end block

    ! solve_vec
    block
      real :: b(3), x(3), x_exp(3)

      ! solvevec
      b     = [-1,2,5]
      x_exp = [-2.416666666666667e+00, -5.416666666666665e-01, 8.333333333333333e-01]

      call d_upper%solve_vec(b, x)
      call tc%assert_eq(x_exp, x, 1e-12, "solve_vec")
    end block

    ! solve_mat
    block
      real :: b(3,3), x(3,3), x_exp(3,3)

      ! solvemat
      b = reshape([8,3,4,1,5,9,6,7,2], [3,3])                     ! magic(3)
      x_exp = reshape([ 8.0,-4.333333333333333e+00,-1.055555555555556e+00,1.0,1.0,0.0,6.0, -1.666666666666667e+00, &
        &               -2.277777777777778e+00], [3,3])

      call d_lower%solve_mat(b, x)
      call tc%assert_eq(x_exp, x, 1e-12, "solve_mat")
    end block

    call tc%finish()
  end subroutine

end submodule
