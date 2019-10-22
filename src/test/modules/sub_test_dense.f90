submodule(test_matrix_m) test_dense_m
  use matrix_m
  implicit none

contains
  module subroutine test_dense()
    type(test_case) :: tc
    print "(1A)", "test_dense"
    call tc%init("dense")

    ! mul_vec
    block
      integer          :: i, j
      real             :: v0(3), v1(5), e0(3), e1(5)
      type(dense_real) :: d

      call d%init(5, ncols = 3)
      do i = 1, 5
        do j = 1, 3
          d%d(i,j) = i * j
        end do
      end do
      do i = 1, 3
        v0(i) = 1.3 - 1.2 * i * i
      end do

      call d%mul_vec(v0, v1)
      e1 = [-35.4, -70.8, -106.2, -141.6, -177.0]
      call tc%assert_eq(e1, v1, 1e-12, "mul_vec: test 1")

      call d%mul_vec(v0, v1, fact_y = -0.5)
      e1 = [-17.7, -35.4, -53.1, -70.8, -88.5]
      call tc%assert_eq(e1, v1, 1e-12, "mul_vec: test 2")

      call d%mul_vec(v1, v0, fact_y = 1.6, trans = 'T')
      e0 = [-973.34, -1952.60, -2935.70]
      call tc%assert_eq(e0, v0, 1e-12, "mul_vec: test 3")
    end block

    ! mul_mat
    block
      integer          :: i, j
      real             :: m0(4,6), m1(3,6), e0(4,6), e1(3,6)
      type(dense_real) :: d

      call d%init(3, ncols = 4)
      do i = 1, 3
        do j = 1, 4
          d%d(i,j) = i * j
        end do
      end do
      do i = 1, 4
        do j = 1, 6
          m0(i,j) = i - 2 * j + 1.3
        end do
      end do

      call d%mul_mat(m0, m1)
      e1 = reshape([ 23.0,   46.0,   69.0,   3.0,    6.0,    9.0,        &
                    -17.0,  -34.0,  -51.0, -37.0,  -74.0, -111.0,        &
                    -57.0, -114.0, -171.0, -77.0, -154.0, -231.0], [3, 6])
      call tc%assert_eq(e1, m1, 1e-12, "mul_mat: test 1")

      call d%mul_mat(m0, m1, fact_y = 0.1)
      e1 = reshape([ 25.3,   50.6,   75.9,   3.3,    6.6,    9.9,        &
                    -18.7,  -37.4,  -56.1, -40.7,  -81.4, -122.1,        &
                    -62.7, -125.4, -188.1, -84.7, -169.4, -254.1], [3, 6])
      call tc%assert_eq(e1, m1, 1e-12, "mul_mat: test 2")

      call d%mul_mat(m1, m0, fact_y = -1.2, trans = 'T')
      e0 = reshape([  353.84,  706.84, 1059.84, 1412.84,   48.24,   93.24,        &
                      138.24,  183.24, -257.36, -520.36, -783.36,-1046.36,        &
                     -562.96,-1133.96,-1704.96,-2275.96, -868.56,-1747.56,        &
                    -2626.56,-3505.56,-1174.16,-2361.16,-3548.16,-4735.16], [4, 6])
      call tc%assert_eq(e0, m0, 1e-12, "mul_mat: test 3")
    end block

    ! factorize, solve_vec, solve_mat
    block
      integer          :: i, j
      real             :: b(4,3), x(4,3), e(4,3)
      type(dense_real) :: d

      call d%init(4)
      do i = 1, 4
        do j = 1, 4
          d%d(i,j) = 1.6 / i - 1.2 * i * i * j
        end do
        d%d(i,i) = d%d(i,i) + 1.0
        do j = 1, 3
          b(i,j) = 6.2 - (2.0 * i) / (j + 31.6)
        end do
      end do

      call d%factorize()

      e = reshape([1.6971935705406791e+00, 2.2974350695435128e+00, 6.7474471559634475e-01,-2.1718802483809423e+00,       &
                   1.6973860158693055e+00, 2.2984958657452248e+00, 6.7528450127420747e-01,-2.1729692073136730e+00,       &
                   1.6975673371905025e+00, 2.2994953442474140e+00, 6.7579308546779970e-01,-2.1739952206433548e+00], [4,3])

      call d%solve_vec(b(:,1), x(:,1))
      call tc%assert_eq(e(:,1), x(:,1), 1e-12, "solve_vec")
      call d%solve_mat(b, x)
      call tc%assert_eq(e, x, 1e-12, "solve_mat")
    end block

    ! ! to_dense, to_sparse
    ! block
    !   integer          :: i, j
    !   type(dense_real) :: d1, d2


    ! end block

    call tc%finish()
  end subroutine
end submodule