submodule (test_matrix_m) test_dense_m

  implicit none

  real, parameter :: rtol = 1e-14
  real, parameter :: atol = 1e-16

contains

  module subroutine test_dense()
    type(test_case) :: tc

    call tc%init("dense")

    ! scale real
    block
      real             :: d_0(5,5)
      real, parameter  :: l = 3.0
      type(dense_real) :: d

      d_0 = reshape([         &
        & 11, 16,  0,  0,  0, &
        &  7, 12, 17,  0,  0, &
        &  3,  8, 13, 18,  0, &
        &  0,  4,  9, 14, 19, &
        &  0,  0,  5, 10, 15  ], [5, 5])

      call d%init(d_0)
      call d%scale(l)

      call tc%assert_eq(d_0 * l, d%d, 0.0, 0.0, "scale")
    end block

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
      call tc%assert_eq(e1, v1, rtol, atol, "mul_vec: test 1")

      call d%mul_vec(v0, v1, fact_y = -0.5)
      e1 = [-17.7, -35.4, -53.1, -70.8, -88.5]
      call tc%assert_eq(e1, v1, rtol, atol, "mul_vec: test 2")

      call d%mul_vec(v1, v0, fact_y = 1.6, trans = 'T')
      e0 = [-973.34, -1952.60, -2935.70]
      call tc%assert_eq(e0, v0, rtol, atol, "mul_vec: test 3")
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
      call tc%assert_eq(e1, m1, rtol, atol, "mul_mat: test 1")

      call d%mul_mat(m0, m1, fact_y = 0.1)
      e1 = reshape([ 25.3,   50.6,   75.9,   3.3,    6.6,    9.9,        &
                    -18.7,  -37.4,  -56.1, -40.7,  -81.4, -122.1,        &
                    -62.7, -125.4, -188.1, -84.7, -169.4, -254.1], [3, 6])
      call tc%assert_eq(e1, m1, rtol, atol, "mul_mat: test 2")

      call d%mul_mat(m1, m0, fact_y = -1.2, trans = 'T')
      e0 = reshape([  353.84,  706.84, 1059.84, 1412.84,   48.24,   93.24,        &
                      138.24,  183.24, -257.36, -520.36, -783.36,-1046.36,        &
                     -562.96,-1133.96,-1704.96,-2275.96, -868.56,-1747.56,        &
                    -2626.56,-3505.56,-1174.16,-2361.16,-3548.16,-4735.16], [4, 6])
      call tc%assert_eq(e0, m0, rtol, atol, "mul_mat: test 3")
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
      call tc%assert_eq(e(:,1), x(:,1), rtol, atol, "solve_vec")
      call d%solve_mat(b, x)
      call tc%assert_eq(e, x, rtol, atol, "solve_mat")
    end block

    ! transpose
    block
      integer           :: i, j
      type(dense_real)  :: d, d2
      real              :: mat(3,3)

      call d%init(3)
      do i = 1, 3
        do j = 1, 3
          d%d(i,j) = i-2*j
          mat(i,j) = i-2*j
        end do
      end do

      call d%transpose
      mat = transpose(mat)
      call tc%assert_eq(d%d, mat, rtol, atol, "transpose")

      call d%transpose(d2)
      mat = transpose(mat)
      call tc%assert_eq(d2%d, mat, rtol, atol, "transpose")
    end block

    ! eig real
    block
      complex          :: e(2), e_exp(2), L(2,2), L_exp(2,2), R(2,2), R_exp(2,2)
      type(dense_real) :: d

      call d%init(real(reshape([&
        &  4, 12, &
        & 12, 11  ], [2, 2], order = [2, 1])))

      e_exp = [-5, 20]
      R_exp = reshape([&
        & -0.8,  0.6,  &
        & -0.6, -0.8   ], [2, 2])
      L_exp = R_exp

      call d%eig(e, R = R, L = L, sort = .true.)

      call tc%assert_eq(e_exp, e, rtol, atol, "eig real eval")
      call tc%assert_eq(R_exp, R, rtol, atol, "eig real R")
      call tc%assert_eq(L_exp, L, rtol, atol, "eig real L")
    end block

    ! eig cmplx
    block
      complex           :: e(2), e_exp(2), L(2,2), L_exp(2,2), R(2,2), R_exp(2,2)
      real, parameter   :: sqrt2 = 1/sqrt(2.0)
      type(dense_cmplx) :: d

      e_exp = [(1, 1), (1, -1)]
      call example_matrix6(d)
      call d%eig(e, R = R, L = L, sort = .true.)

      R_exp = reshape([      &
        & (sqrt2), ( sqrt2), &
        & (sqrt2), (-sqrt2)  ], [2, 2], order = [2, 1])
      L_exp = reshape([      &
        & (sqrt2), (-sqrt2), &
        & (sqrt2), ( sqrt2)  ], [2, 2], order = [2, 1])

      call tc%assert_eq(e_exp, e, rtol, atol, "eig cmplx eval")
      call tc%assert_eq(R_exp, R, rtol, atol, "eig cmplx R")
      call tc%assert_eq(L_exp, L, rtol, atol, "eig cmplx L")
    end block

    ! eye real
    block
      real             :: d_eye(5,5)
      type(dense_real) :: eye

      d_eye = reshape([   &
        &  1, 0, 0, 0, 0, &
        &  0, 1, 0, 0, 0, &
        &  0, 0, 1, 0, 0, &
        &  0, 0, 0, 1, 0, &
        &  0, 0, 0, 0, 1  ], [5, 5], order = [2, 1])

      eye = dense_eye_real(5)

      call tc%assert_eq(d_eye, eye%d, 0.0, 0.0, "eye")
    end block

    call tc%finish()
  end subroutine
end submodule
