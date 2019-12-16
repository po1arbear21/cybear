#include "../../util/assert.f90.inc"

submodule(test_matrix_m) test_dense_m
  use matrix_m
  implicit none

contains
  module subroutine test_dense
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

    ! to_dense
    block
      integer           :: i, j
      type(dense_real)  :: d1, d2
      real              :: val

      call d1%init(3)
      call d2%init(5)
      do i = 1, 5
        do j = 1, 5
          if (i <= 3 .and. j <= 3) d1%d(i,j) = 1.0
          d2%d(i,j) = 2.0
        end do
      end do

      call d1%to_dense(d2, i0=1, j0=3)

      do i = 1, 5
        do j = 1, 5
          val = 2.0
          if (i <= 3 .and. j >= 3) val = 1.0
          call tc%assert_eq(d2%d(i,j), val, 1e-12, "to_dense")
        end do
      end do
    end block

    ! to_sparse
    block
      integer           :: i, j
      type(dense_real)  :: d
      real              :: val, mat(5,5), e(5), v(5), v_exp(5)
      real, allocatable :: tmp(:,:)
      logical, allocatable :: sparsity(:,:)
      type(sparse_real)  :: s
      type(spbuild_real) :: sb

      ! test 1: add dense to empty sparse
      call s%init(5)
      call sb%init(s)

      call d%init(3)
      mat = 0.0
      do i = 1, 3
        do j = 1, 3
          d%d(i,j)    = i-2*j
          mat(i, j+2) = i-2*j
        end do
      end do

      call d%to_sparse(sb, i0=1, j0=3)
      call sb%save

      ! test 2: add dense to partially filled sparse
      do i = 1, 5
        e     = 0.0
        e(i)  = 1.0
        v_exp = matmul(mat, e)
        call s%mul_vec(e, v)
        call tc%assert_eq(v, v_exp, 1e-12, "to_sparse")
      end do

      tmp = reshape([1,3,5,7],[2,2],order=[2,1])
      call d%init(tmp)

      call d%to_sparse(sb, i0=4, j0=2)
      call s%init(5)
      call sb%save
      mat(4:5,2:3) = tmp

      do i = 1, 5
        e     = 0.0
        e(i)  = 1.0
        v_exp = matmul(mat, e)
        call s%mul_vec(e, v)
        call tc%assert_eq(v, v_exp, 1e-12, "to_sparse")
      end do

      ! test 3: use sparsity struct
      tmp = reshape([4,1,0,2, &
                    &8,5,2,0, &
                    &0,9,6,3, &
                    &1,0,10,7 ], [4, 4], order=[2, 1])
      call d%init(tmp)

      allocate(sparsity(4,4), source=.true.)
      sparsity(1,2) = .false.
      sparsity(2,2) = .false.
      sparsity(4,3) = .false.

      mat = 0.0
      mat(1:4,2:5) = tmp
      mat(1,2+1) = 0.0
      mat(2,2+1) = 0.0
      mat(4,3+1) = 0.0

      call s%init(5)
      call sb%init(s)
      call d%to_sparse(sb, i0=1, j0=2, struct=sparsity)
      call sb%save

      do i = 1, 5
        e     = 0.0
        e(i)  = 1.0
        v_exp = matmul(mat, e)
        call s%mul_vec(e, v)
        call tc%assert_eq(v, v_exp, 1e-12, "to_sparse: sparsity")
      end do

      ! test 4: use drop_zeros
      mat = 0.0
      mat(2:5,1:4) = tmp

      call s%init(5)
      call sb%init(s)
      call d%to_sparse(sb, i0=2, j0=1, drop_zeros=.true.)
      call sb%save

      do i = 1, 5
        e     = 0.0
        e(i)  = 1.0
        v_exp = matmul(mat, e)
        call s%mul_vec(e, v)
        call tc%assert_eq(v, v_exp, 1e-12, "to_sparse: drop zeros")
      end do
      call tc%assert_eq(count(s%a == 0.0), 0, "to_sparse: drop zeros")
    end block

    ! transpose
    block
      integer           :: i, j
      type(dense_real)  :: d, d2
      real              :: val, mat(3,3), e(3), v(3), v_exp(3)

      call d%init(3)
      do i = 1, 3
        do j = 1, 3
          d%d(i,j) = i-2*j
          mat(i,j) = i-2*j
        end do
      end do

      call d%transpose
      mat = transpose(mat)
      call tc%assert_eq(d%d, mat, 1e-12, "transpose")

      call d%transpose(d2)
      mat = transpose(mat)
      call tc%assert_eq(d2%d, mat, 1e-12, "transpose")
    end block

    ! add dense
    block
      integer           :: i, j
      type(dense_real)  :: d1, d2, d3
      real              :: val, mat1(3,3), mat2(3,3), mat3(3,3), e(3), v(3), v_exp(3)

      call d1%init(3)
      do i = 1, 3
        do j = 1, 3
          d1%d(i,j) = i-2*j
          mat1(i,j) = i-2*j
        end do
      end do

      call d2%init(3)
      do i = 1, 3
        do j = 1, 3
          d2%d(i,j) = 5*j
          mat2(i,j) = 5*j
        end do
      end do

      call d1%add_dense(d2, fact=2.0)
      mat1 = mat1 + 2*mat2
      call tc%assert_eq(d1%d, mat1, 1e-12, "add dense")

      call d1%add_dense(d2, d3, fact1=3.0, fact2=2.0)
      mat3 = 3*mat1 + 2*mat2
      call tc%assert_eq(d3%d, mat3, 1e-12, "add dense3")
    end block

    ! add sparse
    block
      integer           :: i, j
      type(dense_real)  :: d1, d2, d3
      real              :: val, mat1(3,3), mat2(3,3), mat3(3,3)
      type(sparse_real)  :: s1
      type(spbuild_real) :: sb1

      ! build sparse s1, mat1
      call s1%init(3)
      call sb1%init(s1)

      call d1%init(3)
      do i = 1, 3
        do j = 1, 3
          d1%d(i,j) = i-2*j
        end do
      end do
      mat1 = d1%d

      call d1%to_sparse(sb1)
      call sb1%save

      ! build dense d2, mat2
      call d2%init(3)
      do i = 1, 3
        do j = 1, 3
          d2%d(i,j) = i*i-3*j
        end do
      end do
      mat2 = d2%d

      ! test adding
      call d2%add_sparse(s1, fact=3.0)
      mat2 = mat2 + 3*mat1
      call tc%assert_eq(d2%d, mat2, 1e-12, "add sparse")

      call d3%init(3)
      call d2%add_sparse(s1, d3, fact1=3.0, fact2=2.0)
      mat3 = 3*mat2 + 2*mat1
      call tc%assert_eq(d3%d, mat3, 1e-12, "add sparse3")
    end block

    ! add band
    block
      integer           :: i, j
      real              :: mat1(4,4), mat2(4,4), diags(3,4), mat_b(4,4)
      type(band_real)   :: b
      type(dense_real)  :: d1, d2

      ! build band
      diags = reshape([0,1,2,3, &
                      &4,5,6,7, &
                      &8,9,10,0 ], [3,4], order=[2,1])
      mat_b = reshape([4,1,0,0, &
                      &8,5,2,0, &
                      &0,9,6,3, &
                      &0,0,10,7 ], [4, 4], order=[2, 1])
      call b%init(4, 1, d0=diags)

      ! dense d1
      call d1%init(4)
      do i = 1, 4
        do j = 1, 4
          d1%d(i,j) = i-2*j
        end do
      end do
      mat1 = d1%d

      ! add band
      call d1%add_band(b, fact=2.0)
      mat1 = mat1 + 2*mat_b
      call tc%assert_eq(d1%d, mat1, 1e-12, "add band")

      call d1%add_band(b, d2, fact1=3.0, fact2=2.0)
      mat2 = 3*mat1 + 2*mat_b
      call tc%assert_eq(d2%d, mat2, 1e-12, "add band3")
    end block

    ! mul dense
    block
      integer :: i, j
      real    :: mat1(3,4), mat2(4,2), mat3(3,2)
      type(dense_real) :: d1, d2, d3

      do i = 1, 4
        do j = 1, 4
          if (i < 4) mat1(i,j) = i-2*j
          if (j < 3) mat2(i,j) = i*i + j
        end do
      end do
      call d1%init(mat1)
      call d2%init(mat2)

      call d1%mul_dense(d2, d3)
      mat3 = matmul(mat1, mat2)
      call tc%assert_eq(d3%d, mat3, 1e-12, "mul dense")
    end block

    ! eig
    block
      type(dense_real)  :: d
      real :: mat(3,3)
      complex, dimension(3,3) :: R, R_exp, L, L_exp
      complex, dimension(3) :: e, e_exp
      integer :: i

      mat = 0.0
      mat(1,1) = 1.0
      mat(2,2) = 2.0
      mat(3,3) = 3.0

      e_exp = [1,2,3]
      d = eye_real(3)
      R_exp = d%d
      L_exp = d%d

      call d%init(mat)
      call d%eig(e, R=R, L=L, sort=.true.)

      ! scale Evecs
      do i = 1, 3
        R(:,i) = R_exp(i,i) / R(i,i) * R(:,i)
        L(i,:) = L_exp(i,i) / L(i,i) * L(i,:)
      end do

      call tc%assert_eq(e_exp, e, 1e-12, "eig eval")
      call tc%assert_eq(R_exp, R, 1e-12, "eig R")
      call tc%assert_eq(L_exp, L, 1e-12, "eig L")
    end block

    call tc%finish
  end subroutine
end submodule
