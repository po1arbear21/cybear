submodule (test_matrix_m) test_arith_m
  implicit none

contains

  module subroutine test_arith()
    type(test_case) :: tc
    call tc%init("arith")

    ! fixme test interface `add` here. currently it is included in matrix tests: sub_test_dense, etc.
    call test_diag(tc)
    call test_approx(tc)

    call tc%finish()
  end subroutine

  subroutine test_diag(tc)
    !! tests procedures: matrix -> diagonal array and vice versa.
    !!
    !! matrix:
    !!  real: [1 2 3 4]         complex: [1+2i 3+4i 5+6i 7+8i]
    !!        [5 6 7 8]                  [9    1+3i 5+7i 9+1i]
    !!        [9 0 1 2]                  [4+7i   4i 8+2i 7+2i]
    !!
    !! diagonals
    !!  real: [1 6 1]           complex: [1+2i 1+3i 8+2i]

    type(test_case), intent(inout) :: tc

    complex :: matr_in_c(3,4), matr_out_c(3,3), exp_arr_c(3)
    real    :: matr_in_r(3,4), matr_out_r(3,3), exp_arr_r(3)

    ! setup expected in/output
    block
      integer :: i

      matr_in_r = reshape([(mod(i, 10), i=1,12)], [3, 4], order=[2, 1])
      exp_arr_r = [1, 6, 1]

      matr_in_c = reshape([(1,2), (3,4), (5,6), (7,8), &
        &                   (9,0), (1,3), (5,7), (9,1), &
        &                   (4,7), (0,4), (8,2), (7,2)  ], [3, 4], order=[2, 1])
      exp_arr_c  = [(1,2), (1,3), (8,2)]

      matr_out_r = 0
      matr_out_c = 0
      do i = 1, 3
        matr_out_r(i,i) = exp_arr_r(i)
        matr_out_c(i,i) = exp_arr_c(i)
      end do
    end block

    ! dense: matrix <-> array
    block
      complex           :: arr_c(3)
      real              :: arr_r(3)
      type(dense_cmplx) :: matr_c
      type(dense_real)  :: matr_r

      ! dense real -> array
      call matr_r%init(matr_in_r)
      call matrix_diag(matr_r, arr_r)
      call tc%assert_eq(exp_arr_r, arr_r, 1e-12, "diag: dense real -> array")

      ! dense complex -> array
      call matr_c%init(matr_in_c)
      call matrix_diag(matr_c, arr_c)
      call tc%assert_eq(exp_arr_c, arr_c, 1e-12, "diag: dense complex -> array")

      ! array -> dense real
      call matrix_diag(exp_arr_r, matr_r)
      call tc%assert_eq(matr_out_r, matr_r%d, 1e-12, "diag: array -> dense real")

      ! array -> dense complex
      call matrix_diag(exp_arr_c, matr_c)
      call tc%assert_eq(matr_out_c, matr_c%d, 1e-12, "diag: array -> dense complex")
    end block

    ! sparse: matrix <-> array
    block
      complex             :: arr_c(3)
      integer             :: i, j
      real                :: arr_r(3)
      type(sparse_cmplx)  :: s_c
      type(sparse_real)   :: s_r
      type(spbuild_cmplx) :: sb_c
      type(spbuild_real)  :: sb_r

      ! sparse real -> array
      call s_r%init(3, ncols=4)
      call sb_r%init(s_r)
      do i = 1, 3
        do j = 1, 4
          call sb_r%set(i, j, matr_in_r(i,j))
        end do
      end do
      call sb_r%save()
      call matrix_diag(s_r, arr_r)
      call tc%assert_eq(exp_arr_r, arr_r, 1e-12, "diag: sparse real -> array"   )

      ! sparse complex -> array
      call s_c%init(3, ncols=4)
      call sb_c%init(s_c)
      do i = 1, 3
        do j = 1, 4
          call sb_c%set(i, j, matr_in_c(i,j))
        end do
      end do
      call sb_c%save()
      call matrix_diag(s_c, arr_c)
      call tc%assert_eq(exp_arr_c, arr_c, 1e-12, "diag: sparse complex -> array")

      ! array -> sparse real
      call matrix_diag(exp_arr_r, s_r)
      call tc%assert_eq(exp_arr_r,    s_r%a, 1e-12, "diag: array -> sparse real: a" )
      call tc%assert_eq([(i, i=1,4)], int(s_r%ia),  "diag: array -> sparse real: ia")
      call tc%assert_eq([(i, i=1,3)], s_r%ja,       "diag: array -> sparse real: ja")

      ! array -> sparse complex
      call matrix_diag(exp_arr_c, s_c)
      call tc%assert_eq(exp_arr_c,    s_c%a, 1e-12, "diag: array -> sparse complex: a" )
      call tc%assert_eq([(i, i=1,4)], int(s_c%ia),  "diag: array -> sparse complex: ia")
      call tc%assert_eq([(i, i=1,3)], s_c%ja,       "diag: array -> sparse complex: ja")
    end block
  end subroutine

  subroutine test_approx(tc)
    !! tests procedure: approx

    type(test_case), intent(inout) :: tc

    ! test: sparse_real
    !   input matrix:                   [ 10    2  -3 4]
    !                                   [ 15  -20   7 1]
    !                                   [-19    0  10 2]
    !
    !   expected result for thres=0.25: [ 10    0  -3 4]
    !                                   [ 15  -20   7 0]
    !                                   [-19    0  10 0]
    block
      integer                   :: i
      real                      :: x(4), y(3), d_in(3,4), d_out(3,4)
      type(sparse_real), target :: s_in, s_out
      type(spbuild_real)        :: sb

      ! expected in/output
      d_in  = reshape([ 10,   2, -3, 4, &
        &               15, -20,  7, 1, &
        &              -19,   0, 10, 2  ], [3, 4], order=[2, 1])
      d_out = reshape([ 10,   0, -3, 4, &
        &               15, -20,  7, 0, &
        &              -19,   0, 10, 0  ], [3, 4], order=[2, 1])

      ! build sparse in
      call s_in%init(3, ncols=4)
      call sb%init(s_in)
      do i = 1, 3
        call sb%set_row(i, d_in(i,:))
      end do
      call sb%save()

      ! apply approximation
      call matrix_approx(0.25, s_in, s_out)

      ! check s_out == d_out by matrix multiplication with unit vectors (result should be columns of d_out)
      do i = 1, 4
        x    = 0
        x(i) = 1
        call s_out%mul_vec(x, y)
        call tc%assert_eq(d_out(:,i), y, 1e-12, "approx: sparse")
      end do
    end block

  end subroutine

end submodule
