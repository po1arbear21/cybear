submodule (test_matrix_m) test_arith_m

  implicit none

  real, parameter :: rtol = 1e-14
  real, parameter :: atol = 1e-16

contains

  module subroutine test_arith()
    type(test_case) :: tc
    call tc%init("arith")

    call test_add(tc)
    call test_mul(tc)
    call test_diag(tc)
    call test_approx(tc)

    call tc%finish()
  end subroutine

  subroutine test_add(tc)
    !! test matrix_add
    type(test_case), intent(inout) :: tc

    ! add dense-dense (store result in second or third matrix)
    block
      integer           :: i, j
      type(dense_real)  :: d1, d2, d3
      real              :: mat1(3,3), mat2(3,3), mat3(3,3)

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

      ! d1 <- d1 + fact*d2
      call matrix_add(d2, d1, fact=2.0)
      mat1 = mat1 + 2*mat2
      call tc%assert_eq(d1%d, mat1, rtol, atol, "add dense")

      ! d3 <- fact1*d1 + fact2*d2
      call matrix_add(d1, d2, d3, fact1=3.0, fact2=2.0)
      mat3 = 3*mat1 + 2*mat2
      call tc%assert_eq(d3%d, mat3, rtol, atol, "add dense3")
    end block

    ! add dense-sparse (store result in second or third matrix)
    block
      integer            :: i, j
      type(dense_real)   :: d1, d2, d3
      real               :: mat1(3,3), mat2(3,3), mat3(3,3)
      type(sparse_real)  :: s1
      type(spbuild_real) :: sb1

      ! init matrices
      call s1%init(3)
      call d1%init(3)
      call d3%init(3)

      ! build d1
      do i = 1, 3
        do j = 1, 3
          d1%d(i,j) = i-2*j
        end do
      end do
      mat1 = d1%d

      ! test adding empty sparse matrix to d1
      call matrix_add(s1, d1)
      call tc%assert_eq(d1%d, mat1, rtol, atol, "add empty sparse")
      call matrix_add(s1, d1, d3, fact1=7.8, fact2=2.5)
      call tc%assert_eq(d3%d, 2.5 * mat1, rtol, atol, "add empty sparse3")

      ! build sparse s1
      call sb1%init(s1)
      call matrix_convert(d1, sb1)
      call sb1%save()

      ! build dense d2, mat2
      call d2%init(3)
      do i = 1, 3
        do j = 1, 3
          d2%d(i,j) = i*i-3*j
        end do
      end do
      mat2 = d2%d

      ! test adding: d2 <- d2 + fact*s1
      call matrix_add(s1, d2, fact=3.0)
      mat2 = mat2 + 3*mat1
      call tc%assert_eq(d2%d, mat2, rtol, atol, "add sparse")

      ! d3 <- fact1*s1 + fact2*d2
      call matrix_add(s1, d2, d3, fact1=2.0, fact2=3.0)
      mat3 = 2*mat1 + 3*mat2
      call tc%assert_eq(d3%d, mat3, rtol, atol, "add sparse3")
    end block

    ! add dense-band (store result in second or third matrix)
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
      call matrix_add(b, d1, fact=2.0)
      mat1 = mat1 + 2*mat_b
      call tc%assert_eq(d1%d, mat1, rtol, atol, "add band")

      call matrix_add(b, d1, d2, fact1=2.0, fact2=3.0)
      mat2 = 2*mat_b + 3*mat1
      call tc%assert_eq(d2%d, mat2, rtol, atol, "add band3")
    end block

    ! add sparse-sparse (store result in second matrix)
    block
      integer, allocatable :: ia_exp(:), ja_exp(:)
      real,    allocatable :: a_exp(:)
      type(sparse_real)    :: sA, sB, sE, sE2

      call get_test_matrix2(sA)
      call example_matrix3(sB)
      call get_empty_matrix(sE)

      a_exp  = sA%a
      ia_exp = int(sA%ia)
      ja_exp = sA%ja
      call matrix_add(sE, sA, fact=-1.5)
      call tc%assert_eq(a_exp,  sA%a,  0.0, 0.0, "add_sparse empty 1: a")
      call tc%assert_eq(ia_exp, int(sA%ia),      "add_sparse empty 1: ia")
      call tc%assert_eq(ja_exp, sA%ja,           "add_sparse empty 1: ja")

      call get_test_matrix2(sA)
      a_exp = 3 * sA%a
      call matrix_add(sA, sE, fact=3.0)
      call tc%assert_eq(a_exp,  sE%a,  0.0, 0.0, "add_sparse empty 2: a")
      call tc%assert_eq(ia_exp, int(sE%ia),      "add_sparse empty 2: ia")
      call tc%assert_eq(ja_exp, sE%ja,           "add_sparse empty 2: ja")

      call get_empty_matrix(sE)
      call get_empty_matrix(sE2)
      call matrix_add(sE2, sE, fact = 9.5)
      call tc%assert(sE%is_empty(), "add_sparse empty 3")

      ! A - 2 * B =
      !  -2    -4     5     1
      !   0    -6     0     0
      !   0     0    -2    -3
      !  -1    -2     0   -10
      !
      ! a  = [-2 -4 5 1 -6 -2 -3 -1 -2 -10]
      ! ia = [1 5 6 8 11]
      ! ja = [1 2 3 4 2 3 4 1 2 4]

      a_exp  = [-2,-4,5,1,-6,-2,-3,-1,-2,-10]
      ia_exp = [1,5,6,8,11]
      ja_exp = [1,2,3,4,2,3,4,1,2,4]
      call matrix_add(sB, sA, fact = -2.0)

      call tc%assert_eq(a_exp,  sA%a, 0.0, 0.0, "add_sparse: a")
      call tc%assert_eq(ia_exp, int(sA%ia),     "add_sparse: ia")
      call tc%assert_eq(ja_exp, sA%ja,          "add_sparse: ja")
    end block

    ! add sparse-sparse (store result in third matrix)
    block
      integer, allocatable :: ia_exp(:), ja_exp(:)
      real,    allocatable :: a_exp(:)
      type(sparse_real)    :: sA, sB, sC, sE

      call example_matrix3(sA)
      call get_test_matrix2(sB)
      call get_empty_matrix(sE)

      call matrix_add(sA, sE, sC, fact1 = 3.0, fact2 = -5.0)
      call tc%assert_eq(3.0*sA%a  ,     sC%a  , rtol, atol, "add_sparse3 empty 1: a")
      call tc%assert_eq(3.0*sA%a  ,     sC%a  , rtol, atol, "add_sparse3 empty 1: a")
      call tc%assert_eq(int(sA%ia), int(sC%ia),             "add_sparse3 empty 1: ia")
      call tc%assert_eq(    sA%ja ,     sC%ja ,             "add_sparse3 empty 1: ja")

      call matrix_add(sE, sA, sC, fact1=-5.0, fact2=3.0)
      call tc%assert_eq(3.0*sA%a  ,     sC%a  , rtol, atol, "add_sparse3 empty 2: a")
      call tc%assert_eq(3.0*sA%a  ,     sC%a  , rtol, atol, "add_sparse3 empty 2: a")
      call tc%assert_eq(int(sA%ia), int(sC%ia),             "add_sparse3 empty 2: ia")
      call tc%assert_eq(    sA%ja ,     sC%ja ,             "add_sparse3 empty 2: ja")

      ! 3*B-2*A=
      !  -2    -4    15     3
      !   0    -2     0     0
      !   0     0    -2    -9
      !  -3    -2     0   -10
      !
      ! a  = [-2 -4 15 3 -2 -2 -9 -3 -2 -10]
      ! ia = [1 5 6 8 11]
      ! ja = [1 2 3 4 2 3 4 1 2 4]
      a_exp  = [-2,-4,15,3,-2,-2,-9,-3,-2,-10]
      ia_exp = [1,5,6,8,11]
      ja_exp = [1,2,3,4,2,3,4,1,2,4]
      call matrix_add(sB, sA, sC, fact1 = 3.0, fact2 = -2.0)

      call tc%assert_eq(a_exp,  sC%a, 0.0, 0.0, "add_sparse3: a")
      call tc%assert_eq(ia_exp, int(sC%ia),     "add_sparse3: ia")
      call tc%assert_eq(ja_exp, sC%ja,          "add_sparse3: ja")
    end block

    ! add sparse-band (store result in second matrix)
    block
      integer, allocatable :: ia_exp(:), ja_exp(:)
      real,    allocatable :: a_exp(:)
      type(band_real)      :: B
      type(sparse_real)    :: S

      ! B =
      !  3    -2     0     0
      ! -1     3    -2     0
      !  0    -1     3    -2
      !  0     0    -1     3
      call B%init(4, 1, 1)
      B%d(-1,2:4) = -2.0
      B%d( 0, : ) = 3.0
      B%d(+1,1:3) = -1.0

      ! 2*B =
      !  6    -4     0     0
      ! -2     6    -4     0
      !  0    -2     6    -4
      !  0     0    -2     6
      a_exp = [6, -4, -2, 6, -4, -2, 6, -4, -2, 6]
      ia_exp = [1, 3, 6, 9, 11]
      ja_exp = [1, 2, 1, 2, 3, 2, 3, 4, 3, 4]
      call get_empty_matrix(S)
      call matrix_add(B, S, fact = 2.0)
      call tc%assert_eq(a_exp,  S%a, rtol, atol, "add_band empty: a")
      call tc%assert_eq(ia_exp, int(S%ia),       "add_band empty: ia")
      call tc%assert_eq(ja_exp, S%ja,            "add_band empty: ja")

      call example_matrix3(S)

      ! S+2*B=
      !  7    -2     0     0
      ! -2    10    -4     0
      !  0    -2     7    -4
      !  0     1    -2    11
      a_exp  = [7, -2, -2, 10, -4, -2, 7, -4, 1, -2, 11]
      ia_exp = [1, 3,  6,  9, 12]
      ja_exp = [1, 2, 1, 2, 3, 2, 3, 4, 2, 3, 4]
      call matrix_add(B, S, fact=2.0)
      call tc%assert_eq(a_exp,  S%a, 0.0, 0.0, "add_band: a")
      call tc%assert_eq(ia_exp, int(S%ia),     "add_band: ia")
      call tc%assert_eq(ja_exp, S%ja,          "add_band: ja")
    end block

    ! add sparse-band (store result in third matrix)
    block
      integer, allocatable :: ia_exp(:), ja_exp(:)
      real,    allocatable :: a_exp(:)
      type(band_real)      :: B
      type(sparse_real)    :: S, S2

      ! B =
      !  3    -2     0     0
      ! -1     3    -2     0
      !  0    -1     3    -2
      !  0     0    -1     3
      call B%init(4, 1, 1)
      B%d(-1,2:4) = -2.0
      B%d( 0, : ) = 3.0
      B%d(+1,1:3) = -1.0

      ! 2*B =
      !  6    -4     0     0
      ! -2     6    -4     0
      !  0    -2     6    -4
      !  0     0    -2     6
      a_exp = [6, -4, -2, 6, -4, -2, 6, -4, -2, 6]
      ia_exp = [1, 3, 6, 9, 11]
      ja_exp = [1, 2, 1, 2, 3, 2, 3, 4, 3, 4]
      call get_empty_matrix(S)
      call matrix_add(B, S, S2, fact1 = 2.0, fact2 = 1e99)
      call tc%assert_eq(a_exp,  S2%a, rtol, atol, "add_band3 empty: a")
      call tc%assert_eq(ia_exp, int(S2%ia),       "add_band3 empty: ia")
      call tc%assert_eq(ja_exp, S2%ja,            "add_band3 empty: ja")

      call example_matrix3(S)

      ! -1*S + 2*B=
      !  5    -6     0     0
      ! -2     2    -4     0
      !  0    -2     5    -4
      !  0    -1    -2     1
      a_exp  = [5, -6, -2, 2, -4, -2, 5, -4, -1, -2, 1]
      ia_exp = [1, 3,  6,  9, 12]
      ja_exp = [1, 2, 1, 2, 3, 2, 3, 4, 2, 3, 4]
      call matrix_add(B, S, S2, fact1 = 2.0, fact2 = -1.0)

      call tc%assert_eq(a_exp,  S2%a, rtol, atol, "add_band3: a")
      call tc%assert_eq(ia_exp, int(S2%ia),       "add_band3: ia")
      call tc%assert_eq(ja_exp, S2%ja,            "add_band3: ja")
    end block

    ! add block-block (store result in second or third matrix)
    block
      real                      :: d_exp(6,6), d1_exp(2,4)
      type(dense_real)          :: d_conv
      type(dense_real), pointer :: d1p
      type(block_real)          :: bl1, bl2, bl3

      ! get block matrix
      call example_matrix9(bl1)
      ! set dense pointer to non-owned dense block of example matrix
      select type (p => bl1%b(1,2)%p)
      type is (dense_real)
        d1p => p
      end select
      d1_exp = reshape([0.0, 1.0, 0.0, 0.0, &
                       &0.0, 1.0, 3.0, 0.0], shape(d1_exp), order = [2,1])

      ! add block-block with second matrix empty (result should be the same as copy_deep)
      call bl2%init([2,4])
      call matrix_add(bl1, bl2)
      call matrix_convert(bl2, d_conv)
      d_exp = reshape([1.0, 0.0, 0.0, 1.0, 0.0, 0.0, &
                      &0.0, 1.0, 0.0, 1.0, 3.0, 0.0, &
                      &2.0, 2.0, 1.0, 0.0, 1.0, 0.0, &
                      &0.0, 0.0, 0.0, 1.0, 1.0, 3.0, &
                      &0.0, 0.0, 0.0, 0.0, 1.0, 0.0, &
                      &0.0, 6.0, 0.0, 0.0, 0.0, 1.0], shape(d_exp), order = [2,1])
      call tc%assert_eq(d_conv%d, d_exp,  rtol, atol, "add block: empty matrix")
      call tc%assert_eq(d1p%d,    d1_exp, rtol, atol, "block: non-owned pointer") ! non-owned pointer should not change

      ! add block-block (store result in third matrix)
      call matrix_add(bl1, bl2, bl3, fact2 = 2.0)
      call d_conv%destruct()
      call matrix_convert(bl3, d_conv)
      d_exp = reshape([3.0, 0.0 , 0.0, 3.0, 0.0, 0.0, &
                      &0.0, 3.0 , 0.0, 3.0, 9.0, 0.0, &
                      &6.0, 6.0 , 3.0, 0.0, 3.0, 0.0, &
                      &0.0, 0.0 , 0.0, 3.0, 3.0, 9.0, &
                      &0.0, 0.0 , 0.0, 0.0, 3.0, 0.0, &
                      &0.0, 18.0, 0.0, 0.0, 0.0, 3.0], shape(d_exp), order = [2,1])
      call tc%assert_eq(d_conv%d, d_exp,  rtol, atol, "add block3")
      call tc%assert_eq(d1p%d,    d1_exp, rtol, atol, "block: non-owned pointer") ! non-owned pointer should not change

      ! add block-block (store result in second matrix)
      call matrix_add(bl1, bl3, fact = -1.0)
      call d_conv%destruct()
      call matrix_convert(bl3, d_conv)
      d_exp = reshape([2.0, 0.0 , 0.0, 2.0, 0.0, 0.0, &
                      &0.0, 2.0 , 0.0, 2.0, 6.0, 0.0, &
                      &4.0, 4.0 , 2.0, 0.0, 2.0, 0.0, &
                      &0.0, 0.0 , 0.0, 2.0, 2.0, 6.0, &
                      &0.0, 0.0 , 0.0, 0.0, 2.0, 0.0, &
                      &0.0, 12.0, 0.0, 0.0, 0.0, 2.0], shape(d_exp), order = [2,1])
      call tc%assert_eq(d_conv%d, d_exp,  rtol, atol, "add block2")
      call tc%assert_eq(d1p%d,    d1_exp, rtol, atol, "block: non-owned pointer") ! non-owned pointer should not change
    end block
  end subroutine

  subroutine test_mul(tc)
    !! test matrix_mul
    type(test_case), intent(inout) :: tc

    ! mul dense-dense
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

      call matrix_mul(d1, d2, d3)
      mat3 = matmul(mat1, mat2)
      call tc%assert_eq(d3%d, mat3, rtol, atol, "mul dense")
    end block

    ! mul sparse-sparse
    block
      integer, allocatable :: ia_exp(:), ja_exp(:)
      real,    allocatable :: a_exp(:)
      type(sparse_real)    :: sA, sB, sC, sE
      type(spbuild_real)   :: sbuild

      call get_empty_matrix(sE)
      call example_matrix3(sA)

      ! test 1: small matrices. check each entry
      ! B =
      !  0     0     5     1
      !  0     2     0     0
      !  0     0     0    -3
      ! -1     0     0     0

      call sB%init(4)
      call sbuild%init(sB)

      call sbuild%add(1, 3,  5.0)
      call sbuild%add(1, 4,  1.0)
      call sbuild%add(2, 2,  2.0)
      call sbuild%add(3, 4, -3.0)
      call sbuild%add(4, 1, -1.0)

      call sbuild%save()

      ! C = A*E
      call matrix_mul(sA, sE, sC)
      call tc%assert(sC%is_empty(), "mul sparse empty 1")

      ! C = E*A
      call matrix_mul(sE, sA, sC)
      call tc%assert(sC%is_empty(), "mul sparse empty 2")

      ! C = A*B =
      !  0     4     5     1
      !  0     8     0     0
      !  0     0     0    -3
      ! -5     2     0     0
      !
      ! a  = [4 5 1 8 -3 -5 2]
      ! ia = [1 4 5 6 8]
      ! ja = [2 3 4 2 4 1 2]
      !
      a_exp  = [4, 5, 1, 8, -3, -5, 2]
      ia_exp = [1, 4, 5, 6, 8]
      ja_exp = [2, 3, 4, 2, 4, 1, 2]

      call matrix_mul(sA, sB, sC)

      call tc%assert_eq(a_exp,  sC%a, 0.0, 0.0, "mul sparse 1: a")
      call tc%assert_eq(ia_exp, int(sC%ia),     "mul sparse 1: ia")
      call tc%assert_eq(ja_exp, sC%ja,          "mul sparse 1: ja")

      ! test 2: large matrices. check sum of entries
      ! C <- A*A
      block
        integer, parameter :: n=500

        call get_test_matrix3(n, sA)
        call get_test_matrix3(n, sB)
        call matrix_mul(sA, sB, sC)

        call tc%assert_eq(8.913712864525120e+14, sum(sC%a), rtol, atol, "mul sparse 2")
      end block

      ! test 3: large matrices in extern files
      ! C <- A*B
      call sA%input(file = "S1.test")
      call sB%input(file = "S2.test")
      call sE%input(file = "S3.test") ! expected S1 * S2
      call matrix_mul(sA, sB, sC)
      call tc%assert_eq(int(sE%ia), int(sC%ia),  "mul sparse 3: ia")
      call tc%assert_eq(sE%ja, sC%ja,            "mul sparse 3: ja")
      call tc%assert_eq(sE%a,  sC%a, rtol, atol, "mul sparse 3: a")
    end block
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
      call tc%assert_eq(exp_arr_r, arr_r, rtol, atol, "diag: dense real -> array")

      ! dense complex -> array
      call matr_c%init(matr_in_c)
      call matrix_diag(matr_c, arr_c)
      call tc%assert_eq(exp_arr_c, arr_c, rtol, atol, "diag: dense complex -> array")

      ! array -> dense real
      call matrix_diag(exp_arr_r, matr_r)
      call tc%assert_eq(matr_out_r, matr_r%d, rtol, atol, "diag: array -> dense real")

      ! array -> dense complex
      call matrix_diag(exp_arr_c, matr_c)
      call tc%assert_eq(matr_out_c, matr_c%d, rtol, atol, "diag: array -> dense complex")
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
      call tc%assert_eq(exp_arr_r, arr_r, rtol, atol, "diag: sparse real -> array"   )

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
      call tc%assert_eq(exp_arr_c, arr_c, rtol, atol, "diag: sparse complex -> array")

      ! array -> sparse real
      call matrix_diag(exp_arr_r, s_r)
      call tc%assert_eq(exp_arr_r,    s_r%a, rtol, atol, "diag: array -> sparse real: a" )
      call tc%assert_eq([(i, i=1,4)], int(s_r%ia),  "diag: array -> sparse real: ia")
      call tc%assert_eq([(i, i=1,3)], s_r%ja,       "diag: array -> sparse real: ja")

      ! array -> sparse complex
      call matrix_diag(exp_arr_c, s_c)
      call tc%assert_eq(exp_arr_c,    s_c%a, rtol, atol, "diag: array -> sparse complex: a" )
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
        call tc%assert_eq(d_out(:,i), y, rtol, atol, "approx: sparse")
      end do
    end block
  end subroutine

end submodule
