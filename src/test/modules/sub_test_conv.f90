#include "../../util/macro.f90.inc"

submodule (test_matrix_m) test_conv_m

  implicit none

contains

  module subroutine test_conv()
    type(test_case) :: tc

    call tc%init("conv")

    call test_to_band(tc)
    call test_to_dense(tc)
    call test_to_sparse(tc)
    call test_to_real(tc)
    call test_to_cmplx(tc)

    call tc%finish()
  end subroutine

  subroutine test_to_band(tc)
    type(test_case), intent(inout) :: tc

    ! sparse
    block
      type(band_real)   :: b
      type(sparse_real) :: s

      call example_matrix3(s)
      !
      ! test 1: set band equal sparse, not insertion as block anywhere
      !
      call b%init(4, 2, nupper=1)
      call matrix_convert(s, b)

      ! just for testing purposes
      block
        real, allocatable :: d_exp(:,:)

        allocate (d_exp(4,4))
        d_exp(1,:) = [0, 2, 0, 0]
        d_exp(2,:) = [1, 4, 1, 5]
        d_exp(3,:) = 0
        d_exp(4,:) = [0, 1, 0, 0]

        call tc%assert_eq(d_exp, b%d, 1e-12, "to_band: no insertion arguments")
      end block

      !
      ! test 2: insert sparse into band
      !
      call b%init(7, 3, nupper=0)
      b%d = -1                          ! just for testing purposes: to see which data will get overwritten
      call matrix_convert(s, b, i0=3, j0=2)

      ! just for testing purposes
      block
        integer           :: i, j
        real, allocatable :: d_exp(:,:)
        type(dense_real)  :: d

        allocate (d_exp(7,7), source=0.0)

        ! set default value "-1" for bands
        do i = 1, 7
          do j = 1, 7
            if (j > i) exit
            if (j < i-3) cycle
            d_exp(i,j) = -1
          end do
        end do

        ! insert values from test matrix
        d_exp(3  ,2:3) = [1, 2]
        d_exp(4  ,3  ) =     4
        d_exp(5  ,4  ) =        1
        d_exp(6  ,3  ) =     1
        d_exp(6  ,5  ) =          5

        call d%init(7)
        call matrix_convert(b, d)
        call tc%assert_eq(d_exp, d%d, 1e-12, "to_band: insert sparse into band")
      end block
    end block
  end subroutine

  subroutine test_to_dense(tc)
    type(test_case), intent(inout) :: tc

    ! band
    block
      type(band_real)   :: b
      type(dense_real)  :: d, b2d
      real, allocatable :: d_exp(:,:)

      ! band = diag(7:10,1) + diag(11:15) + diag(16:19,-1) + diag(3:5,2)
      !  11    7    3    0    0
      !  16   12    8    4    0
      !   0   17   13    9    5
      !   0    0   18   14   10
      !   0    0    0   19   15

      ! setup band matrix
      call matrix2(b)

      ! conversion
      call d%init(7)
      call matrix_convert(b, d, i0=2, j0=2)

      ! expected matrix
      call matrix2(b2d)
      allocate (d_exp(7,7), source=0.0)
      d_exp(2:6,2:6) = b2d%d

      call tc%assert_eq(d_exp, d%d, 0.0, "convert: band -> dense")
    end block

    ! block
    block
      real, allocatable        :: d_exp(:,:)
      type(block_real), target :: bl0
      type(dense_real)         :: d0, d

      call matrix1(bl0)
      call d%init(11)
      call matrix_convert(bl0, d, i0=2, j0=3)

      call matrix1(d0)
      allocate (d_exp(11,11), source=0.0)
      d_exp(2:10,3:11) = d0%d

      call tc%assert_eq(d_exp, d%d, 0.0, "convert: block -> dense")
    end block

    ! dense
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

      call matrix_convert(d1, d2, i0=1, j0=3)

      do i = 1, 5
        do j = 1, 5
          val = 2.0
          if (i <= 3 .and. j >= 3) val = 1.0
          call tc%assert_eq(d2%d(i,j), val, 1e-12, "convert: dense -> dense")
        end do
      end do
    end block

  end subroutine

  subroutine test_to_sparse(tc)
    type(test_case), intent(inout) :: tc

    integer, allocatable :: ia_exp(:), ja_exp(:)
    real,    allocatable :: a_exp(:)
    type(sparse_real)    :: s
    type(spbuild_real)   :: sb

    ! band
    block
      type(band_real) :: b

      call matrix2(b)
      call s%init(7)
      call sb%init(s)
      call matrix_convert(b, sb, 2, 2)
      call sb%save()

      ! s(2:6,2:6) =
      !  11    7    3    0    0
      !  16   12    8    4    0
      !   0   17   13    9    5
      !   0    0   18   14   10
      !   0    0    0   19   15
      a_exp  = [11, 7, 3, 16, 12, 8, 4, 17, 13, 9, 5, 18, 14, 10, 19, 15]
      ia_exp = [1, 1, 4, 8, 12, 15, 17, 17]
      ja_exp = [2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 5, 6]

      call tc%assert_eq(a_exp,      s%a, 1e-16, "to sparse: a" )
      call tc%assert_eq(ia_exp, int(s%ia),      "to sparse: ia")
      call tc%assert_eq(ja_exp,     s%ja,       "to sparse: ja")
    end block

    ! block
    block
      type(block_real) :: b

      call matrix1(b)
      call s%init(9)
      call sb%init(s)
      call matrix_convert(b, sb, drop_zeros=.true.)
      call sb%save()

      ! M =
      !     1     2     0     1     0     0     0     0     0
      !     0     3     4     0     3     0     0     0     0
      !     0     5     6     1     3     2     0     0     0
      !     3     0     0     4     2     0     1     2     3
      !     0     5     0     0     5     3     0     4     5
      !     0     0     7     0     0     6     0     0     6
      !     0     0     0     1     2     3     1     0     0
      !     0     0     0     4     5     6     0     1     0
      !     0     0     0     0     7     8     0     0     1
      a_exp = [1, 2, 1,          &
        &      3, 4, 3,          &
        &      5, 6, 1, 3, 2,    &
        &      3, 4, 2, 1, 2, 3, &
        &      5, 5, 3, 4, 5,    &
        &      7, 6, 6,          &
        &      1, 2, 3, 1,       &
        &      4, 5, 6, 1,       &
        &      7, 8, 1           ]
      ia_exp = [1, 4, 7, 12, 18, 23, 26, 30, 34, 37]
      ja_exp = [1, 2, 4,          &
        &       2, 3, 5,          &
        &       2, 3, 4, 5, 6,    &
        &       1, 4, 5, 7, 8, 9, &
        &       2, 5, 6, 8, 9,    &
        &       3, 6, 9,          &
        &       4, 5, 6, 7,       &
        &       4, 5, 6, 8,       &
        &       5, 6, 9           ]

      call tc%assert_eq(a_exp,      s%a, 1e-16, "to sparse: a" )
      call tc%assert_eq(ia_exp, int(s%ia),      "to sparse: ia")
      call tc%assert_eq(ja_exp,     s%ja,       "to sparse: ja")
    end block

    ! dense
    block
      integer              :: i, j
      logical, allocatable :: sparsity(:,:)
      real                 :: mat(5,5), e(5), v(5), v_exp(5)
      real,    allocatable :: tmp(:,:)
      type(dense_real)     :: d

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

      call matrix_convert(d, sb, i0=1, j0=3)
      call sb%save()

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

      call matrix_convert(d, sb, i0=4, j0=2)
      call s%init(5)
      call sb%save()
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
      call matrix_convert(d, sb, i0=1, j0=2, struct=sparsity)
      call sb%save()

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
      call matrix_convert(d, sb, i0=2, j0=1, drop_zeros=.true.)
      call sb%save()

      do i = 1, 5
        e     = 0.0
        e(i)  = 1.0
        v_exp = matmul(mat, e)
        call s%mul_vec(e, v)
        call tc%assert_eq(v, v_exp, 1e-12, "to_sparse: drop zeros")
      end do
      call tc%assert_eq(count(s%a == 0.0), 0, "to_sparse: drop zeros")
    end block
  end subroutine

  subroutine test_to_real(tc)
    type(test_case), intent(inout) :: tc

    ! fixme
    IGNORE(tc)
  end subroutine

  subroutine test_to_cmplx(tc)
    type(test_case), intent(inout) :: tc

    ! fixme
    IGNORE(tc)
  end subroutine

end submodule
