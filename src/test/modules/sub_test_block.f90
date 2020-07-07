submodule(test_matrix_m) test_block_m
  !! TODO fixme
  !!    test if a block of a block matrix could be another block matrix (instead of the usual dense or so)
  use matrix_m
  implicit none

contains

  module subroutine test_block
    type(test_case)   :: tc
    type(block_real)  :: M
    real, allocatable :: d0(:,:)
    integer :: i, j

    print "(A)", "test_block"
    call tc%init("block")

    ! test init routine
    block
      type(block_real) :: A
      integer          :: row_dim(3), col_dim(4), i0(3), i1(3), j0(4), j1(4)

      row_dim = [2, 3, 1]
      col_dim = [2, 2, 4, 3]

      i0 = [1, 3, 6]
      i1 = [2,5,6]
      j0 = [1,3,5,9]
      j1 = [2,4,8,11]

      call A%init(row_dim, col_dim=col_dim)

      call tc%assert_eq(i0, A%i0, "init: i0")
      call tc%assert_eq(i1, A%i1, "init: i1")
      call tc%assert_eq(j0, A%j0, "init: j0")
      call tc%assert_eq(j1, A%j1, "init: j1")
    end block

    ! create block matrix M for futher tests: mulvec, etc.
    block
      ! Example
      !        / A B C \
      !    M = | D E F | for A,..,I is (3x3)-matrix
      !        \ G H I /
      !
      !    => nbrows=3, nbcols=3, i0=[1,4,7], i1=[3,6,9], j0=[1,4,7], j1=[3,6,9]
      call M%init([3,3,3], tridiag=.true.)

      ! set dense A
      block
        type(dense_real), pointer :: A

        ! get pointer into M's blocks
        call M%get(1, 1, A)

        ! A=dense=
        !   1 2 0
        !   0 3 4
        !   0 5 6
        d0 = reshape([1,0,0,2,3,5,0,4,6], [3,3])
        call A%init(d0)
      end block

      ! set sparse B
      block
        type(sparse_real), pointer :: B
        type(spbuild_real)         :: sb

        ! get pointer into M's blocks
        call M%get(1, 2, B)

        ! B = sparse =
        !   1     0     0
        !   0     3     0
        !   1     3     2
        d0 = reshape([1,0,1,0,3,3,0,0,2], [3,3])

        call B%init(3)
        call sb%init(B)

        do i = 1 ,3
        do j = 1 ,3
          if (d0(i,j) == 0.0) cycle
          call sb%add(i, j, d0(i,j))
        end do
        end do
        call sb%save
      end block

      ! set D: band diag matrix
      block
        type(band_real), pointer :: D

        ! get pointer into M's blocks
        call M%get(2, 1, D)

        ! D =
        !   3     0     0
        !   0     5     0
        !   0     0     7
        d0 = reshape([3,5,7], [1,3])
        call D%init(3, 0, d0=d0)
      end block

      ! set E: band matrix
      block
        type(band_real), pointer :: E

        ! get pointer into M's blocks
        call M%get(2, 2, E)

        ! E =
        !   4     2     0
        !   0     5     3
        !   0     0     6
        d0 = reshape([(i, i=1,6)], [2,3], order=[2,1])
        call E%init(3, 0, nupper=1, d0=d0)
      end block

      ! set F: triangular matrix
      block
        type(triang_real) :: F

        ! F =
        !   1 2 3
        !   0 4 5
        !   0 0 6
        call F%init(3, .true.)
        F%d = reshape([1,2,3,0,4,5,0,0,6], [3, 3], order=[2, 1])

        call M%set(2, 3, F)
      end block

      ! set H: Hessenberg matrix
      block
        type(hessenberg_real) :: H

        ! H =
        !   1 2 3
        !   4 5 6
        !   0 7 8
        call H%init(3, .true.)
        H%d = reshape([1,2,3,4,5,6,0,7,8], [3, 3], order=[2, 1])

        call M%set(3, 2, H)
      end block

      ! set I: eye matrix
      block
        ! I =
        !   1 0 0
        !   0 1 0
        !   0 0 1
        call M%set(3, 3, band_eye_real(3))
      end block
    end block

    ! mul_vec
    block
      real :: x(9), y(9), y_exp(9), fact = -3

      x     = [(i, i=1, 9)]
      y     = [(i, i=11, 19)]
      y_exp = [-24,-3,20,37,85,63,-12,31,35]

      call M%mul_vec(x, y, fact_y=fact)

      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec")
    end block

    ! to dense
    block
      type(dense_real) :: d
      real             :: d_exp(9,9)

      call d%init(9)
      call M%to_dense(d)

      d_exp = reshape([&
        & 1, 2, 0, 1, 0, 0, 0, 0, 0, &
        & 0, 3, 4, 0, 3, 0, 0, 0, 0, &
        & 0, 5, 6, 1, 3, 2, 0, 0, 0, &
        & 3, 0, 0, 4, 2, 0, 1, 2, 3, &
        & 0, 5, 0, 0, 5, 3, 0, 4, 5, &
        & 0, 0, 7, 0, 0, 6, 0, 0, 6, &
        & 0, 0, 0, 1, 2, 3, 1, 0, 0, &
        & 0, 0, 0, 4, 5, 6, 0, 1, 0, &
        & 0, 0, 0, 0, 7, 8, 0, 0, 1  ], [9,9], order=[2,1])

      call tc%assert_eq(d_exp, d%d, 1e-12, "to dense")
    end block

    ! to sparse
    block
      type(sparse_real)    :: s
      type(spbuild_real)   :: sb
      real,    allocatable :: a_exp(:)
      integer, allocatable :: ia_exp(:), ja_exp(:)

      call s%init(9)
      call sb%init(s)

      call M%to_sparse(sb, drop_zeros=.true.)
      call sb%save

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

      call tc%assert_eq(a_exp,  s%a, 1e-12, "to sparse: a")
      call tc%assert_eq(ia_exp, s%ia,       "to sparse: ia")
      call tc%assert_eq(ja_exp, s%ja,       "to sparse: ja")
    end block

    call tc%finish
  end subroutine

end submodule
