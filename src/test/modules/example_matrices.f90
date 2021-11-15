module example_matrices_m

  use matrix_m, only: band_real, band_eye_real, block_real, dense_real, dense_cmplx, sparse_real, sparse_cmplx, spbuild_real, spbuild_cmplx

  implicit none

  private
  public matrix1
  public matrix2
  public example_matrix3
  public example_matrix4
  public example_matrix5
  public example_matrix6
  public example_matrix7
  public example_matrix8

  interface matrix1
    !! full 9x9 matrix.
    module procedure :: matrix1_block_real
    module procedure :: matrix1_dense_real
  end interface

  interface matrix2
    !! 5x5 matrix
    !!  [ 11    7    3    0    0
    !!    16   12    8    4    0
    !!     0   17   13    9    5
    !!     0    0   18   14   10
    !!     0    0    0   19   15 ]
    !! octave:
    !!    = diag(3:5, 2) + diag(7:10, 1) + diag(11:15) + diag(16:19, -1)
    module procedure :: matrix2_band_real
    module procedure :: matrix2_dense_real
    module procedure :: matrix2_sparse_real
  end interface

  interface example_matrix3
    !! 4x4 matrix
    !!  1     2     0     0
    !!  0     4     0     0
    !!  0     0     1     0
    !!  0     1     0     5
    !! octave:
    !!    = [1 2 0 0; 0 4 0 0; 0 0 1 0; 0 1 0 5]
    module procedure :: example_matrix3_sparse_real
  end interface

  interface example_matrix4
    !! 600x600 hermitian complex matrix
    !! 3000 entries -> sparse
    !! taken from FEAST example folder
    module procedure :: example_matrix4_sparse_cmplx
  end interface

  interface example_matrix5
    !! 4x4 matrix
    !! 1+i   0    -2-3i   0
    !!  0   2+0i    0     0
    !!  0    0     0-i    0
    !!  0    0      0    0+i
    !! octave:
    !!    = [(1,1) 0 (-2,-3) 0; 0 (2,0) 0 0; 0 0 (0,-1) 0; 0 0 0 (0,1)]
    module procedure :: example_matrix5_sparse_cmplx
  end interface

  interface example_matrix6
    ! 2x2 matrix
    ! 1  i
    ! i  1
    module procedure :: example_matrix6_dense_cmplx
  end interface

  interface example_matrix7
    !! 3x3 Matrix
    module procedure :: example_matrix7_block_real
  end interface

  interface example_matrix8
    !! 3x3 Matrix
    module procedure :: example_matrix8_block_real
  end interface

contains

  subroutine matrix1_block_real(bl)
    !! create block matrix M for futher tests: mulvec, etc.
    !! Example
    !!        / A B C \
    !!    M = | D E F | for A,..,I is (3x3)-matrix
    !!        \ G H I /
    !!
    !!    => nbrows=3, nbcols=3, i0=[1,4,7], i1=[3,6,9], j0=[1,4,7], j1=[3,6,9]

    type(block_real), intent(out), target :: bl

    integer           :: i, j
    real, allocatable :: d0(:,:)

    call bl%init([3, 3, 3], tridiag=.true.)

    ! set dense A
    block
      type(dense_real), pointer :: A

      ! get pointer into bl's blocks
      call bl%get(1, 1, A)

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

      ! get pointer into bl's blocks
      call bl%get(1, 2, B)

      ! B = sparse =
      !   1     0     0
      !   0     3     0
      !   1     3     2
      d0 = reshape([1,0,1,0,3,3,0,0,2], [3,3])

      call B%init(3)
      call sb%init(B)

      do i = 1, 3
        do j = 1, 3
          if (d0(i,j) == 0.0) cycle
          call sb%add(i, j, d0(i,j))
        end do
      end do
      call sb%save
    end block

    ! set D: band diag matrix
    block
      type(band_real), pointer :: D

      ! get pointer into bl's blocks
      call bl%get(2, 1, D)

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

      ! get pointer into bl's blocks
      call bl%get(2, 2, E)

      ! E =
      !   4     2     0
      !   0     5     3
      !   0     0     6
      d0 = reshape([(i, i=1,6)], [2,3], order=[2,1])
      call E%init(3, 0, nupper=1, d0=d0)
    end block

    ! set F: dense matrix
    block
      type(dense_real) :: F

      ! F =
      !   1 2 3
      !   0 4 5
      !   0 0 6
      call F%init(3)
      F%d = reshape([1,2,3,0,4,5,0,0,6], [3, 3], order=[2, 1])

      call bl%set(2, 3, F)
    end block

    ! set H: dense matrix
    block
      type(dense_real) :: H

      ! H =
      !   1 2 3
      !   4 5 6
      !   0 7 8
      call H%init(3)
      H%d = reshape([1,2,3,4,5,6,0,7,8], [3, 3], order=[2, 1])

      call bl%set(3, 2, H)
    end block

    ! set I: eye matrix
    block
      ! I =
      !   1 0 0
      !   0 1 0
      !   0 0 1
      call bl%set(3, 3, band_eye_real(3))
    end block
  end subroutine

  subroutine matrix1_dense_real(d)
    type(dense_real), intent(out) :: d

    call d%init(real(reshape([     &
      & 1, 2, 0, 1, 0, 0, 0, 0, 0, &
      & 0, 3, 4, 0, 3, 0, 0, 0, 0, &
      & 0, 5, 6, 1, 3, 2, 0, 0, 0, &
      & 3, 0, 0, 4, 2, 0, 1, 2, 3, &
      & 0, 5, 0, 0, 5, 3, 0, 4, 5, &
      & 0, 0, 7, 0, 0, 6, 0, 0, 6, &
      & 0, 0, 0, 1, 2, 3, 1, 0, 0, &
      & 0, 0, 0, 4, 5, 6, 0, 1, 0, &
      & 0, 0, 0, 0, 7, 8, 0, 0, 1  ], [9, 9], order=[2, 1])))
  end subroutine

  subroutine matrix2_band_real(b)
    type(band_real), intent(out) :: b

    integer :: i

    call b%init(5, 1, nupper=2, d0=reshape([(real(i), i=1,20)], [4, 5], order=[2, 1]))
  end subroutine

  subroutine matrix2_dense_real(d)
    type(dense_real), intent(out) :: d

    call d%init(real(reshape([  &
      & 11, 16,  0,  0,  0,     &
      &  7, 12, 17,  0,  0,     &
      &  3,  8, 13, 18,  0,     &
      &  0,  4,  9, 14, 19,     &
      &  0,  0,  5, 10, 15      ], [5, 5])))
  end subroutine

  subroutine matrix2_sparse_real(s, fact_lodiag)
    type(sparse_real), intent(out) :: s
    real, optional,    intent(in) :: fact_lodiag
      !! scaling factor for lower diagonal

    type(spbuild_real) :: sb

    !  [ 11       7       3       0       0
    !    16*fact   12     8       4       0
    !     0      17*fact 13       9       5
    !     0       0      18*fact 14      10
    !     0       0       0      19*fact 15 ]

    call s%init(5)
    call sb%init(s)

    call sb%add(1, 1, 11.0)
    call sb%add(1, 2,  7.0)
    call sb%add(1, 3,  3.0)

    call sb%add(2, 1, 16.0*fact_lodiag)
    call sb%add(2, 2, 12.0)
    call sb%add(2, 3,  8.0)
    call sb%add(2, 4,  4.0)

    call sb%add(3, 2, 17.0*fact_lodiag)
    call sb%add(3, 3, 13.0)
    call sb%add(3, 4,  9.0)
    call sb%add(3, 5,  5.0)

    call sb%add(4, 3, 18.0*fact_lodiag)
    call sb%add(4, 4, 14.0)
    call sb%add(4, 5, 10.0)

    call sb%add(5, 4, 19.0*fact_lodiag)
    call sb%add(5, 5, 15.0)

    call sb%save()
  end subroutine

  subroutine example_matrix3_sparse_real(S)
    type(sparse_real), intent(out) :: S

    type(spbuild_real) :: sbuild

    ! S =
    !  1     2     0     0
    !  0     4     0     0
    !  0     0     1     0
    !  0     1     0     5
    !
    ! a  = [1 2 4 1 1 5]
    ! ia = [1 3 4 5 7]
    ! ja = [1 2 2 3 2 4]

    call S%init(4)
    call sbuild%init(S)

    call sbuild%add(1, 1, 1.0)
    call sbuild%add(1, 2, 2.0)
    call sbuild%add(2, 2, 4.0)
    call sbuild%add(3, 3, 1.0)
    call sbuild%add(4, 2, 1.0)
    call sbuild%add(4, 4, 5.0)

    call sbuild%save()
  end subroutine

  subroutine example_matrix4_sparse_cmplx(s)
    !! load example matrix 4 as complex sparse matrix
    type(sparse_cmplx), intent(out) :: s

    integer             :: iounit, i, n, row, rows, col, cols
    real                :: re, im
    type(spbuild_cmplx) :: sb

    open (newunit=iounit, file='src/test/example_matrix4.dat', action='read')
    read (iounit, *)
    read (iounit, *) rows, cols, n
    call s%init(rows)
    call sb%init(s)

    ! load values and save into spbuild
    do i = 1, n
      read (iounit, *) row, col, re, im
      call sb%add(row, col, cmplx(re, y=im))
    end do

    close (unit=iounit)
    call sb%save()
  end subroutine

  subroutine example_matrix5_sparse_cmplx(s)
    !! load example matrix 5 as complex sparse matrix
    type(sparse_cmplx), intent(out) :: s

    type(spbuild_cmplx) :: sb

    ! 4x4 matrix
    ! 1+i   0    -2-3i   0
    !  0   2+0i    0     0
    !  0    0     0-i    0
    !  0    0      0    0+i
    !
    ! a  = [(1,1), (-2,-3), (2,0), (0,-1), (1,0)]
    ! ia = [1, 3, 2, 3, 4]
    ! ja = [1, 3, 4, 5, 6]

    call s%init(4)
    call sb%init(s)

    call sb%add(1, 1, ( 1,  1))
    call sb%add(1, 3, (-2, -3))
    call sb%add(2, 2, ( 2,  0))
    call sb%add(3, 3, ( 0, -1))
    call sb%add(4, 4, ( 0,  1))
    call sb%save()
  end subroutine

  subroutine example_matrix6_dense_cmplx(d)
    !! load example matrix 6 as complex dense matrix
    type(dense_cmplx), intent(out) :: d

    ! 2x2 matrix
    ! 1  i
    ! i  1

    call d%init(reshape([&
      & (1, 0), (0, 1),  &
      & (0, 1), (1,  0)  ], [2, 2], order = [2, 1]))
  end subroutine

  subroutine example_matrix7_block_real(M)
    !! create block matrix M for some tests: mulvec, etc.
    !! Example
    !!        / A B \
    !!    M = \ C D / for A is (1,1), B is (1,2),
    !!                    C is (2,1) and D is (2,2)-matrix
    !!        / -1  1  1 \
    !!    M = |  1  0  1 |
    !!        \  1  1  0 /
    !!
    !!    => nbrows=[1, 2], nbcols=[1, 2], i0=[1,2], i1=[1,3], j0=[1,2], j1=[1,3]
    type(block_real), intent(out) :: M

    call M%init([1, 2])

    block
      type(dense_real), pointer :: A

      call M%get(1, 1, A)
      call A%init(real(reshape([-1], [1, 1])))
    end block

    block
      type(dense_real), pointer :: B

      call M%get(1, 2, B)
      call B%init(real(reshape([1, 1], [1, 2])))
    end block

    block
      type(dense_real), pointer :: C

      call M%get(2, 1, C)
      call C%init(real(reshape([1, 1], [2, 1])))
    end block

    block
      type(dense_real), pointer :: D

      call M%get(2, 2, D)
      call D%init(real(reshape([ &
                          &0, 1, &
                          &1, 0  ], [2, 2], order = [2, 1])))
    end block
  end subroutine

  subroutine example_matrix8_block_real(M)
    !! create block matrix M for some tests: mulvec, etc.
    !! Example
    !!        / A 0 \
    !!    M = \ 0 B / for A and B are (2,2) matrices
    !!
    !!        /  1  2  0  0 \
    !!    M = |  2  1  0  0 |
    !!        |  0  0 -1  3 |
    !!        \  0  0  3 -1 /
    !!
    !!
    type(block_real), intent(out) :: M

    call M%init([2, 2], diag = .true.)

    block
      type(dense_real), pointer :: A

      call M%get(1, 1, A)
      call A%init(real(reshape([&
                        & 1, 2, &
                        & 2, 1  ], [2, 2])))
    end block

    block
      type(dense_real), pointer :: B

      call M%get(2, 2, B)
      call B%init(real(reshape([&
                      & -1,  3, &
                      &  3, -1  ], [2, 2])))
    end block
  end subroutine

end module
