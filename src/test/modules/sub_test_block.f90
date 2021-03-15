submodule(test_matrix_m) test_block_m
  !! TODO fixme
  !!    test if a block of a block matrix could be another block matrix (instead of the usual dense or so)
  use matrix_m
  implicit none

contains

  module subroutine test_block()
    type(test_case) :: tc

    integer                  :: i
    real, allocatable        :: d0(:,:)
    type(block_real), target :: M

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
    call matrix1(M)

    ! mul_vec
    block
      real :: x(9), y(9), y_exp(9), fact = -3

      x     = [(i, i=1, 9)]
      y     = [(i, i=11, 19)]
      y_exp = [-24,-3,20,37,85,63,-12,31,35]

      call M%mul_vec(x, y, fact_y=fact)

      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec")
    end block

    ! delete_block
    block
      type(block_real)           :: A
      type(spbuild_real)         :: Asb
      type(sparse_real)          :: Asp
      type(sparse_real), pointer :: Abl

      !     //  \       \
      ! A = |\  /       |
      !     |/1  \ /1  \|
      !     \\  1/ \  1//
      !
      !                        /           \
      ! 2x delete_block -> A = |           |
      !                        |/1  \      |
      !                        \\  1/      /

      call A%init([2, 2])

      call A%get(1, 1, Abl)   ! A allocates block
      call A%get(2, 1, Abl)   ! A allocates block
      Abl = sparse_eye_real(2)
      call A%get(2, 2, Abl)   ! A allocates block
      Abl = sparse_eye_real(2)

      call A%delete_block(1, 1)      ! necessary, otherwise error
      call A%delete_block(2, 2)      ! just here for testing

      call Asp%init(4)
      call Asb%init(Asp)
      call matrix_convert(A, Asb)
      call Asb%save()

      call tc%assert_eq([1.0, 1.0],      Asp%a, 1e-12, "delete block: a" )
      call tc%assert_eq([1, 1, 1, 2, 3], int(Asp%ia),  "delete block: ia")
      call tc%assert_eq([1, 2],          Asp%ja,       "delete block: ja")
    end block

    ! set_ptr
    block
      type(dense_real), target :: ext_A

      d0 = reshape([1,0,0,2,3,5,0,4,6], [3,3])
      call ext_A%init(d0)
      call M%set_ptr(1, 1, ext_A)
      call M%set_ptr(1, 2, ext_A)
      call M%set_ptr(2, 1, ext_A)

      ! destruct should not delete ext_A
      call M%destruct()

      call tc%assert_eq(3, ext_A%nrows, "set_ptr")
    end block

    call tc%finish
  end subroutine

end submodule
