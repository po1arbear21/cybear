submodule (test_matrix_m) test_block_m
  !! TODO fixme
  !!    test if a block of a block matrix could be another block matrix (instead of the usual dense or so)

  implicit none

  real, parameter :: rtol = 1e-14
  real, parameter :: atol = 1e-16

contains

  module subroutine test_block()
    type(test_case) :: tc

    real, allocatable :: d0(:,:)
    type(block_real)  :: M1, M2

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

      call A%init(row_dim, col_dim = col_dim)

      call tc%assert_eq(i0, A%i0, "init: i0")
      call tc%assert_eq(i1, A%i1, "init: i1")
      call tc%assert_eq(j0, A%j0, "init: j0")
      call tc%assert_eq(j1, A%j1, "init: j1")
    end block

    ! create block matrix M1 for futher tests: mulvec, etc.
    call example_matrix7(M1)

    ! scale, get
    block
      real, parameter           :: l = 3.0
      real                      :: A_exp(1,1), B_exp(1,2), C_exp(2,1), D_exp(2,2)
      type(dense_real), pointer :: A, B, C, D
      type(block_real)          :: M

      call example_matrix7(M)

      A_exp = -1
      B_exp =  1
      C_exp =  1
      D_exp = reshape([&
                        & 0, 1, &
                        & 1, 0  ], [2, 2], order = [2, 1])

      call M%scale(l)

      call M%get(1, 1, A)
      call M%get(1, 2, B)
      call M%get(2, 1, C)
      call M%get(2, 2, D)

      call tc%assert_eq(A_exp * l, A%d, 0.0, 0.0, "scale: A")
      call tc%assert_eq(B_exp * l, B%d, 0.0, 0.0, "scale: B")
      call tc%assert_eq(C_exp * l, C%d, 0.0, 0.0, "scale: C")
      call tc%assert_eq(D_exp * l, D%d, 0.0, 0.0, "scale: D")
    end block

    ! mul_vec
    block
      real :: x(3), y(3), y_exp(3)

      x     = [1, 3, 1]
      y_exp = [3, 2, 4]
      call M1%mul_vec(x, y)
      call tc%assert_eq(y_exp, y, 0.0, 0.0, "mul_vec")

      call M1%mul_vec(x, y, trans = "T")
      call tc%assert_eq(y_exp, y, 0.0, 0.0, "mul_vec T")

      y     = [1,  1, 0]
      y_exp = [0, -1, 4]
      call M1%mul_vec(x, y, fact_y = -3.0)
      call tc%assert_eq(y_exp, y, 0.0, 0.0, "mul_vec fact")

      y = [1, 1, 0]
      call M1%mul_vec(x, y, fact_y = -3.0, trans = "T")
      call tc%assert_eq(y_exp, y, 0.0, 0.0, "mul_vec fact T")
    end block

    ! mul_mat
    block
      real :: x(3,3), y(3,3), y_exp(3,3)

      x = reshape([&
        & 1, 0, 1, &
        & 3, 1, 3, &
        & 1, 0, 1  ], [3, 3], order = [2, 1])

      y_exp = reshape([&
        & 3, 1, 3, &
        & 2, 0, 2, &
        & 4, 1, 4  ], [3, 3], order = [2, 1])

      call M1%mul_mat(x, y)
      call tc%assert_eq(y_exp, y, rtol, atol, "mul_mat")

      call M1%mul_mat(x, y, trans = "T")
      call tc%assert_eq(y_exp, y, rtol, atol, "mul_mat T")

      call M1%mul_mat(x, y, fact_y = -2.0)
      call tc%assert_eq(-1 * y_exp, y, rtol, atol, "mul_mat fact")

      call M1%mul_mat(x, y, fact_y = -2.0)
      call tc%assert_eq(3 * y_exp, y, rtol, atol, "mul_mat fact T")
    end block

    call example_matrix8(M2)
    call M2%factorize()

    ! factorize, solve_vec
    block
      real :: x(4), x_exp(4), y(4)

      y     = [-1,  1,  4, -4]
      x_exp = [ 1, -1, -1,  1]
      call M2%solve_vec(y, x)
      call tc%assert_eq(x_exp, x, rtol, atol, "solve_vec")

      call M2%solve_vec(y, x, trans = "T")
      call tc%assert_eq(x_exp, x, rtol, atol, "solve_vec T")
    end block

    ! factorize, solve_mat
    block
      real :: x(4,2), x_exp(4,2), y(4,2)

      y     = reshape([&
             &  3,  2, &
             &  3,  1, &
             &  2, -3, &
             &  2,  1  ], [4, 2], order = [2, 1])
      x_exp = reshape([&
             &  1,  0, &
             &  1,  1, &
             &  1,  0, &
             &  1, -1  ], [4, 2], order = [2, 1])
      call M2%solve_mat(y, x)
      call tc%assert_eq(x_exp, x, rtol, atol, "solve_mat")

      call M2%solve_mat(y, x, trans = "T")
      call tc%assert_eq(x_exp, x, rtol, atol, "solve_mat T")
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

      call tc%assert_eq([1.0, 1.0],      Asp%a, rtol, atol, "delete block: a" )
      call tc%assert_eq([1, 1, 1, 2, 3], int(Asp%ia),       "delete block: ia")
      call tc%assert_eq([1, 2],          Asp%ja,            "delete block: ja")
    end block

    ! get_global_idx, get_local_idx
    block
      integer :: ib, jb, ig, jg, il, jl

      call M1%get_global_idx(2, 2, 1, 2, ig, jg)
      call tc%assert_eq(2, ig, "global row")
      call tc%assert_eq(3, jg, "global col")

      call M1%get_local_idx(3, 2, ib, jb, il, jl)
      call tc%assert_eq(2, ib, "block row")
      call tc%assert_eq(2, jb, "block col")
      call tc%assert_eq(2, il, "local row")
      call tc%assert_eq(1, jl, "local col")
    end block

    ! set
    block
      type(dense_real)          :: A
      type(dense_real), pointer :: A_p

      call A%init(real(reshape([1], [1, 1])))
      call M1%set(1, 1, A)
      call M1%get(1, 1, A_p)
      call tc%assert_eq(A%d,   A_p%d, 0.0, 0.0, "set")

      call M1%set(1, 1, A, fact = 2.0)
      call M1%get(1, 1, A_p)
      call tc%assert_eq(2*A%d, A_p%d, 0.0, 0.0, "set fact")
    end block

    block
      type(block_real)          :: A
      type(dense_real), pointer :: A_p, M1_p

      call M1%copy(A)
      call tc%assert_eq(M1%nbrows, A%nbrows, "copy: nbrows")
      call tc%assert_eq(M1%nbcols, A%nbcols, "copy: nbcols")
      call tc%assert_eq(M1%i0,     A%i0,     "copy: i0")
      call tc%assert_eq(M1%i1,     A%i1,     "copy: i1")
      call tc%assert_eq(M1%j0,     A%j0,     "copy: j0")
      call tc%assert_eq(M1%j1,     A%j1,     "copy: j1")

      call M1%get(1, 1, M1_p)
      call A%get( 1, 1, A_p)
      call tc%assert_eq(M1_p%d, A_p%d, 0.0, 0.0, "copy: val (1,1)")

      call M1%get(1, 2, M1_p)
      call A%get( 1, 2, A_p)
      call tc%assert_eq(M1_p%d, A_p%d, 0.0, 0.0, "copy: val (1,2)")

      call M1%get(2, 1, M1_p)
      call A%get( 2, 1, A_p)
      call tc%assert_eq(M1_p%d, A_p%d, 0.0, 0.0, "copy: val (2,1)")

      call M1%get(2, 2, M1_p)
      call A%get( 2, 2, A_p)
      call tc%assert_eq(M1_p%d, A_p%d, 0.0, 0.0, "copy: val (2,2)")
    end block

    ! set_ptr, destruct
    block
      type(dense_real) :: ext_A

      d0 = reshape([1,0,0,2,3,5,0,4,6], [3,3])
      call ext_A%init(d0)
      call M1%set_ptr(1, 1, ext_A)
      call M1%set_ptr(1, 2, ext_A)
      call M1%set_ptr(2, 1, ext_A)

      ! destruct should not delete ext_A
      call M1%destruct()

      call tc%assert_eq(3, ext_A%nrows, "set_ptr")
    end block

    call tc%finish()
  end subroutine

end submodule
