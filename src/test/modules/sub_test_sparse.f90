submodule(test_matrix_m) test_sparse_m
  use matrix_m
  implicit none

contains

  module subroutine test_sparse()
    type(test_case)   :: tc

    print "(1A)", "test_sparse"
    call tc%init("sparse")

    ! sparse builder: check if sparse builder created correct sparse matrix
    block
      integer              :: i, j
      type(sparse_real)    :: sA
      real,    allocatable :: a_exp(:)
      integer, allocatable :: ia_exp(:), ja_exp(:)


      call get_test_matrix(sA)

      a_exp  = [1, 2, 4, 1, 1, 5]
      ia_exp = [1, 3, 4, 5, 7]
      ja_exp = [1, 2, 2, 3, 2, 4]

      call tc%assert_eq(a_exp,  sA%a, 1e-12, "sparse_builder to sparse: a")
      call tc%assert_eq(ia_exp, sA%ia,       "sparse_builder to sparse: ia")
      call tc%assert_eq(ja_exp, sA%ja,       "sparse_builder to sparse: ja")
    end block

    ! sparse builder: reset row + set
    block
      type(spbuild_real)   :: sb
      type(sparse_real)    :: sp
      real,    allocatable :: a_exp(:)
      integer, allocatable :: ia_exp(:), ja_exp(:)

      ! A =
      !  1     2     0     0
      !  0     4     0     0
      !  0     0     1     0
      !  0     1     0     5
      !
      ! a  = [1 2 4 1 1 5]
      ! ia = [1 3 4 5 7]
      ! ja = [1 2 2 3 2 4]

      call sp%init(4)
      call sb%init(sp)

      call sb%add(1, 1, 1.0)
      call sb%add(1, 2, 2.0)
      call sb%add(2, 2, 4.0)
      call sb%add(3, 3, 1.0)
      call sb%add(4, 2, 1.0)
      call sb%add(4, 4, 5.0)

      !
      ! test 1: reset row 1
      !
      ! a  = [4 1 1 5]
      ! ia = [1 1 2 3 5]
      ! ja = [2 3 2 4]

      call sb%reset_row(1)

      call sb%save

      a_exp  = [4, 1, 1, 5]
      ia_exp = [1, 1, 2, 3, 5]
      ja_exp = [2, 3, 2, 4]

      call tc%assert_eq(a_exp,  sp%a, 1e-12, "sparse_builder%resest_row to sparse: a")
      call tc%assert_eq(ia_exp, sp%ia,       "sparse_builder%resest_row to sparse: ia")
      call tc%assert_eq(ja_exp, sp%ja,       "sparse_builder%resest_row to sparse: ja")

      !
      ! test 2: insert 1 on diagonal
      !
      ! a  = [1 4 1 1 5]
      ! ia = [1 2 3 4 6]
      ! ja = [1 2 3 2 4]
      call sb%set(1, 1, 1.0)

      call sb%save

      a_exp  = [1, 4, 1, 1, 5]
      ia_exp = [1, 2, 3, 4, 6]
      ja_exp = [1, 2, 3, 2, 4]

      call tc%assert_eq(a_exp,  sp%a, 1e-12, "sparse_builder%set to sparse: a")
      call tc%assert_eq(ia_exp, sp%ia,       "sparse_builder%set to sparse: ia")
      call tc%assert_eq(ja_exp, sp%ja,       "sparse_builder%set to sparse: ja")
    end block

    ! mul_vec
    block
      type(sparse_real) :: sA
      real, allocatable :: x(:), y(:), y_exp(:)

      x     = [1,-2,3,-4]
      y_exp = [-3, -8, 3, -22]
      allocate(y(4))

      ! set y to nan
      y = 0.0 / 0.0

      call get_test_matrix(sA)
      call sA%mul_vec(x, y)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec")
    end block

    ! mul_mat
    block
      type(sparse_real) :: sA
      real, allocatable :: x(:,:), y(:,:), y_exp(:,:)

      ! x = A(:,1:3) =
      !  1     2     0
      !  0     4     0
      !  0     0     1
      !  0     1     0
      x = reshape([1, 0, 0, 0, 2, 4, 0, 1, 0, 0, 1, 0], [4,3])

      ! y = A*A(:,1:3) =
      !  1    10     0
      !  0    16     0
      !  0     0     1
      !  0     9     0
      y_exp = reshape([1, 0, 0, 0, 10, 16, 0, 9, 0, 0, 1, 0], [4,3])
      allocate(y(4,3))

      ! set y to nan
      y = 0.0 / 0.0

      call get_test_matrix(sA)
      call sA%mul_mat(x, y)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_mat")
    end block

    ! mul_sparse
    block
      type(sparse_real)    :: sA, sB, sC
      type(spbuild_real)   :: sbuild
      real,    allocatable :: a_exp(:)
      integer, allocatable :: ia_exp(:), ja_exp(:)

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

      call sbuild%save

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

      call get_test_matrix(sA)
      call sA%mul_sparse(sB, sC)

      call tc%assert_eq(a_exp,  sC%a, 1e-12, "mul sparse: a")
      call tc%assert_eq(ia_exp, sC%ia,       "mul sparse: ia")
      call tc%assert_eq(ja_exp, sC%ja,       "mul sparse: ja")
    end block

    ! solve_vec
    block
      real, allocatable :: b(:), x(:), x_exp(:)
      type(sparse_real) :: sA

      b     = [1, 3, 0, -2]
      x_exp = [-0.5, 0.75, 0.0, -0.55]
      allocate(x(4))

      call get_test_matrix(sA)
      call sA%factorize
      call sA%solve_vec(b, x)

      call tc%assert_eq(x_exp, x, 1e-12, "solve_vec")
    end block

    ! solve_mat
    block
      real, allocatable :: b(:,:), x(:,:), x_exp(:,:)
      type(sparse_real) :: sA

      ! b =
      !   1    -2     3
      !  -4     5    -6
      !   7    -8     9
      ! -10    11   -12
      !
      ! A\b =
      !   3.0000   -4.5000    6.0000
      !  -1.0000    1.2500   -1.5000
      !   7.0000   -8.0000    9.0000
      !  -1.8000    1.9500   -2.1000

      b     = reshape([1,-4,7,-10,-2,5,-8,11,3,-6,9,-12], [4,3])
      x_exp = reshape([3.0, -1.0, 7.0, -1.8, -4.5, 1.25, -8.0, 1.95, 6.0, -1.5, 9.0, -2.1], [4,3])
      allocate(x(4,3))

      call get_test_matrix(sA)
      call sA%factorize
      call sA%solve_mat(b, x)

      call tc%assert_eq(x_exp, x, 1e-12, "solve_mat")
    end block

    ! add_sparse: A <- A + fact * B
    block
      type(sparse_real)    :: A, B
      real,    allocatable :: a_exp(:)
      integer, allocatable :: ia_exp(:), ja_exp(:)

      call get_test_matrix2(A)
      call get_test_matrix(B)

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
      call A%add_sparse(B, fact=-2.0)

      call tc%assert_eq(a_exp,  A%a, 1e-12, "add_sparse: a")
      call tc%assert_eq(ia_exp, A%ia,       "add_sparse: ia")
      call tc%assert_eq(ja_exp, A%ja,       "add_sparse: ja")
    end block

    ! add_sparse3: C <- fact1 * A + fact2 * B
    block
      type(sparse_real)    :: A, B, C
      real,    allocatable :: a_exp(:)
      integer, allocatable :: ia_exp(:), ja_exp(:)

      call get_test_matrix(A)
      call get_test_matrix2(B)

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
      call B%add_sparse(A, C, fact1=3.0, fact2=-2.0)

      call tc%assert_eq(a_exp,  C%a, 1e-12, "add_sparse: a")
      call tc%assert_eq(ia_exp, C%ia,       "add_sparse: ia")
      call tc%assert_eq(ja_exp, C%ja,       "add_sparse: ja")
    end block

    ! add_band: S <- S + fact * B
    block
      type(sparse_real)    :: S
      type(band_real)      :: B
      real,    allocatable :: a_exp(:)
      integer, allocatable :: ia_exp(:), ja_exp(:)

      call get_test_matrix(S)

      ! B =
      !  3    -2     0     0
      ! -1     3    -2     0
      !  0    -1     3    -2
      !  0     0    -1     3
      call B%init(4, 1, 1)
      B%d(-1,2:4) = -2.0
      B%d( 0, : ) = 3.0
      B%d(+1,1:3) = -1.0

      ! S+2*B=
      !  7    -2     0     0
      ! -2    10    -4     0
      !  0    -2     7    -4
      !  0     1    -2    11
      a_exp  = [7, -2, -2, 10, -4, -2, 7, -4, 1, -2, 11]
      ia_exp = [1, 3,  6,  9, 12]
      ja_exp = [1, 2, 1, 2, 3, 2, 3, 4, 2, 3, 4]
      call S%add_band(B, 2.0)

      call tc%assert_eq(a_exp,  S%a, 1e-12, "add_band: a")
      call tc%assert_eq(ia_exp, S%ia,       "add_band: ia")
      call tc%assert_eq(ja_exp, S%ja,       "add_band: ja")
    end block

    ! add_band3: S2 <- fact1 * S + fact2 * B
    block
      type(sparse_real)    :: S, S2
      type(band_real)      :: B
      real,    allocatable :: a_exp(:)
      integer, allocatable :: ia_exp(:), ja_exp(:)

      call get_test_matrix(S)

      ! B =
      !  3    -2     0     0
      ! -1     3    -2     0
      !  0    -1     3    -2
      !  0     0    -1     3
      call B%init(4, 1, 1)
      B%d(-1,2:4) = -2.0
      B%d( 0, : ) = 3.0
      B%d(+1,1:3) = -1.0

      ! -1*S + 2*B=
      !  5    -6     0     0
      ! -2     2    -4     0
      !  0    -2     5    -4
      !  0    -1    -2     1
      a_exp  = [5, -6, -2, 2, -4, -2, 5, -4, -1, -2, 1]
      ia_exp = [1, 3,  6,  9, 12]
      ja_exp = [1, 2, 1, 2, 3, 2, 3, 4, 2, 3, 4]
      call S%add_band(B, S2, fact1 = -1.0, fact2 = 2.0)

      call tc%assert_eq(a_exp,  S2%a, 1e-12, "add_band3: a")
      call tc%assert_eq(ia_exp, S2%ia,       "add_band3: ia")
      call tc%assert_eq(ja_exp, S2%ja,       "add_band3: ja")
    end block

    call tc%finish()
  end subroutine

  subroutine get_test_matrix(S)
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

    call sbuild%save
  end subroutine

  subroutine get_test_matrix2(S)
    type(sparse_real), intent(out) :: S

    type(spbuild_real) :: sbuild

    ! S =
    !  0     0     5     1
    !  0     2     0     0
    !  0     0     0    -3
    ! -1     0     0     0

    call S%init(4)
    call sbuild%init(S)

    call sbuild%add(1, 3,  5.0)
    call sbuild%add(1, 4,  1.0)
    call sbuild%add(2, 2,  2.0)
    call sbuild%add(3, 4, -3.0)
    call sbuild%add(4, 1, -1.0)

    call sbuild%save
  end subroutine

end submodule
