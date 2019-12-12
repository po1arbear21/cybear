submodule(test_matrix_m) test_sparse_m
  use matrix_m
  implicit none

contains
  module subroutine test_sparse()
    type(test_case)   :: tc
    type(sparse_real) :: sA

    ! s = A =
    !  1     2     0     0
    !  0     4     0     0
    !  0     0     1     0
    !  0     1     0     5
    ! will be set in first test (block)

    print "(1A)", "test_sparse"
    call tc%init("sparse")

    ! sparse builder: check if sparse builder created correct sparse matrix
    block
      integer            :: i, j
      type(spbuild_real) :: sbuild
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

      call sA%init(4)
      call sbuild%init(sA)

      call sbuild%add(1, 1, 1.0)
      call sbuild%add(1, 2, 2.0)
      call sbuild%add(2, 2, 4.0)
      call sbuild%add(3, 3, 1.0)
      call sbuild%add(4, 2, 1.0)
      call sbuild%add(4, 4, 5.0)

      call sbuild%save

      a_exp  = [1, 2, 4, 1, 1, 5]
      ia_exp = [1, 3, 4, 5, 7]
      ja_exp = [1, 2, 2, 3, 2, 4]

      call tc%assert_eq(a_exp,  sA%a, 1e-12, "sparse_builder to sparse: a")
      call tc%assert_eq(ia_exp, sA%ia,       "sparse_builder to sparse: ia")
      call tc%assert_eq(ja_exp, sA%ja,       "sparse_builder to sparse: ja")
    end block

    ! add sparse: B <- B + fact * A
    block
      type(sparse_real)    :: sB
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

      ! B-2*A=
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
      call sB%add_sparse(sA, fact=-2.0)

      call tc%assert_eq(a_exp,  sB%a, 1e-12, "add_sparse: a")
      call tc%assert_eq(ia_exp, sB%ia,       "add_sparse: ia")
      call tc%assert_eq(ja_exp, sB%ja,       "add_sparse: ja")
    end block

    ! add_sparse3: C <- fact1 * A + fact2 * B
    block
      type(sparse_real)    :: sB, sC
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
      call sB%add_sparse(sA, sC, fact1=3.0, fact2=-2.0)

      call tc%assert_eq(a_exp,  sC%a, 1e-12, "add_sparse: a")
      call tc%assert_eq(ia_exp, sC%ia,       "add_sparse: ia")
      call tc%assert_eq(ja_exp, sC%ja,       "add_sparse: ja")
    end block

    ! mat vec
    block
      real, allocatable, dimension(:) :: x, y, y_exp

      x     = [1,-2,3,-4]
      y_exp = [-3, -8, 3, -22]
      allocate(y(4))

      call sA%mul_vec(x, y)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec")
    end block

    ! mul mat
    block
      real, allocatable, dimension(:,:) :: x, y, y_exp

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

      call sA%mul_mat(x, y)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_mat")
    end block

    ! mul_sparse
    block
      type(sparse_real)    :: sB, sC
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

      call sA%mul_sparse(sB, sC)

      call tc%assert_eq(a_exp,  sC%a, 1e-12, "mul sparse: a")
      call tc%assert_eq(ia_exp, sC%ia,       "mul sparse: ia")
      call tc%assert_eq(ja_exp, sC%ja,       "mul sparse: ja")
    end block

    ! solve vec
    block
      real, allocatable, dimension(:) :: b, x, x_exp
      type(sparse_real)               :: s2

      b     = [1, 3, 0, -2]
      x_exp = [-0.5, 0.75, 0.0, -0.55]
      allocate(x(4))

      s2 = sA
      call s2%factorize
      call s2%solve_vec(b, x)

      call tc%assert_eq(x_exp, x, 1e-12, "solve_vec")
    end block

    ! solve mat
    block
      real, allocatable, dimension(:,:) :: b, x, x_exp
      type(sparse_real)                 :: s2

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

      s2 = sA
      call s2%factorize
      call s2%solve_mat(b, x)

      call tc%assert_eq(x_exp, x, 1e-12, "solve_mat")
    end block


    call tc%finish()
  end subroutine
end submodule
