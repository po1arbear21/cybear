submodule(test_matrix_m) test_sparse_m

  use ieee_arithmetic
  use matrix_m
  use math_m,         only: PI

  implicit none

contains

  module subroutine test_sparse()
    type(test_case) :: tc

    print "(A)", "test_sparse"
    call tc%init("sparse")

    ! sparse builder: check if sparse builder created correct sparse matrix
    block
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

    ! sparse builder: set -> reset row
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

    ! sparse builder: set -> reset row -> set_row
    block
      integer              :: row, j0, j1
      real, allocatable    :: vals(:)
      type(dense_real)     :: d_exp, d
      type(sparse_real)    :: sp
      type(spbuild_real)   :: sb

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
      call d_exp%init(4)
      call d%init(4)

      call sb%add(1, 1, 1.0)
      call sb%add(1, 2, 2.0)
      call sb%add(2, 2, 4.0)
      call sb%add(3, 3, 1.0)
      call sb%add(4, 2, 1.0)
      call sb%add(4, 4, 5.0)

      call sb%save
      call sp%to_dense(d_exp)

      !
      ! test 1: reset row 1 -> set_row(1)
      !
      ! A =
      !  6     0     7     8
      !  0     4     0     0
      !  0     0     1     0
      !  0     1     0     5
      row  = 1
      vals = [6, 0, 7, 8]
      ! sb%keep_struct(row) = .false.
      call sb%reset_row(row)
      call sb%set_row(row, vals)
      call sb%save
      call sp%to_dense(d)
      d_exp%d(row,:) = vals

      call tc%assert_eq(d_exp%d,  d%d, 1e-12, "sparse_builder: set -> reset_row -> set_row: 1: full row")

      !
      ! test 2: reset row 3 -> set_row using column arguments
      !
      ! A =
      !  6     0     7     8
      !  0     4     0     0
      !  0     9    10     0
      !  0     1     0     5
      row  = 3
      j0   = 2
      j1   = 3
      vals = [9, 10]
      sb%keep_struct(row) = .false.     ! needed if resetting only slice of row and there was data before
      call sb%reset_row(row)
      call sb%set_row(row, vals, j0=j0, j1=j1)
      call sb%save
      call sp%to_dense(d)
      d_exp%d(row,j0:j1) = vals

      call tc%assert_eq(d_exp%d,  d%d, 1e-12, "sparse_builder: set -> reset_row -> set_row: 2: slice of row")
    end block

    ! mul_vec
    block
      integer           :: i, n
      type(sparse_real) :: sA
      real, allocatable :: x(:), y(:), y_exp(:)

      x= [ 1, -2, 3,  -4]
      allocate (y(4))

      ! test: y <- A*x
      y     = ieee_value(y, ieee_quiet_nan)   ! set y to nan: in case y(i) is not overwritten then nan will throw error while testing
      y_exp = [-3, -8, 3, -22]

      call get_test_matrix(sA)
      call sA%mul_vec(x, y)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec 1")

      ! test: y <- A*x -2y
      y     = [ 1,  0, 0,  -1]
      y_exp = [-5, -8, 3, -20]

      call sA%mul_vec(x, y, fact_y=-2.0)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec 2")

      ! test: y <- (A^t)*x
      y     = ieee_value(y, ieee_quiet_nan)   ! set y to nan: in case y(i) is not overwritten then nan will throw error while testing
      y_exp = [1, -10, 3, -20]

      call sA%mul_vec(x, y, trans='T')
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec 3")

      ! test: y <- (A^t)*x -2y
      y     = [ 1,   0, 0,  -1]
      y_exp = [-1, -10, 3, -18]

      call sA%mul_vec(x, y, fact_y=-2.0, trans='T')
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec 4")


      ! tests using large matrix
      ! expected solutions are too long to include here. just compare the sum
      n = 500
      call get_test_matrix3(500, sA)
      deallocate (x, y)
      x = sin([( (4*PI*i)/n, i=1,n )])    ! x(1)~0, x(n)=sin(4PI)=0. two full periods
      allocate (y(n))

      ! test: y <- A*x
      call sA%mul_vec(x, y)
      call tc%assert_eq(-6.418833765673056e+05, sum(y), 1e-7, "mul_vec 5")

      ! test: y <- A*x -3y
      y = [(i, i=1,n)]
      call sA%mul_vec(x, y, fact_y=-3.0)
      call tc%assert_eq(-1.017633376567306e+06, sum(y), 1e-6, "mul_vec 6")

      ! test: y <- (A^t)*x
      y = ieee_value(y, ieee_quiet_nan)   ! set y to nan: in case y(i) is not overwritten then nan will throw error while testing
      call sA%mul_vec(x, y, trans='T')
      call tc%assert_eq(1.601097642095075e+08,  sum(y), 1e-4, "mul_vec 7")

      ! test: y <- (A^t)*x +4y
      y = [(2*i, i=1,n)]
      call sA%mul_vec(x, y, fact_y=4.0, trans='T')
      call tc%assert_eq(1.611117642095075e+08,  sum(y), 1e-4, "mul_vec 8")
    end block

    ! mul_mat
    block
      type(sparse_real) :: sA
      real, allocatable :: x(:,:), y(:,:), y_exp(:,:)

      ! test 1
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
      y = ieee_value(y, ieee_quiet_nan)

      call get_test_matrix(sA)
      call sA%mul_mat(x, y)
      call tc%assert_eq(y_exp, y, 1e-12, "mul_mat 1")

      ! test 2: using trans, fact_y option for small matrices
      ! y <- transpose(A)*B -2*A = transpose(A)*x-2*y
      block
        type(dense_real)  :: dA, dB
        type(sparse_real) :: sB

        call get_test_matrix2(sB)
        call dA%init(4)
        call dB%init(4)
        call sA%to_dense(dA)
        call sB%to_dense(dB)

        ! y = transpose(A)*B -2*A =
        ! -2    -4     5     1
        ! -1     0    10     2
        !  0     0    -2    -3
        ! -5    -2     0   -10
        x     = dB%d
        y     = dA%d
        y_exp = reshape([-2,-1, 0,-5,-4, 0, 0,-2, 5,10,-2, 0, 1, 2,-3,-10], [4,4])
        call sA%mul_mat(x, y, fact_y=-2.0, trans='T')

        call tc%assert_eq(y_exp, y, 1e-12, "mul_mat 2")
      end block

      ! test 3: using trans, fact_y option for large matrices
      ! y <- transpose(A)*A -10^5 * A
      block
        type(dense_real)   :: dA
        integer, parameter :: n=500

        call get_test_matrix3(n, sA)
        call dA%init(n)
        call sA%to_dense(dA)

        x = dA%d
        y = dA%d
        call sA%mul_mat(x, y, fact_y=-1e5, trans='T')

        call tc%assert_eq(1.684838714750912e+15, sum(y), 1e2, "mul_mat 3")
      end block
    end block

    ! mul_sparse
    block
      type(sparse_real)    :: sA, sB, sC
      type(spbuild_real)   :: sbuild
      real,    allocatable :: a_exp(:)
      integer, allocatable :: ia_exp(:), ja_exp(:)

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

      call tc%assert_eq(a_exp,  sC%a, 1e-12, "mul sparse 1: a")
      call tc%assert_eq(ia_exp, sC%ia,       "mul sparse 1: ia")
      call tc%assert_eq(ja_exp, sC%ja,       "mul sparse 1: ja")

      ! test 2: large matrices. check sum of entries
      ! C <- A*A
      block
        integer, parameter :: n=500

        call get_test_matrix3(n, sA)
        call get_test_matrix3(n, sB)
        call sA%mul_sparse(sB, sC)

        call tc%assert_eq(8.913712864525120e+14, sum(sC%a), 1e2, "mul sparse 2")
      end block
    end block

    ! solve_vec
    block
      real, allocatable  :: b(:), x(:), x_exp(:)
      type(sparse_real)  :: sA

      b     = [1, 3, 0, -2]
      x_exp = [-0.5, 0.75, 0.0, -0.55]
      allocate (x(4))

      ! test 1: using default solver
      call get_test_matrix(sA)
      call sA%factorize
      call sA%solve_vec(b, x)
      call sA%destruct
      call tc%assert_eq(x_exp, x, 1e-12, "solve_vec: default solver")

      ! test 2: using pardiso solver
      call get_test_matrix(sA)
      sA%solver=SOLVER_PARDISO
      call sA%factorize
      call sA%solve_vec(b, x)
      call sA%destruct
      call tc%assert_eq(x_exp, x, 1e-12, "solve_vec: pardiso solver")

      ! test 3: using ilupack solver
      call get_test_matrix(sA)
      sA%solver=SOLVER_ILUPACK
      call sA%factorize
      call sA%solve_vec(b, x)
      call sA%destruct
      call tc%assert_eq(x_exp, x, 1e-12, "solve_vec: ilupack solver")

      ! test 4: using ilupack solver, changing ilupack parameters
      block
        use ilupack_m, only: get_ilupack_handle_ptr, ilupack_handle

        type(ilupack_handle), pointer :: ilu

        call get_test_matrix(sA)
        sA%solver=SOLVER_ILUPACK
        call sA%init_solver

        call get_ilupack_handle_ptr(sA%solver_handle, ilu)
        ilu%elbow    = 10
        ilu%droptol  = 3e-2
        ilu%droptolS = ilu%droptol / 10
        ilu%matching = 0
        ilu%restol   = 1e-14

        call sA%factorize
        call sA%solve_vec(b, x)
        call sA%destruct
        call tc%assert_eq(x_exp, x, 1e-12, "solve_vec: ilupack solver")
      end block
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
      call sA%destruct

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

    ! diag
    block
      type(sparse_real) :: S
      real              :: d_exp(4), d(4)

      call get_test_matrix(S)
      d_exp = [1,4,1,5]
      call S%diag(d)
      call tc%assert_eq(d_exp, d, 1e-12, "diag")

      call get_test_matrix2(S)
      d_exp = [0,2,0,0]
      call S%diag(d)
      call tc%assert_eq(d_exp, d, 1e-12, "diag")
    end block

    ! zero
    block
      integer, parameter :: n = 3
      type(sparse_real)  :: z_sp
      type(dense_real)   :: z_d
      real, allocatable  :: z_exp(:,:)

      z_sp = sparse_zero_real(n)

      call z_d%init(n)
      call z_sp%to_dense(z_d)

      allocate (z_exp(n,n), source=0.0)

      call tc%assert_eq(z_exp, z_d%d, epsilon(1.0), "zero(3)")
    end block

    ! nnz
    block
      type(sparse_real)  :: s
      type(spbuild_real) :: spb

      call get_test_matrix(s)
      call tc%assert_eq(6, s%nnz(), "nnz")

      call get_test_matrix(s)
      call tc%assert_eq(6, s%nnz(only_nonzeros=.true.), "nnz")

      call get_test_matrix2(s)
      call tc%assert_eq(5, s%nnz(), "nnz")

      call get_test_matrix2(s)
      call tc%assert_eq(5, s%nnz(only_nonzeros=.true.), "nnz")

      ! test matrix also includes 0 entries
      ! S = [ 0  1  2]
      !     [-3  0  4]
      !     [ 5 -6  0]
      call s%init(3)
      call spb%init(s)
      call spb%add(1, 2,  1.0)
      call spb%add(1, 3,  2.0)
      call spb%add(2, 1, -3.0)
      call spb%add(2, 2,  0.0)   ! 0 entry !!!
      call spb%add(2, 3,  4.0)
      call spb%add(3, 1,  5.0)
      call spb%add(3, 2, -6.0)
      call spb%add(3, 3, -0.0)   ! 0 entry !!!
      call spb%save

      call tc%assert_eq(6, s%nnz(), "nnz")
    end block

    ! to_band
    block
      type(band_real)   :: b
      type(sparse_real) :: s

      !
      ! test 1: set band equal sparse, not insertion as block anywhere
      !
      call get_test_matrix(s)
      call b%init(4, 2, nupper=1)
      call s%to_band(b)

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
      call get_test_matrix(s)
      call b%init(7, 3, nupper=0)
      b%d = -1                          ! just for testing purposes: to see which data will get overwritten
      call s%to_band(b, i0=3, j0=2)

      ! just for testing purposes
      block
        integer :: i, j
        type(dense_real)  :: d
        real, allocatable :: d_exp(:,:)

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
        call b%to_dense(d)
        call tc%assert_eq(d_exp, d%d, 1e-12, "to_band: insert sparse into band")
      end block
    end block

    call tc%finish
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

  subroutine get_test_matrix3(n, S)
    !! some large sparse matrix

    integer,           intent(in)  :: n
      !! matrix dimension
    type(sparse_real), intent(out) :: S

    integer :: i, j
    type(spbuild_real) :: sbuild

    ! S_ij = {i+2*j-i**2    if    (i+j)%31 == 0
    !        {0             else

    call S%init(n)
    call sbuild%init(S)

    do i = 1, n
      do j = 1, n
        if (mod(i+j, 31) == 0) call sbuild%add(i, j, real(i+2*j-i**2))
      end do
    end do

    call sbuild%save
  end subroutine

end submodule
