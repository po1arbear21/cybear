m4_include(../../util/macro.f90.inc)

submodule (test_matrix_m) test_sparse_m

  use ieee_arithmetic
  use math_m,         only: PI

  implicit none

contains

  module subroutine test_sparse()
    type(test_case) :: tc

    call tc%init("sparse")

    ! sparse builder: check if sparse builder created correct sparse matrix
    block
      integer, allocatable :: ia_exp(:), ja_exp(:)
      real,    allocatable :: a_exp(:)
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

      call sb%add(1, 1, 1.0)
      call sb%add(1, 2, 2.0)
      call sb%add(2, 2, 4.0)
      call sb%add(3, 3, 1.0)
      call sb%add(4, 2, 1.0)
      call sb%add(4, 4, 5.0)
      call sb%add(1, 2, 3.0)

      call sb%save()

      a_exp  = [1, 5, 4, 1, 1, 5]
      ia_exp = [1, 3, 4, 5, 7]
      ja_exp = [1, 2, 2, 3, 2, 4]

      call tc%assert_eq(a_exp,  sp%a, 0.0,  "sparse_builder to sparse: a")
      call tc%assert_eq(ia_exp, int(sp%ia), "sparse_builder to sparse: ia")
      call tc%assert_eq(ja_exp, sp%ja,      "sparse_builder to sparse: ja")
      deallocate (a_exp)
    end block

    ! sparse builder: check if sparse builder created correct sparse matrix
    block
      integer, allocatable :: ia_exp(:), ja_exp(:)
      real,    allocatable :: a_exp(:)
      type(sparse_real)    :: sA

      call example_matrix3(sA)

      a_exp  = [1, 2, 4, 1, 1, 5]
      ia_exp = [1, 3, 4, 5, 7]
      ja_exp = [1, 2, 2, 3, 2, 4]

      call tc%assert_eq(a_exp,  sA%a,  0.0, "sparse_builder to sparse: a")
      call tc%assert_eq(ia_exp, int(sA%ia), "sparse_builder to sparse: ia")
      call tc%assert_eq(ja_exp, sA%ja,      "sparse_builder to sparse: ja")
    end block

    ! sparse builder: set -> reset row
    block
      integer, allocatable :: ia_exp(:), ja_exp(:)
      real,    allocatable :: a_exp(:)
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

      call sb%save()

      a_exp  = [4, 1, 1, 5]
      ia_exp = [1, 1, 2, 3, 5]
      ja_exp = [2, 3, 2, 4]

      call tc%assert_eq(a_exp,  sp%a, 0.0,  "sparse_builder%resest_row to sparse: a")
      call tc%assert_eq(ia_exp, int(sp%ia), "sparse_builder%resest_row to sparse: ia")
      call tc%assert_eq(ja_exp, sp%ja,      "sparse_builder%resest_row to sparse: ja")

      !
      ! test 2: insert 1 on diagonal
      !
      ! a  = [1 4 1 1 5]
      ! ia = [1 2 3 4 6]
      ! ja = [1 2 3 2 4]
      call sb%set(1, 1, 1.0)

      call sb%save()

      a_exp  = [1, 4, 1, 1, 5]
      ia_exp = [1, 2, 3, 4, 6]
      ja_exp = [1, 2, 3, 2, 4]

      call tc%assert_eq(a_exp,  sp%a, 0.0,  "sparse_builder%set to sparse: a")
      call tc%assert_eq(ia_exp, int(sp%ia), "sparse_builder%set to sparse: ia")
      call tc%assert_eq(ja_exp, sp%ja,      "sparse_builder%set to sparse: ja")
    end block

    ! sparse builder: set -> reset row -> set_row
    block
      integer            :: row, j0, j1
      real, allocatable  :: vals(:)
      type(dense_real)   :: d_exp, d
      type(sparse_real)  :: sp
      type(spbuild_real) :: sb

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

      call sb%save()
      call matrix_convert(sp, d_exp)

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
      call sb%save()
      call matrix_convert(sp, d)
      d_exp%d(row,:) = vals

      call tc%assert_eq(d_exp%d,  d%d, 0.0, "sparse_builder: set -> reset_row -> set_row: 1: full row")

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
      call sb%set_row(row, vals, j0 = j0, j1 = j1)
      call sb%save()
      call matrix_convert(sp, d)
      d_exp%d(row,j0:j1) = vals

      call tc%assert_eq(d_exp%d, d%d, 0.0, "sparse_builder: set -> reset_row -> set_row: 2: slice of row")
    end block

    ! sparse builder: load
    block
      integer            :: ia_exp(5), ja_exp(6)
      real               :: a_exp(6)
      type(sparse_real)  :: sA
      type(spbuild_real) :: sb

      call example_matrix3(sA)
      call sb%init(sA)
      call sb%load()
      call sb%save()

      a_exp  = [1, 2, 4, 1, 1, 5]
      ja_exp = [1, 2, 2, 3, 2, 4]
      ia_exp = [1, 3, 4, 5, 7]

      call tc%assert_eq(a_exp,  sA%a,  0.0, "load: a")
      call tc%assert_eq(ja_exp, sA%ja,      "load: ja")
      call tc%assert_eq(ia_exp, int(sA%ia), "load: ia")
    end block

    ! scale
    block
      real, parameter   :: l = 3.0
      real              :: val(6)
      type(sparse_real) :: sA

      call example_matrix3(sA)
      call sA%scale(l)

      val = [1, 2, 4, 1, 1, 5]

      call tc%assert_eq(val * l, sA%a, 0.0, "scale")
    end block

    ! mul_vec real
    block
      integer           :: i, n
      real, allocatable :: x(:), y(:), y_exp(:)
      type(sparse_real) :: sA, sE

      x= [1, -2, 3, -4]
      allocate (y(4))

      ! get matrices
      call get_empty_matrix(sE)
      call example_matrix3(sA)

      ! test: y <- A*x
      y     = ieee_value(y, ieee_quiet_nan)
      y_exp = [0, 0, 0, 0]
      call sE%mul_vec(x, y)
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec real empty 1")
      y     = ieee_value(y, ieee_quiet_nan)
      y_exp = [-3, -8, 3, -22]
      call sA%mul_vec(x, y)
      call tc%assert_eq(y_exp, y, 1e-15, "mul_vec real 1")

      ! test: y <- A*x -2y
      y     = [ 1, 0, 0, -1]
      y_exp = [-2, 0, 0,  2]
      call sE%mul_vec(x, y, fact_y = -2.0)
      call tc%assert_eq(y_exp, y, 1e-15, "mul_vec real empty 2")
      y     = [ 1,  0, 0, - 1]
      y_exp = [-5, -8, 3, -20]
      call sA%mul_vec(x, y, fact_y = -2.0)
      call tc%assert_eq(y_exp, y, 1e-15, "mul_vec real 2")

      ! test: y <- (A^t)*x
      y     = ieee_value(y, ieee_quiet_nan)
      y_exp = [0, 0, 0, 0]
      call sE%mul_vec(x, y, trans = 'T')
      call tc%assert_eq(y_exp, y, 1e-15, "mul_vec real empty 3")
      y     = ieee_value(y, ieee_quiet_nan)
      y_exp = [1, -10, 3, -20]
      call sA%mul_vec(x, y, trans = 'T')
      call tc%assert_eq(y_exp, y, 1e-15, "mul_vec real 3")

      ! test: y <- (A^t)*x -2y
      y     = [ 1, 0, 0, -1]
      y_exp = [-2, 0, 0,  2]
      call sE%mul_vec(x, y, fact_y = -2.0, trans = 'T')
      call tc%assert_eq(y_exp, y, 1e-12, "mul_vec real empty 4")
      y     = [ 1,   0, 0,  -1]
      y_exp = [-1, -10, 3, -18]
      call sA%mul_vec(x, y, fact_y = -2.0, trans = 'T')
      call tc%assert_eq(y_exp, y, 1e-15, "mul_vec real 4")

      ! tests using large matrix
      ! expected solutions are too long to include here. just compare the sum
      n = 500
      call get_test_matrix3(500, sA)
      deallocate (x, y)
      x = sin([( (4*PI*i)/n, i = 1, n )])    ! x(1)~0, x(n)=sin(4PI)=0. two full periods
      allocate (y(n))

      ! test: y <- A*x
      call sA%mul_vec(x, y)
      call tc%assert_eq(-6.418833765673056e+05, sum(y), 1e-7, "mul_vec real 5")

      ! test: y <- A*x -3y
      y = [(i, i=1,n)]
      call sA%mul_vec(x, y, fact_y = -3.0)
      call tc%assert_eq(-1.017633376567306e+06, sum(y), 1e-6, "mul_vec real 6")

      ! test: y <- (A^t)*x
      y = ieee_value(y, ieee_quiet_nan)   ! set y to nan: in case y(i) is not overwritten then nan will throw error while testing
      call sA%mul_vec(x, y, trans = 'T')
      call tc%assert_eq( 1.601097642095075e+08, sum(y), 1e-4, "mul_vec real 7")

      ! test: y <- (A^t)*x +4y
      y = [(2*i, i=1,n)]
      call sA%mul_vec(x, y, fact_y = 4.0, trans = 'T')
      call tc%assert_eq( 1.611117642095075e+08, sum(y), 1e-4, "mul_vec real 8")
    end block

    ! mul_vec cmplx
    block
      complex, allocatable :: x(:), y(:), y_exp(:)
      integer              :: i, n
      type(sparse_cmplx)   :: sA, sE

      x = [(1, 1), (-2, 0), (3, -1), (-4, 2)]
      allocate (y(4))

      ! get matrices
      call sE%init(4)           ! empty matrix
      call example_matrix5(sA)

      ! test: y <- A*x
      y     = ieee_value(1.0, ieee_quiet_nan)
      y_exp = [0, 0, 0, 0]
      call sE%mul_vec(x, y)
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec cmplx empty 1")
      y     = ieee_value(1.0, ieee_quiet_nan)
      y_exp = [(-9, -5), (-4, 0), (-1, -3), (-2, -4)]
      call sA%mul_vec(x, y)
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec cmplx 1")

      ! test: y <- A*x -2y
      y     = [( 2, 0), (0,  1), ( 1,  1), ( 1, -2)]
      y_exp = [(-4, 0), (0, -2), (-2, -2), (-2,  4)]
      call sE%mul_vec(x, y, fact_y = (-2, 0))
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec cmplx empty 2")
      y     = [( 2,   0), ( 0,  1), ( 1,  1), ( 1, -2)]
      y_exp = [(-13, -5), (-4, -2), (-3, -5), (-4,  0)]
      call sA%mul_vec(x, y, fact_y = (-2, 0))
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec cmplx 2")

      ! test: y <- (A^t)*x
      y     = ieee_value(1.0, ieee_quiet_nan)
      y_exp = [0, 0, 0, 0]
      call sE%mul_vec(x, y, trans = 'T')
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec cmplx empty 3")
      y     = ieee_value(1.0, ieee_quiet_nan)
      y_exp = [(0, 2), (-4, 0), (0, -8), (-2, -4)]
      call sA%mul_vec(x, y, trans = 'T')
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec cmplx 3")

      ! test: y <- (A^t)*x -2y
      y     = [( 2, 0), (0,  1), ( 1,  1), ( 1, -2)]
      y_exp = [(-4, 0), (0, -2), (-2, -2), (-2,  4)]
      call sE%mul_vec(x, y, fact_y = (-2, 0), trans = 'T')
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec cmplx empty 4")
      y     = [( 2, 0), ( 0,  1), ( 1,   1), ( 1, -2)]
      y_exp = [(-4, 2), (-4, -2), (-2, -10), (-4,  0)]
      call sA%mul_vec(x, y, fact_y = (-2, 0), trans = 'T')
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec cmplx 4")

      ! tests using large matrix
      ! expected solutions are too long to include here. just compare the sum
      n = 600
      call example_matrix4(sA)
      deallocate (x, y)
      x = [(cmplx(i, n-2*i), i=0, n-1)]    ! [(0, 600), (1, 598)...]
      allocate (y(n))
      ! test: y <- A*x
      call sA%mul_vec(x, y)
      call tc%assert_eq((46419.56269819, 690392.61484622), sum(y), 1e-8, "mul_vec cmplx 5")

      ! test: y <- A*x -3y
      y = [(cmplx(i,-i), i=0, n-1)]
      call sA%mul_vec(x, y, fact_y = (3, 0))
      call tc%assert_eq((585519.56269819, 151292.61484622), sum(y), 1e-8, "mul_vec cmplx 6")

      ! test: y <- (A^t)*x
      y = ieee_value(1.0, ieee_quiet_nan)   ! set y to nan: in case y(i) is not overwritten then nan will throw error while testing
      call sA%mul_vec(x, y, trans='T')
      call tc%assert_eq((46203.56269819, 690284.61484622), sum(y), 1e-8, "mul_vec cmplx 7")

      ! test: y <- (A^t)*x +4y
      y = [(cmplx(i,-i), i=0, n-1)]
      call sA%mul_vec(x, y, fact_y = (-1, 1), trans='T')
      call tc%assert_eq((46203.56269819, 1049684.61484622), sum(y), 1e-8, "mul_vec cmplx 8")
    end block

    ! mul_vec_slice
    block
      real, allocatable :: x(:), y(:), y_exp(:)
      type(sparse_real) :: sA

      x = [1, -2, 3, -4]
      allocate (y(2))

      ! get matrix
      call example_matrix3(sA)

      ! testing i0, i1. y <- A*x-2y
      y_exp = [-5, -8]
      y     = [ 1,  0]
      call sA%mul_vec_slice(x, y, fact_y = -2.0, i1 = 2)
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec_slice")
      y     = [1, 0]
      call sA%mul_vec_slice(x, y, fact_y = -2.0, i0 = 1, i1 = 2)
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec_slice")

      y_exp = [-8, 3]
      y     = [ 0, 0]
      call sA%mul_vec_slice(x, y, fact_y = -2.0, i0 = 2, i1 = 3)
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec_slice")

      y_exp = [3, -20]
      y     = [0, - 1]
      call sA%mul_vec_slice(x, y, fact_y = -2.0, i0 = 3, i1 = 4)
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec_slice")
      y     = [0, -1]
      call sA%mul_vec_slice(x, y, fact_y = -2.0, i0 = 3)
      call tc%assert_eq(y_exp, y, 0.0, "mul_vec_slice")
    end block

    ! mul_mat
    block
      real, allocatable :: x(:,:), y(:,:), y_exp(:,:)
      type(sparse_real) :: sA, sE

      ! test 1
      ! x = A(:,1:3) =
      !  1     2     0
      !  0     4     0
      !  0     0     1
      !  0     1     0
      x = reshape([1, 0, 0, 0, 2, 4, 0, 1, 0, 0, 1, 0], [4,3])

      ! set y to nan
      allocate(y(4,3))
      y = ieee_value(y, ieee_quiet_nan)

      allocate (y_exp(4,3), source = 0.0)
      call get_empty_matrix(sE)
      call sE%mul_mat(x, y)
      call tc%assert_eq(y_exp, y, 0.0, "mul_mat empty")

      ! y = A*A(:,1:3) =
      !  1    10     0
      !  0    16     0
      !  0     0     1
      !  0     9     0
      y_exp = reshape([1, 0, 0, 0, 10, 16, 0, 9, 0, 0, 1, 0], [4,3])

      call example_matrix3(sA)
      call sA%mul_mat(x, y)
      call tc%assert_eq(y_exp, y, 0.0, "mul_mat 1")

      ! test 2: using trans, fact_y option for small matrices
      ! y <- transpose(A)*B -2*A = transpose(A)*x-2*y
      block
        type(dense_real)  :: dA, dB
        type(sparse_real) :: sB

        call get_test_matrix2(sB)
        call dA%init(4)
        call dB%init(4)
        call matrix_convert(sA, dA)
        call matrix_convert(sB, dB)

        ! y = transpose(A)*B -2*A =
        ! -2    -4     5     1
        ! -1     0    10     2
        !  0     0    -2    -3
        ! -5    -2     0   -10
        x     = dB%d
        y     = dA%d
        y_exp = reshape([-2,-1, 0,-5,-4, 0, 0,-2, 5,10,-2, 0, 1, 2,-3,-10], [4,4])
        call sA%mul_mat(x, y, fact_y = -2.0, trans = 'T')

        call tc%assert_eq(y_exp, y, 0.0, "mul_mat 2")
      end block

      ! test 3: using trans, fact_y option for large matrices
      ! y <- transpose(A)*A -10^5 * A
      block
        type(dense_real)   :: dA
        integer, parameter :: n=500

        call get_test_matrix3(n, sA)
        call dA%init(n)
        call matrix_convert(sA, dA)

        x = dA%d
        y = dA%d
        call sA%mul_mat(x, y, fact_y = -1e5, trans = 'T')

        call tc%assert_eq(1.684838714750912e+15, sum(y), 1e2, "mul_mat 3")
      end block
    end block

    ! mul_sparse
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
      call sA%mul_sparse(sE, sC)
      call tc%assert(sC%is_empty(), "mul sparse empty 1")

      ! C = E*A
      call sE%mul_sparse(sA, sC)
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

      call sA%mul_sparse(sB, sC)

      call tc%assert_eq(a_exp,  sC%a, 0.0,  "mul sparse 1: a")
      call tc%assert_eq(ia_exp, int(sC%ia), "mul sparse 1: ia")
      call tc%assert_eq(ja_exp, sC%ja,      "mul sparse 1: ja")

      ! test 2: large matrices. check sum of entries
      ! C <- A*A
      block
        integer, parameter :: n=500

        call get_test_matrix3(n, sA)
        call get_test_matrix3(n, sB)
        call sA%mul_sparse(sB, sC)

        call tc%assert_eq(8.913712864525120e+14, sum(sC%a), 1e2, "mul sparse 2")
      end block

      ! test 3: large matrices in extern files
      ! C <- A*B
      call sA%input(file = "src/test/S1.test")
      call sB%input(file = "src/test/S2.test")
      call sE%input(file = "src/test/S3.test") ! expected S1 * S2
      call sA%mul_sparse(sB, sC)
      call tc%assert_eq(int(sE%ia), int(sC%ia), "mul sparse 3: ia")
      call tc%assert_eq(sE%ja, sC%ja,           "mul sparse 3: ja")
      call tc%assert_eq(sE%a,  sC%a, 1e-14,     "mul sparse 3: a")
    end block

    ! solve_vec real
    block
      real              :: b(4), x(4), x_exp(4)
      type(sparse_real) :: sA

      b     = [1, 3, 0, -2]
      x_exp = [-0.5, 0.75, 0.0, -0.55]

      ! test 1: using default solver
      call example_matrix3(sA)
      call sA%factorize()
      call sA%solve_vec(b, x)
      call tc%assert_eq(x_exp, x, 1e-15, "solve_vec real: default solver")

      ! test 2: using pardiso solver
      call example_matrix3(sA)
      sA%solver=SOLVER_PARDISO
      call sA%factorize()
      call sA%solve_vec(b, x)
      call tc%assert_eq(x_exp, x, 1e-15, "solve_vec real: pardiso solver")

      m4_divert(m4_ifdef({m4_ilupack},0,-1))
        ! test 3: using ilupack solver
        call example_matrix3(sA)
        sA%solver=SOLVER_ILUPACK
        call sA%factorize()
        call sA%solve_vec(b, x)
        call tc%assert_eq(x_exp, x, 1e-15, "solve_vec real: ilupack solver")

        ! test 4: using ilupack solver, changing ilupack parameters
        block
          use ilupack_m, only: get_ilupack_handle_ptr, ilupack_handle

          type(ilupack_handle), pointer :: ilu

          call example_matrix3(sA)
          sA%solver=SOLVER_ILUPACK
          call sA%init_solver()

          call get_ilupack_handle_ptr(sA%solver_handle, ilu)
          ilu%elbow    = 10
          ilu%droptol  = 3e-2
          ilu%droptolS = ilu%droptol / 10
          ilu%matching = 0
          ilu%restol   = 1e-14

          call sA%factorize()
          call sA%solve_vec(b, x)
          call tc%assert_eq(x_exp, x, 1e-15, "solve_vec real: ilupack solver")
        end block
      m4_divert(0)
    end block

    ! solve_vec cmplx
    block
      complex            :: b(4), x(4), x_exp(4)
      type(sparse_cmplx) :: sA

      b     = [(-1, -2), (2, 0), (0, -1), (0, 1)]
      x_exp = [  1,       1,      1,       1]

      ! test 1: using default solver
      call example_matrix5(sA)
      call sA%factorize()
      call sA%solve_vec(b, x)
      call tc%assert_eq(x_exp, x, 1e-15, "solve_vec cmplx: default solver")

      ! test 2: using pardiso solver
      call example_matrix5(sA)
      sA%solver=SOLVER_PARDISO
      call sA%factorize()
      call sA%solve_vec(b, x)
      call tc%assert_eq(x_exp, x, 1e-15, "solve_vec cmplx: pardiso solver")

      m4_divert(m4_ifdef({m4_ilupack},0,-1))
        ! test 3: using ilupack solver
        call example_matrix5(sA)
        sA%solver=SOLVER_ILUPACK
        call sA%factorize()
        call sA%solve_vec(b, x)
        call tc%assert_eq(x_exp, x, 1e-15, "solve_vec cmplx: ilupack solver")

        ! test 4: using ilupack solver, changing ilupack parameters
        block
          use ilupack_m, only: get_ilupack_handle_ptr, ilupack_handle

          type(ilupack_handle), pointer :: ilu

          call example_matrix5(sA)
          sA%solver=SOLVER_ILUPACK
          call sA%init_solver()

          call get_ilupack_handle_ptr(sA%solver_handle, ilu)
          ilu%elbow    = 10
          ilu%droptol  = 3e-2
          ilu%droptolS = ilu%droptol / 10
          ilu%matching = 0
          ilu%restol   = 1e-14

          call sA%factorize()
          call sA%solve_vec(b, x)
          call tc%assert_eq(x_exp, x, 1e-15, "solve_vec cmplx: ilupack solver")
        end block
      m4_divert(0)
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

      b     = reshape([1,   -4,   7,  -10,   -2,   5,    -8,  11,    3,   -6,   9,  -12],   [4, 3])
      x_exp = reshape([3.0, -1.0, 7.0, -1.8, -4.5, 1.25, -8.0, 1.95, 6.0, -1.5, 9.0, -2.1], [4, 3])
      allocate (x(4,3))

      call example_matrix3(sA)
      call sA%factorize()
      call sA%solve_mat(b, x)

      call tc%assert_eq(x_exp, x, 1e-15, "solve_mat")
    end block

    ! add_sparse: A <- A + fact * B
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
      call tc%assert_eq(a_exp,  sA%a,  0.0, "add_sparse empty 1: a")
      call tc%assert_eq(ia_exp, int(sA%ia),      "add_sparse empty 1: ia")
      call tc%assert_eq(ja_exp, sA%ja,      "add_sparse empty 1: ja")

      call get_test_matrix2(sA)
      a_exp = 3 * sA%a
      call matrix_add(sA, sE, fact=3.0)
      call tc%assert_eq(a_exp,  sE%a,  0.0, "add_sparse empty 2: a")
      call tc%assert_eq(ia_exp, int(sE%ia), "add_sparse empty 2: ia")
      call tc%assert_eq(ja_exp, sE%ja,      "add_sparse empty 2: ja")

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

      call tc%assert_eq(a_exp,  sA%a, 0.0, "add_sparse: a")
      call tc%assert_eq(ia_exp, int(sA%ia),     "add_sparse: ia")
      call tc%assert_eq(ja_exp, sA%ja,     "add_sparse: ja")
    end block

    ! add_sparse3: C <- fact1 * A + fact2 * B
    block
      integer, allocatable :: ia_exp(:), ja_exp(:)
      real,    allocatable :: a_exp(:)
      type(sparse_real)    :: sA, sB, sC, sE

      call example_matrix3(sA)
      call get_test_matrix2(sB)
      call get_empty_matrix(sE)

      call matrix_add(sA, sE, sC, fact1 = 3.0, fact2 = -5.0)
      call tc%assert_eq(3.0*sA%a  ,     sC%a  , 1e-15, "add_sparse3 empty 1: a")
      call tc%assert_eq(3.0*sA%a  ,     sC%a  , 1e-15, "add_sparse3 empty 1: a")
      call tc%assert_eq(int(sA%ia), int(sC%ia),        "add_sparse3 empty 1: ia")
      call tc%assert_eq(    sA%ja ,     sC%ja ,        "add_sparse3 empty 1: ja")

      call matrix_add(sE, sA, sC, fact1=-5.0, fact2=3.0)
      call tc%assert_eq(3.0*sA%a  ,     sC%a  , 1e-15, "add_sparse3 empty 2: a")
      call tc%assert_eq(3.0*sA%a  ,     sC%a  , 1e-15, "add_sparse3 empty 2: a")
      call tc%assert_eq(int(sA%ia), int(sC%ia),        "add_sparse3 empty 2: ia")
      call tc%assert_eq(    sA%ja ,     sC%ja ,        "add_sparse3 empty 2: ja")

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

      call tc%assert_eq(a_exp,  sC%a, 0.0,  "add_sparse3: a")
      call tc%assert_eq(ia_exp, int(sC%ia), "add_sparse3: ia")
      call tc%assert_eq(ja_exp, sC%ja,      "add_sparse3: ja")
    end block

    ! add_band: S <- S + fact * B
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
      call tc%assert_eq(a_exp,  S%a, 1e-16, "add_band empty: a")
      call tc%assert_eq(ia_exp, int(S%ia),  "add_band empty: ia")
      call tc%assert_eq(ja_exp, S%ja,       "add_band empty: ja")

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
      call tc%assert_eq(a_exp,  S%a, 0.0,  "add_band: a")
      call tc%assert_eq(ia_exp, int(S%ia), "add_band: ia")
      call tc%assert_eq(ja_exp, S%ja,      "add_band: ja")
    end block

    ! add_band3: S2 <- fact1 * S + fact2 * B
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
      call tc%assert_eq(a_exp,  S2%a, 1e-16, "add_band3 empty: a")
      call tc%assert_eq(ia_exp, int(S2%ia),  "add_band3 empty: ia")
      call tc%assert_eq(ja_exp, S2%ja,       "add_band3 empty: ja")

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

      call tc%assert_eq(a_exp,  S2%a, 1e-16, "add_band3: a")
      call tc%assert_eq(ia_exp, int(S2%ia),  "add_band3: ia")
      call tc%assert_eq(ja_exp, S2%ja,       "add_band3: ja")
    end block

    ! zero
    block
      integer, parameter   :: n = 3
      real                 :: z_exp(n,n)
      type(dense_real)     :: z_d
      type(sparse_real)    :: z_sp

      z_sp = sparse_zero_real(n)

      z_exp = 0

      call z_d%init(n)
      call matrix_convert(z_sp, z_d)

      call tc%assert_eq(z_exp, z_d%d, 0.0, "zero(3)")
    end block

    ! nnz
    block
      type(sparse_real)  :: s
      type(spbuild_real) :: spb

      call get_empty_matrix(s)
      call tc%assert_eq(0, int(s%nnz()), "nnz 1")

      call example_matrix3(s)
      call tc%assert_eq(6, int(s%nnz()), "nnz 2")

      call example_matrix3(s)
      call tc%assert_eq(6, int(s%nnz(only_nonzeros = .true.)), "nnz 3")

      call get_test_matrix2(s)
      call tc%assert_eq(5, int(s%nnz()), "nnz 4")

      call get_test_matrix2(s)
      call tc%assert_eq(5, int(s%nnz(only_nonzeros = .true.)), "nnz 5")

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
      call spb%save()

      call tc%assert_eq(6, int(s%nnz()), "nnz 6")
    end block

    ! eye
    block
      integer            :: i, ncols
      integer, parameter :: nrows = 10
      real, allocatable  :: x(:), y(:), y_exp(:)
      type(sparse_real)  :: S

      ! test 1: only supply nrows
      S = sparse_eye_real(nrows)

      allocate (x(nrows), y(nrows), y_exp(nrows))
      do i = 1, nrows
        x     = 0
        x(i)  = 1
        y_exp = x
        call S%mul_vec(x, y)

        call tc%assert_eq(y_exp, y, 0.0, "eye: only nrows")
      end do
      deallocate (x, y, y_exp)

      ! test 2: nrows > ncols
      ncols = nrows-2
      S = sparse_eye_real(nrows, ncols=ncols)

      allocate (x(ncols), y(nrows), y_exp(nrows))
      do i = 1, ncols
        x    = 0
        x(i) = 1
        y_exp(:ncols)   = x
        y_exp(ncols+1:) = 0
        call S%mul_vec(x, y)

        call tc%assert_eq(y_exp, y, 0.0, "eye: nrows>ncols")
      end do
      deallocate (x, y, y_exp)

      ! test 2: nrows < ncols
      ncols = nrows+2
      S = sparse_eye_real(nrows, ncols=ncols)

      allocate (x(ncols), y(nrows), y_exp(nrows))
      do i = 1, ncols
        x     = 0
        x(i)  = 1
        y_exp = x(:nrows)
        call S%mul_vec(x, y)

        call tc%assert_eq(y_exp, y, 0.0, "eye: nrows<ncols")
      end do
      deallocate (x, y, y_exp)
    end block

    call tc%finish()
  end subroutine

  subroutine get_empty_matrix(S)
    type(sparse_real), intent(out) :: S

    call S%init(4)
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

    call sbuild%save()
  end subroutine

  subroutine get_test_matrix3(n, S)
    !! some large sparse matrix

    integer,           intent(in)          :: n
      !! matrix dimension
    type(sparse_real), intent(out) :: S

    integer            :: i, j
    type(spbuild_real) :: sbuild

    ! S_ij = /i+2*j-i**2    if    (i+j)%31 == 0
    !        \0             else

    call S%init(n)
    call sbuild%init(S)

    do i = 1, n
      do j = 1, n
        if (mod(i+j, 31) == 0) call sbuild%add(i, j, real(i+2*j-i**2))
      end do
    end do

    call sbuild%save()
  end subroutine

end submodule
