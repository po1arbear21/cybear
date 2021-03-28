#ifdef USE_MUMPS
module test_mumps_m
  use mumps_m
  use sparse_idx_m
  use test_case_m
  implicit none

  private
  public test_mumps

contains

  subroutine test_mumps()
    type(test_case) :: tc

    call tc%init("mumps")

    ! test real mumps
    block
      integer                          :: m, n, n_rhs, i_rhs
      integer,             allocatable :: ja(:)
      integer(SPARSE_IDX), allocatable :: ia(:)
      real,                allocatable :: a(:), b(:,:), x_exp(:,:), x(:)

      !
      ! build matrix
      !     / 2 0 0 1 \
      ! M = | 5 8 0 0 |
      !     | 0 0 3 0 |
      !     \ 0 8 0 5 /
      !
      n  = 4
      a  = [ 2.0, 1.0, 5.0, 8.0, 3.0, 8.0, 5.0 ]
      ja = [   1,   4,   1,   2,   3,   2,   4 ]
      ia = int([   1,   3,   5,   6,   8 ], kind = SPARSE_IDX)

      !
      ! rhs
      !
      n_rhs = 5
      allocate (b(n,n_rhs))
      b(:,1) = [ 1.0, 2.0, 3.0, 4.0 ]
      b(:,2) = [ 1.0, 3.0, 3.0, 4.0 ]
      b(:,3) = [ 9.0, 2.0, 3.0, 4.0 ]
      b(:,4) = [ 3.0, 1.0, 5.0, 6.0 ]
      b(:,5) = [ 1.0, 3.0, 7.0, 2.0 ]

      !
      ! sols
      !
      allocate (x_exp(n,n_rhs), x(n))
      x_exp(:,1) = [ 0.200000000000000,  0.125000000000000, 1.000000000000000, 0.600000000000000 ]
      x_exp(:,2) = [ 0.266666666666667,  0.208333333333333, 1.000000000000000, 0.466666666666667 ]
      x_exp(:,3) = [ 2.866666666666666, -1.541666666666667, 1.000000000000000, 3.266666666666667 ]
      x_exp(:,4) = [ 0.666666666666667, -0.291666666666667, 1.666666666666667, 1.666666666666667 ]
      x_exp(:,5) = [ 0.400000000000000,  0.125000000000000, 2.333333333333333, 0.200000000000000 ]

      !
      ! mumps: init -> factor -> solve -> delete
      !
      m = create_mumps_handle_r()
      call mumps_factorize(m, ia, ja, a)

      do i_rhs = 1, n_rhs
        call mumps_solve(m, b(:,i_rhs), x)
        call tc%assert_eq(x_exp(:,i_rhs), x, 1e-13, "real solving")
      end do

      call destruct_mumps_handle_r(m)
    end block

    ! test complex mumps
    block
      integer                          :: m, n, n_rhs, i_rhs
      integer,             allocatable :: ja(:)
      integer(SPARSE_IDX), allocatable :: ia(:)
      complex,             allocatable :: a(:), b(:,:), x_exp(:,:), x(:)

      !
      ! build matrix
      !     / 2i 0 0  1+i \
      ! M = | 5  8 0  0   |
      !     | 0  0 3i 0   |
      !     \ 0  8 0  5   /
      !
      n  = 4
      a  = [ (0.0,2.0), (1.0,1.0), (5.0,0.0), (8.0,0.0), (0.0,3.0), (8.0,0.0), (5.0,0.0) ]
      ja = [   1,   4,   1,   2,   3,   2,   4 ]
      ia = int([   1,   3,   5,   6,   8 ], kind = SPARSE_IDX)

      !
      ! rhs
      !
      n_rhs = 3
      allocate (b(n,n_rhs))
      b(:,1) = 0
      b(:,2) = [ (1.0,1.0), (2.0,1.0), (0.0,1.0), (1.0,-2.0) ]
      b(:,3) = [ (1.0,0.0), (2.0,1.0), (3.0,1.0), (4.0, 0.0) ]

      !
      ! sols
      !
      allocate (x_exp(n,n_rhs), x(n))
      x_exp(:,1) = 0
      x_exp(:,2) = [ ( 0.600000000000000, 0.000000000000000), (-0.125000000000000,0.125000000000000), (0.333333333333333, 0.000000000000000), (0.400000000000000,-0.600000000000000) ]
      x_exp(:,3) = [ (-0.020000000000000,-0.140000000000000), ( 0.262500000000000,0.212500000000000), (0.333333333333333,-1.000000000000000), (0.380000000000000,-0.340000000000000) ]

      !
      ! mumps: init -> factor -> solve -> delete
      !
      m = create_mumps_handle_c()
      call mumps_factorize(m, ia, ja, a)

      do i_rhs = 1, n_rhs
        call mumps_solve(m, b(:,i_rhs), x)
        call tc%assert_eq(x_exp(:,i_rhs), x, 1e-13, "complex solving")
      end do

      call destruct_mumps_handle_c(m)
    end block

    call tc%finish
  end subroutine

end module
#endif
