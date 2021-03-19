#ifdef USE_ILUPACK

module test_ilupack_m
  use ilupack_m
  use sparse_idx_m, only: SPARSE_IDX
  use test_case_m

  implicit none

  private
  public test_ilupack

contains

  subroutine test_ilupack()
    type(test_case) :: tc

    print "(A)", "test_ilupack"
    call tc%init("ilupack")

    ! test ilupack_real
    block
      integer                           :: ilu, n, n_rhs, i_rhs
      integer(SPARSE_IDX),  allocatable :: ia(:)
      integer,              allocatable :: ja(:)
      real,                 allocatable :: a(:), b(:,:), x_exp(:,:), x(:)
      type(ilupack_handle), pointer     :: ilu_h

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
      ia = [   1,   3,   5,   6,   8 ]

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
      ! ilu: init -> factor -> solve -> delete
      !
      ilu = create_ilupack_handle(ia, ja, a)
      call get_ilupack_handle_ptr(ilu, ilu_h)
      call ilupack_factorize(ilu)

      call tc%assert_eq(8, ilu_h%nnz(), "real nnz")

      do i_rhs = 1, n_rhs
        call ilupack_solve(ilu, b(:,i_rhs), x)
        call tc%assert_eq(x_exp(:,i_rhs), x, 1e-13, "real solving")
      end do

      call destruct_ilupack_handle(ilu)
    end block

    ! test ilupack_cmplx
    block
      integer                           :: ilu, n, n_rhs, i_rhs
      integer(SPARSE_IDX),  allocatable :: ia(:)
      integer,              allocatable :: ja(:)
      complex,              allocatable :: a(:), b(:,:), x_exp(:,:), x(:)

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
      ia = [   1,   3,   5,   6,   8 ]

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
      ! ilu: init -> factor -> solve -> delete
      !
      ilu = create_ilupack_handle(ia, ja, a)
      call ilupack_factorize(ilu)

      do i_rhs = 1, n_rhs
        call ilupack_solve(ilu, b(:,i_rhs), x)
        call tc%assert_eq(x_exp(:,i_rhs), x, 1e-13, "complex solving")
      end do

      call destruct_ilupack_handle(ilu)
    end block

    call tc%finish()
  end subroutine

end module

#endif
