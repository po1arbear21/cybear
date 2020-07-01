module test_ilupack_m

  use test_case_m
  use ilupack_m

  implicit none

contains

  subroutine test_ilupack
    type(test_case)      :: tc
    type(ilupack)        :: ilu
    integer              :: n, n_rhs, i_rhs
    integer, allocatable :: ia(:), ja(:)
    real,    allocatable :: a(:), b(:,:), x_exp(:,:), x(:)

    print "(1A)", "test_ilupack"
    call tc%init("ilupack")

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
    call ilu%init(ia, ja, a)
    call ilu%factor

    do i_rhs = 1, n_rhs
      x = x_exp(1,i_rhs)
      call ilu%solve(b(:,i_rhs), x)
      call tc%assert_eq(x_exp(:,i_rhs), x, 1e-13, "solving")
    end do

    call ilu%delete

    call tc%finish
  end subroutine

end module
