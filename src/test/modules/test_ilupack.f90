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
    a  = [ 2,1,5,8,3,8,5 ]
    ja = [ 1,4,1,2,3,2,4 ]
    ia = [ 1,3,5,6,8 ]

    !
    ! rhs
    !
    n_rhs = 5
    allocate (b(n,n_rhs))
    b(:,1) = [1,2,3,4]
    b(:,2) = [1,3,3,4]
    b(:,3) = [9,2,3,4]
    b(:,4) = [10,2,3,4]
    b(:,5) = [11,2,3,4]

    !
    ! sols
    !
    allocate (x_exp(n,n_rhs), x(n))
    x_exp(:,1) = [0.2,                0.125,             1.0,               0.6]
    x_exp(:,2) = [0.266666666666667,  0.208333333333333, 1.000000000000000, 0.466666666666667]
    x_exp(:,3) = [2.866666666666666, -1.541666666666667, 1.000000000000000, 3.266666666666666]
    x_exp(:,4) = [3.200000000000000, -1.750000000000000, 1.000000000000000, 3.600000000000000]
    x_exp(:,5) = [3.533333333333333, -1.958333333333333, 1.000000000000000, 3.933333333333333]

    !
    ! ilu: init -> factor -> solve -> delete
    !
    call ilu%init(  a, ia, ja)
    call ilu%factor(a, ia, ja)

    do i_rhs = 1, n_rhs
      call ilu%solve(a, ia, ja, b(:,i_rhs), x)

      call tc%assert_eq(x_exp(:,i_rhs), x, 1e-13, "solving")
    end do

    call ilu%delete

    call tc%finish
  end subroutine

end module
