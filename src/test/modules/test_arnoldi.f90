module test_arnoldi_m
  use arnoldi_m
  use matop_m
  use matrix_m
  use test_case_m
  implicit none

  type(test_case) :: tc

contains

  subroutine test_arnoldi
    ! type(arnoldi_int) :: vec
    real, allocatable :: d0(:,:), b(:), Q_exp(:,:), H_exp(:,:)
    integer :: brkd

    type(band_real) :: A
    type(dense_real) :: H, Q
    type(single_matop_real) :: matop_A

    print "(A)", "test_arnoldi"
    call tc%init("arnoldi")

    ! A = diag(..)
    allocate(d0(1,3))
    d0(1,:) = [1, 2, 3]
    call A%init(3,0,d0=d0)
    call matop_A%init(A)

    ! b vec
    b = [1,1,1]

    ! H, Q
    call Q%init(3, ncols=3)
    call H%init(3, ncols=2)

    call arnoldi(matop_A, b, H, Q, brkd)

    ! expected
    Q_exp = reshape([ 0.577350269189626,  -0.707106781186548,   0.408248290463863, &
      &               0.577350269189626,  -0.000000000000000,  -0.816496580927726, &
      &               0.577350269189626,   0.707106781186547,   0.408248290463864   ], [3, 3], order=[2, 1] )

    H_exp = reshape([ 2.000000000000000,   0.816496580927726, &
      &               0.816496580927726,   2.000000000000000, &
      &               0.0              ,   0.577350269189626  ], [3, 2], order=[2, 1])


    call tc%assert_eq(H%d, H_exp, 1e-13, "mat H")
    call tc%assert_eq(Q%d, Q_exp, 1e-13, "mat Q")

    call tc%finish
  end subroutine

end module
