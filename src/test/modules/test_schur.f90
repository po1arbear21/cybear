module test_schur_m
  use test_case_m
  use schur_m

  implicit none

contains
  subroutine test_schur
    type(test_case) :: tc
    real, allocatable, dimension(:,:) :: d0, Q_exp, U_exp
    type(dense_real) :: A, U, Q

    print "(1A)", "test_schur"
    call tc%init("schur")

    d0 = reshape([3,2,1,4,2,1,4,4,0], [3, 3], order=[2, 1])
    call A%init(d0)

    Q_exp = reshape([-0.498573734016578, -0.764694425952027,  0.408248290463864, &
      &              -0.574051724153223, -0.061627520881869, -0.816496580927726, &
      &              -0.649529714289869,  0.641439384188292,  0.408248290463862  ], [3, 3], order=[2, 1])

    U_exp = reshape([6.605551275463988,  4.490731195102494,  0.826321946833888, &
      &              0.0              , -0.605551275463992, -1.072625458169801, &
      &              0.0              ,  0.0              , -0.999999999999999  ], [3, 3], order=[2, 1])

    call schur(A, Q, U)

    call tc%assert_eq(Q%d, Q_exp, 1e-10, "mat Q")
    call tc%assert_eq(U%d, U_exp, 1e-10, "mat U")

    call tc%finish
  end subroutine
end module
