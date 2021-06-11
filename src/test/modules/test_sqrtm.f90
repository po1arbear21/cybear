module test_sqrtm_m

  use matrix_m,    only: dense_real, hessenberg_real
  use sqrtm_m,     only: sqrtm
  use test_case_m, only: test_case

  implicit none

  private
  public test_sqrtm

contains

  subroutine test_sqrtm
    type(test_case)       :: tc
    real, allocatable     :: d0(:,:), Asqrt_exp(:,:)
    type(dense_real)      :: Asqrt

    call tc%init("sqrtm")

    !
    ! test 1: sqrtm for dense matrix A
    !
    block
      type(dense_real) :: A

      ! magic(3)^2 =
      ! 91    67    67
      ! 67    91    67
      ! 67    67    91
      d0 = reshape([91,67,67,67,91,67,67,67,91], [3,3])
      call A%init(d0)

      ! sqrtm(magic(3)*magic(3)) =
      !  8.265986323710905   3.367006838144549   3.367006838144548
      !  3.367006838144549   8.265986323710909   3.367006838144549
      !  3.367006838144548   3.367006838144549   8.265986323710905
      Asqrt_exp = reshape([ 8.265986323710905, 3.367006838144549, 3.367006838144548, &
        &                   3.367006838144549, 8.265986323710909, 3.367006838144549, &
        &                   3.367006838144548, 3.367006838144549, 8.265986323710905], [3,3])

      call sqrtm(A, Asqrt)

      call tc%assert_eq(Asqrt%d, Asqrt_exp, 1e-14, "sqrtm(dense)")
    end block

    !
    ! test 2: sqrtm for upper Hessenberg matrix A
    !
    block
      type(hessenberg_real) :: A

      ! A =
      ! 91    67    67
      ! 67    91    67
      ! 0    67    91
      d0 = reshape([91,67,0,67,91,67,67,67,91], [3,3])
      call A%init(3, .true.)
      A%d = d0

      ! sqrtm(A) =
      ! 8.998743273758659   3.165852095879277   3.165852095879274
      ! 4.099763788192307   8.064831581445633   3.165852095879278
      ! -0.933911692313032   4.099763788192310   8.998743273758668
      Asqrt_exp = reshape([ 8.998743273758659, 4.099763788192307, -0.933911692313032, &
        &                   3.165852095879277, 8.064831581445633,  4.099763788192310, &
        &                   3.165852095879274, 3.165852095879278,  8.998743273758668], [3,3])

      call sqrtm(A, Asqrt)

      call tc%assert_eq(Asqrt%d, Asqrt_exp, 1e-14, "sqrtm(hessenberg)")
    end block

    call tc%finish
  end subroutine

end module
