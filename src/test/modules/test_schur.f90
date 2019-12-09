module test_schur_m
  use test_case_m
  use schur_m

  implicit none

contains
  subroutine test_schur
    type(test_case)       :: tc
    real, allocatable     :: d0(:,:), Q_exp(:,:), U_exp(:,:)
    type(dense_real)      :: Q
    type(hessenberg_real) :: U

    print "(1A)", "test_schur"
    call tc%init("schur")

    !
    ! test 1: schur for dense matrix A
    !
    block
      type(dense_real) :: A

      d0 = reshape([3,2,1,4,2,1,4,4,0], [3, 3], order=[2, 1])
      call A%init(d0)

      Q_exp = reshape([-0.498573734016578, -0.764694425952027,  0.408248290463864, &
        &              -0.574051724153223, -0.061627520881869, -0.816496580927726, &
        &              -0.649529714289869,  0.641439384188292,  0.408248290463862  ], [3, 3], order=[2, 1])

      U_exp = reshape([6.605551275463988,  4.490731195102494,  0.826321946833888, &
        &              0.0              , -0.605551275463992, -1.072625458169801, &
        &              0.0              ,  0.0              , -0.999999999999999  ], [3, 3], order=[2, 1])

      call schur(A, Q, U)

      call tc%assert_eq(Q%d, Q_exp, 1e-12, "schur(dense): mat Q")
      call tc%assert_eq(U%d, U_exp, 1e-12, "schur(dense): mat U")
    end block

    !
    ! test 2: schur for upper Hessenberg matrix A
    !
    block
      type(hessenberg_real) :: A

      ! A =
      ! 1     3     6    10
      ! 2     4     7    11
      ! 0     5     8    12
      ! 0     0     9    13
      d0 = reshape([1,2,0,0,3,4,5,0,6,7,8,9,10,11,12,13], [4,4])
      call A%init(4, .true.)
      A%d = d0

      ! [Q,U] = schur(A)
      ! Q =
      !     0.4481   -0.3332    0.8294    0.0131
      !    -0.1640    0.8724    0.4367    0.1462
      !     0.7375    0.2077   -0.3238    0.5551
      !    -0.4779   -0.2912    0.1283    0.8187
      ! U =
      !    -0.8882    2.8818    1.6758    6.7350
      !          0    1.8142    2.3300    6.3520
      !          0   -0.0830    1.8142   12.9752
      !          0         0         0   23.2597

      Q_exp = reshape([ 0.448145049702056, -0.163994585129952,  0.737478720162248, -0.477908911596602, &
        &              -0.333183220301448,  0.872367440909075,  0.207735739631395, -0.291221311430168, &
        &               0.829447116463941,  0.436715428401495, -0.323797104959479,  0.128267495514697, &
        &               0.013131531774723,  0.146152176466542,  0.555091375007314,  0.818743347807988], [4,4])

      U_exp = reshape([-0.888229159163929, 0.0,                0.0,                0.0, &
        &               2.881836001894352, 1.814249536124565, -0.083046415070944,  0.0, &
        &               1.675754382182381, 2.330007343752856,  1.814249536124565,  0.0, &
        &               6.735020782689861, 6.351984917558429, 12.975200296321226, 23.259730086914800], [4,4])

      call schur(A, Q, U)

      call tc%assert_eq(Q%d, Q_exp, 1e-12, "schur(Hessenberg): mat Q")
      call tc%assert_eq(U%d, U_exp, 1e-12, "schur(Hessenberg): mat U")
    end block

    call tc%finish
  end subroutine
end module
