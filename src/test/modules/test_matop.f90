module test_matop_m

  use matop_m
  use matrix_m
  use test_case_m

  implicit none

  private
  public test_matop

contains

  subroutine test_matop()
    type(test_case) :: tc

    print "(A)", "test_matop"
    call tc%init("matop")

    ! test: single_matop_real
    block
      real                    :: x1(2), y1(2), y1_exp(2), x2(2,2), y2(2,2), y2_exp(2,2)
      type(dense_real)        :: d
      type(single_matop_real) :: sop

      ! dense matrix:
      ! d = [1 0]
      !     [1 1]
      call d%init(reshape(real([1, 1, 0, 1]), [2, 2]))
      call sop%init(d)

      !
      ! test 1: check matrix vector multiplications
      !
      ! d * x1 = y1 ==? y1_exp

      x1     = [1, 0]
      y1_exp = [1, 1]
      call sop%exec(x1, y1)
      call tc%assert_eq(y1_exp, y1, 1e-13, "single_matop_real: exec1 == mul_vec")

      x1     = [0, 1]
      y1_exp = [0, 1]
      call sop%exec(x1, y1)
      call tc%assert_eq(y1_exp, y1, 1e-13, "single_matop_real: exec1 == mul_vec")

      !
      ! test 2: check matrix matrix multiplications
      !
      ! d * x = [1 0] [1 2] = [1 2]
      !         [1 1] [3 4]   [4 6]
      x2     = reshape([1, 3, 2, 4], [2, 2])
      y2_exp = reshape([1, 4, 2, 6], [2, 2])
      call sop%exec(x2, y2)
      call tc%assert_eq(y2_exp, y2, 1e-13, "single_matop_real: exec2 == mul_mat")

      !
      ! test 3: check solving for a vector
      !
      ! d * x = b = [1]   note: d^-1 = [ 1 0]
      !             [2]                [-1 1]
      call d%factorize
      call sop%init(d, inv=.true.)

      ! solve: d y1 = x1
      x1     = [1, 2]
      y1_exp = [1, 1]
      call sop%exec(x1, y1)
      call tc%assert_eq(y1_exp, y1, 1e-13, "single_matop_real: exec1 == solve_vec")

      !
      ! test 4: check solving for a matrix
      !
      ! d * x = b = [1 3]   note: d^-1 = [ 1 0]
      !             [2 4]                [-1 1]

      ! solve: d y2 = x2
      x2     = reshape([1, 2, 3, 4], [2, 2])
      y2_exp = reshape([1, 1, 3, 1], [2, 2])
      call sop%exec(x2, y2)
      call tc%assert_eq(y2_exp, y2, 1e-13, "single_matop_real: exec2 == solve_mat")
    end block

    ! test: single_matop_cmplx
    block
      complex                  :: x1(2), y1(2), y1_exp(2), x2(2,2), y2(2,2), y2_exp(2,2)
      type(dense_cmplx)        :: d
      type(single_matop_cmplx) :: sop

      ! dense matrix:
      ! d = [    i  0  ]
      !     [0.5+i   -i]
      call d%init(reshape([(0,1), (0.5,1), (0,0), (0,-1)], [2, 2]))
      call sop%init(d)

      !
      ! test 1: check matrix vector multiplications
      !
      ! d * x1 = y1 ==? y1_exp

      x1     = [( 1,1), ( 0,  0  )]
      y1_exp = [(-1,1), (-0.5,1.5)]
      call sop%exec(x1, y1)
      call tc%assert_eq(y1_exp, y1, 1e-13, "single_matop_cmplx: exec1 == mul_vec")

      x1     = [(0,0), (1, 1)]
      y1_exp = [(0,0), (1,-1)]
      call sop%exec(x1, y1)
      call tc%assert_eq(y1_exp, y1, 1e-13, "single_matop_cmplx: exec1 == mul_vec")

      !
      ! test 2: check matrix matrix multiplications
      !
      ! d * x = [    i  0  ] [1    2+i] = [i       -1+2i  ]
      !         [0.5+i   -i] [3-i  4  ]   [-0.5-2i   -1.5i]
      x2     = reshape([(1,0), ( 3,  -1), ( 2,1), (4, 0  )], [2, 2])
      y2_exp = reshape([(0,1), (-0.5,-2), (-1,2), (0,-1.5)], [2, 2])
      call sop%exec(x2, y2)
      call tc%assert_eq(y2_exp, y2, 1e-13, "single_matop_cmplx: exec2 == mul_mat")

      !
      ! test 3: check solving for a vector
      !
      ! d * x = b = [1  ]   note: d^-1 = [ -i     0]
      !             [2+i]                [-0.5-i  i]
      call d%factorize
      call sop%init(d, inv=.true.)

      ! solve: d y1 = x1
      x1     = [(1, 0), ( 2,  1)]
      y1_exp = [(0,-1), (-1.5,1)]
      call sop%exec(x1, y1)
      call tc%assert_eq(y1_exp, y1, 1e-13, "single_matop_cmplx: exec1 == solve_vec")

      !
      ! test 4: check solving for a matrix
      !
      ! d * x = b = [1   3-i ]   note: d^-1 = [ -i     0]
      !             [2+i 4+5i]                [-0.5-i  i]

      ! solve: d y2 = x2
      x2     = reshape([(1, 0), ( 2,  1), ( 3,-1), ( 4,  5  )], [2, 2])
      y2_exp = reshape([(0,-1), (-1.5,1), (-1,-3), (-7.5,1.5)], [2, 2])
      call sop%exec(x2, y2)
      call tc%assert_eq(y2_exp, y2, 1e-13, "single_matop_cmplx: exec2 == solve_mat")
    end block

    call test_matop_c2r(tc)

    call tc%finish()
  end subroutine

  subroutine test_matop_c2r(tc)
    !! tests type matop_c2r

    type(test_case), intent(inout) :: tc

    type(dense_cmplx)        :: d
    type(single_matop_cmplx) :: sop
    type(matop_c2r)          :: c2r

    ! dense matrix:
    ! d = [    i  0  ]
    !     [0.5+i   -i]
    call d%init(reshape([(0,1), (0.5,1), (0,0), (0,-1)], [2, 2]))
    call sop%init(d)
    call c2r%init(sop)

    ! test exec1
    block
      ! y=d*x for x=[1+2i] => x2r = [1, 2, 3, 4]^T
      !             [3+4i]
      !  => y = [-2  +i] => y2r = [-2, 1, 2.5, -1]^T
      !         [ 2.5-i]

      real :: x2r(4), y2r(4), y2r_exp(4)

      x2r = [1, 2, 3, 4]
      call c2r%exec(x2r, y2r)
      y2r_exp = [-2.0, 1.0, 2.5, -1.0]

      call tc%assert_eq(y2r_exp, y2r, 1e-13, "matop_c2r: exec1")
    end block

    ! test exec2
    block
      ! y=d*x for x=[1+2i, 5+6i] => x2r = [1, 2, 3, 4]^T
      !             [3+4i, 7+8i]          [5, 6, 7, 8]
      !  => y = [-2  +i, -6  +5i] => y2r = [-2, 1, 2.5, -1]^T
      !         [ 2.5-i,  4.5+ i]          [-6, 5, 4.5,  1]

      real :: x2r(4,2), y2r(4,2), y2r_exp(4,2)

      x2r(:,1) = [1, 2, 3, 4]
      x2r(:,2) = [5, 6, 7, 8]
      call c2r%exec(x2r, y2r)
      y2r_exp(:,1) = [-2.0, 1.0, 2.5, -1.0]
      y2r_exp(:,2) = [-6.0, 5.0, 4.5,  1.0]

      call tc%assert_eq(y2r_exp, y2r, 1e-13, "matop_c2r: exec2")
    end block
  end subroutine

end module
