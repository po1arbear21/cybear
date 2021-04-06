module test_mkl_ilu_m

  use example_matrices_m, only: matrix2
  use gmres_m,            only: gmres
  use matop_m,            only: single_matop_real
  use matrix_m,           only: sparse_real
  use mkl_ilu_m,          only: mkl_ilu0, mkl_ilut
  use test_case_m,        only: test_case
  use util_m,             only: int2str

  implicit none

  private
  public test_mkl_ilu

contains

  subroutine test_mkl_ilu()
    type(test_case) :: tc

    call tc%init("mkl_ilu")

    call test_ilu0(tc)
    call test_ilut(tc)

    call tc%finish()
  end subroutine

  subroutine test_ilu0(tc)
    !! tests the ilu0 preconditioner
    !!
    !! scaled lower diagonal of exact matrix.
    !!  scaling results into smaller element of lu-decomposition outside of sparsity strucutre of exact A.
    !!  thus, hopefully better convergence of ilu.
    !! still solution is very bad: even though almost exact solution is init guess we only get the expected solution
    !!  upto absolute threshold of 1e-6.
    type(test_case), intent(inout) :: tc

    real, parameter :: b(*)     = [1, 2, 3, 4, 5], &
      &                x_exp(*) = [3.034413115792796e-03, 1.505796570568668e-01, -2.914538122392946e-02, &
      &                            1.014944360842705e-01, 2.797668253999683e-01]

    integer                 :: i
    real                    :: x(5)
    real, allocatable       :: res(:)
    type(mkl_ilu0)          :: ilu
    type(single_matop_real) :: A
    type(sparse_real)       :: Asp

    ! init exact matop
    call matrix2(Asp, fact_lodiag=1.0/2.4)
    call A%init(Asp)

    ! test 1, 2
    do i = 1, 2
      ! compute ilu
      if (i == 1) then
        ! test 1: save pointer to sparse matrix
        call ilu%init(Asp)

      else
        ! test 2: copy exact matrix' values
        block
          type(sparse_real) :: Asp_tmp

          call matrix2(Asp_tmp, fact_lodiag=1.0/2.4)
          call ilu%init(Asp_tmp, copy=.true.)
        end block
      end if

      ! initial guess is exact solution plus small derivation.
      x     = x_exp
      x(1)  = x(1) * (1 + 1e-2)

      ! solve
      call gmres(b, A, x, precon=ilu, residual=res)
      call tc%assert_eq(x_exp, x, 1e-6, "ilu0: test i:"//int2str(i))
    end do
  end subroutine

  subroutine test_ilut(tc)
    !! tests the ilut preconditioner
    !!
    !! scaled lower diagonal of exact matrix.
    !!  scaling results into smaller element of lu-decomposition outside of sparsity strucutre of exact A.
    !!  thus, hopefully better convergence of ilu.
    !! still solution is very bad: even though almost exact solution is init guess we only get the expected solution
    !!  upto absolute threshold of 1e-6.
    type(test_case), intent(inout) :: tc

    integer, parameter :: maxfil = 1
    real,    parameter :: b(*)     = [1, 2, 3, 4, 5],                                                       &
      &                   x_exp(*) = [3.034413115792796e-03, 1.505796570568668e-01, -2.914538122392946e-02, &
      &                               1.014944360842705e-01, 2.797668253999683e-01],                        &
      &                   tol      = 1e-11

    integer                 :: i
    real                    :: x(5)
    real, allocatable       :: res(:)
    type(mkl_ilut)          :: ilu
    type(single_matop_real) :: A
    type(sparse_real)       :: Asp

    ! init exact matop
    call matrix2(Asp, fact_lodiag=1.0/2.4)
    call A%init(Asp)

    ! test 1, 2
    do i = 1, 2
      ! compute ilu
      if (i == 1) then
        ! test 1: save pointer to sparse matrix
        call ilu%init(Asp, tol, maxfil)

      else
        ! test 2: copy exact matrix' values
        block
          type(sparse_real) :: Asp_tmp

          call matrix2(Asp_tmp, fact_lodiag=1.0/2.4)
          call ilu%init(Asp_tmp, tol, maxfil, copy=.true.)
        end block
      end if

      ! initial guess is exact solution plus small derivation.
      x     = x_exp
      x(1)  = x(1) * (1 + 1e-2)

      ! solve
      call gmres(b, A, x, precon=ilu, residual=res)
      call tc%assert_eq(x_exp, x, 1e-6, "ilut: test i:"//int2str(i))
    end do
  end subroutine
end module
