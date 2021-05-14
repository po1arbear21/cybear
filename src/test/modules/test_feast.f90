#ifdef USE_FEAST

module test_feast_m
  use feast_m,     only: feast, feast_option
  use test_case_m, only: test_case

  implicit none

  private
  public test_feast

contains

  subroutine test_feast()
    type(feast_option) :: opt
    type(test_case)    :: tc

    integer, parameter   :: N = 10
    complex              :: A(N,N), Azinv(N,N)
    complex, allocatable :: eval(:), evecR(:,:)

    call tc%init("feast")

    ! build A, Azinv
    block
      integer :: i

      A      = 0
      A(1,1) = (1, -1)
      A(2,2) = (1,  1)
      A(3,3) = -2
      A(4,4) = -10
      do i = 5, 10
        A(i,i) = 100+i
      end do

      Azinv  = 0
    end block

    call opt%init(N, 4)
    opt%mid = (1, 2)
    opt%rad = 1.5

    call feast(opt, fact_Az, solve_Az, mulvec_A, eval, evecR=evecR)

    block
      complex :: evecR_exp(N,1) = 0

      call tc%assert_eq([(1, 1)], eval, 1e-15, "eval")
      evecR_exp(2,1) = 1
      call tc%assert_eq(evecR_exp, evecR / evecR(2,1), 1e-15, "evecR")
    end block

    call tc%finish()

  contains

    subroutine fact_Az(Ze, ctrans)
      !! Factorize Ze*B-A
      complex,           intent(in) :: Ze
      logical, optional, intent(in) :: ctrans
        !! factorize (Ze*B-A)^H. default: .false.

      integer :: i

      ! Az^H can be solved using factorization of Az
      if (present(ctrans)) then
        if (ctrans) return
      end if

      do i = 1, 10
        Azinv(i,i) = 1/(Ze - A(i,i))
      end do
    end subroutine

    subroutine solve_Az(xb, ctrans)
      !! Solve (Ze*B-A)*x=b
      complex,           intent(inout) :: xb(:,:)
        !! input:  b
        !! output: x
      logical, optional, intent(in)    :: ctrans
        !! solve (Ze*B-A)^H*x=b. default: .false.

      logical :: ctrans_

      ctrans_ = .false.
      if (present(ctrans)) ctrans_ = ctrans

      if (ctrans_) then
        xb = matmul(conjg(Azinv), xb)
      else
        xb = matmul(      Azinv,  xb)
      end if
    end subroutine

    subroutine mulvec_A(x, y, ctrans)
      !! Matrix-vector multiplication: y <- A*x
      complex,           intent(in)  :: x(:,:)
        !! input vector x
      complex,           intent(out) :: y(:,:)
        !! output vector y
      logical, optional, intent(in)  :: ctrans
        !! y <- A^H*x. default: .false.

      logical :: ctrans_

      ctrans_ = .false.
      if (present(ctrans)) ctrans_ = ctrans

      if (ctrans_) then
        y = matmul(conjg(A), x)
      else
        y = matmul(      A,  x)
      end if
    end subroutine
  end subroutine

end module

#endif
