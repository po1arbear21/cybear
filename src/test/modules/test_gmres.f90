module test_gmres_m

  use gmres_m
  use matrix_m
  use test_case_m
  use matop_m,        only: single_matop_real

  implicit none

  private
  public test_gmres

contains

  subroutine test_gmres
    integer                 :: n
    real, allocatable       :: x(:), b(:), x0(:), x_exp(:)
    type(band_real)         :: Ab
    type(gmres_options)     :: opts
    type(single_matop_real) :: Aop
    type(test_case)         :: tc

    print "(A)", "test_gmres"
    call tc%init("gmres")

    ! LES: Ax=b
    !
    ! A: tridiagonal, diagonally dominant matrix
    n = 10
    call Ab%init(n, 1)
    call Ab%set_diag(-1.0, k=-1)
    call Ab%set_diag(10.0      )
    call Ab%set_diag(-1.0, k= 1)
    call Aop%init(Ab)

    ! rhs, sol, expected sol
    b = [9,8,8,8,8,8,8,8,8,9]
    allocate (x_exp(n), source=1.0)
    allocate (x(n), x0(n))

    ! init sol: expected solution and some noise
    x0 = x_exp * (1 + 1e-2 * [1, 5, -5, 10, -3, 1, 5, -5, 10, -3])

    !
    ! test: w/o preconditioner
    !
    block
      x = x0
      call opts%init( x, b)
      call opts%check(x, b)
      opts%ipar( 9) = 1
      opts%ipar(10) = 0
      opts%ipar(12) = 1
      opts%dpar( 1) = 1e-13
      call gmres(opts, b, Aop, x)

      call tc%assert_eq(x_exp, x, 1e-13, "gmres: w/o preconditioner")
    end block

    !
    ! test: w/ preconditioner
    !
    block
      type(band_real)         :: prec_b
      type(single_matop_real) :: prec_op


      ! preconditioner: diagonal part of A
      !   P = diag(A)
      ! works b.c. A is a diagonally dominant matrix
      call prec_b%init(n, 0)
      call prec_b%set_diag(10.0)
      call prec_b%factorize
      call prec_op%init(prec_b, inv=.true.)

      x = x0
      call opts%init( x, b)
      call opts%check(x, b)
      opts%ipar( 9) = 1
      opts%ipar(10) = 0
      opts%ipar(11) = 1
      opts%ipar(12) = 1
      opts%dpar( 1) = 1e-13
      call gmres(opts, b, Aop, x, precon=prec_op)

      call tc%assert_eq(x_exp, x, 1e-13, "gmres: w/ preconditioner")
    end block

    call tc%finish
  end subroutine

end module
