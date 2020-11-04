module test_gmres_m

  use gmres_m
  use matrix_m
  use test_case_m
  use matop_m,        only: single_matop_real, single_matop_cmplx

  implicit none

  private
  public test_gmres

contains

  subroutine test_gmres()
    type(test_case) :: tc

    print "(A)", "test_gmres"
    call tc%init("gmres")

    call test_gmres_real(tc)
    call test_gmres_cmplx(tc)

    call tc%finish()
  end subroutine

  subroutine test_gmres_real(tc)
    !! tests gmres for real arguments

    type(test_case), intent(inout) :: tc

    integer                 :: n
    real, allocatable       :: x(:), b(:), x0(:), x_exp(:)
    type(band_real), target :: Ab
    type(gmres_options)     :: opts
    type(single_matop_real) :: A

    ! LES: Ax=b
    !
    ! A: tridiagonal, diagonally dominant matrix
    n = 10
    call Ab%init(n, 1)
    call Ab%set_diag(-1.0, k=-1)
    call Ab%set_diag(10.0      )
    call Ab%set_diag(-1.0, k= 1)
    call A%init(Ab)

    ! rhs, sol, expected sol
    b = [9,8,8,8,8,8,8,8,8,9]
    allocate (x_exp(n), source=1.0)
    allocate (x(n), x0(n))

    ! init sol
    x0 = 0

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
      call gmres(opts, b, A, x)

      call tc%assert_eq(x_exp, x, 1e-13, "gmres: w/o preconditioner")
    end block

    !
    ! test: w/ preconditioner
    !
    block
      type(band_real), target :: prec_b
      type(single_matop_real) :: prec

      ! preconditioner: diagonal part of A
      !   P = diag(A)
      ! works b.c. A is a diagonally dominant matrix
      call prec_b%init(n, 0)
      call prec_b%set_diag(10.0)
      call prec_b%factorize()
      call prec%init(prec_b, inv=.true.)

      x = x0
      call opts%init( x, b)
      call opts%check(x, b)
      opts%ipar( 9) = 1
      opts%ipar(10) = 0
      opts%ipar(11) = 1
      opts%ipar(12) = 1
      opts%dpar( 1) = 1e-13
      call gmres(opts, b, A, x, precon=prec)

      call tc%assert_eq(x_exp, x, 1e-13, "gmres: w/ preconditioner")
    end block
  end subroutine

  subroutine test_gmres_cmplx(tc)
    !! tests gmres for complex arguments

    type(test_case), intent(inout) :: tc

    complex, allocatable     :: x(:), b(:), x0(:), x_exp(:)
    integer                  :: n
    type(band_cmplx), target :: Ab
    type(gmres_options)      :: opts
    type(single_matop_cmplx) :: A

    ! LES: Ax=b
    !
    ! A: tridiagonal, diagonally dominant matrix
    ! matlab: A=diag((10-5i)*ones(10,1)) + diag((-1+2i)*ones(9,1),k=-1) + diag((-1+2i)*ones(9,1),k=1);
    n = 10
    call Ab%init(n, 1)
    call Ab%set_diag((-1,  2), k=-1)
    call Ab%set_diag((10, -5)      )
    call Ab%set_diag((-1,  2), k= 1)
    call A%init(Ab)

    ! rhs, sol, expected sol
    !     b := A*x_exp
    allocate (b(n), source=(19,22))
    b(1) = (27,21)
    b(n) = (27,21)
    allocate (x_exp(n), source=(2,3))
    allocate (x(n), x0(n))

    ! init sol: expected solution and some noise
    x0 = 0

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
      call gmres(opts, b, A, x)

      call tc%assert_eq(x_exp, x, 1e-13, "gmres cmplx: w/o preconditioner")
    end block

    !
    ! test: w/ preconditioner
    !
    block
      type(band_cmplx), target :: prec_b
      type(single_matop_cmplx) :: prec

      ! preconditioner: diagonal part of A
      !   P = diag(A)
      ! works b.c. A is a diagonally dominant matrix
      call prec_b%init(n, 0)
      call prec_b%set_diag((10, -5))
      call prec_b%factorize()
      call prec%init(prec_b, inv=.true.)

      x = x0
      call opts%init( x, b)
      call opts%check(x, b)
      opts%ipar( 9) = 1
      opts%ipar(10) = 0
      opts%ipar(11) = 1
      opts%ipar(12) = 1
      opts%dpar( 1) = 1e-13
      call gmres(opts, b, A, x, precon=prec)

      call tc%assert_eq(x_exp, x, 1e-13, "gmres cmplx: w/ preconditioner")
    end block
  end subroutine

end module
