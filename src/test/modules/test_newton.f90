module test_newton_m

  use matrix_m,    only: dense_real, matrix_real
  use newton_m,    only: newton, newton_opt, newton1D, newton1D_opt
  use test_case_m, only: test_case

  implicit none

  private
  public test_newton

contains

  subroutine test_newton()
    type(test_case) :: tc

    call tc%init("newton")

    ! test newton1D
    block
      real               :: x, dxdp(2), xe, dxedp(2)
      type(newton1D_opt) :: opt

      real,    parameter :: p(2)    = [3, 1]

      ! expected values for f = (x - rp1*ip1) * (4*x**2 + rp2)
      xe       = - 2.0 * p(1)
      dxedp(1) = - 2.0
      dxedp(2) = 0

      call opt%init(atol = 1e-16, rtol = 1e-14, xmin = -10*abs(xe), xmax = 10*abs(xe), log = .true., msg = "poly1D: ")
      call newton1D(poly1D, p, opt, -1.0, x, dxdp = dxdp)

      call tc%assert_eq(   xe,    x, opt%rtol, opt%atol, "newton1D: x"   )
      call tc%assert_eq(dxedp, dxdp, opt%rtol, opt%atol, "newton1D: dxdp")
    end block

    ! test multidimensional newton
    call test_poly2D(tc)

    call tc%finish()
  end subroutine

  subroutine test_poly2D(tc)
    type(test_case), intent(inout) :: tc

    real                     :: x(2), dxdp(2,2), xe(2), dxedp(2,2)
    type(dense_real), target :: jaco, prec
    type(newton_opt)         :: opt

    ! init matrices
    call jaco%init(2)
    call prec%init(2)

    ! expected values
    xe(1) = 1
    xe(2) = 5
    dxedp(1,1) = 1
    dxedp(1,2) = 0
    dxedp(2,1) = 7
    dxedp(2,2) = 1

    ! direct solver
    call opt%init(2, atol = 1e-15, rtol = 1e-14, log = .true., msg = "poly2D: dir sol: ")
    call newton(poly2D, [1.0, 4.0], opt, [0.5, 6.0], x, dxdp = dxdp)

    call tc%assert_eq(   xe,    x, opt%rtol(1), opt%atol(1), "newton: dir sol: x"   )
    call tc%assert_eq(dxedp, dxdp, opt%rtol(1), opt%atol(1), "newton: dir sol: dxdp")

    ! iterative solver with preconditioner
    call opt%init(2, atol = 1e-15, rtol = 1e-14, it_solver = .true., log = .true., &
      &           msg = "poly2D: it sol: ")
    call newton(poly2D, [1.0, 4.0], opt, [0.5, 6.0], x, dxdp = dxdp)

    call tc%assert_eq(   xe,    x, opt%rtol(1), opt%atol(1),    "newton: it sol: x"   )
    call tc%assert_eq(dxedp, dxdp, opt%rtol(1), opt%atol(1)*10, "newton: it sol: dxdp")  ! fixme gmres doesnt achieve tolerance?

    ! destruct matrices
    call jaco%destruct()
    call prec%destruct()

  contains

    subroutine poly2D(x, p, f, dfdx, dfdx_prec, dfdp)
      real,                                  intent(in)  :: x(:)
        !! arguments
      real,                                  intent(in)  :: p(:)
        !! parameters
      real,               optional,          intent(out) :: f(:)
        !! output function values
      class(matrix_real), optional, pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
      class(matrix_real), optional, pointer, intent(out) :: dfdx_prec
        !! optional output pointer to preconditioner jacobian of f wrt x
      real,               optional,          intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p

      ! set residual f
      if (present(f)) then
        f(1) = x(1)**3 + p(2)*x(1) - x(2)
        f(2) = p(1)*x(1)**2*x(2) + p(1)*p(2)*x(2) - x(2)**2
      end if

      ! set jacobian dfdx
      if (present(dfdx)) then
        jaco%d(1,1) = 3*x(1)**2 + p(2)
        jaco%d(1,2) = -1
        jaco%d(2,1) = 2*p(1)*x(1)*x(2)
        jaco%d(2,2) = p(1)*x(1)**2 + p(1)*p(2) - 2*x(2)
        dfdx => jaco
      end if

      ! set preconditioner dfdx_prec
      if (present(dfdx_prec)) then
        prec%d(1,1) = 2*x(1)**2 + p(2)
        prec%d(1,2) = -3
        prec%d(2,1) = p(1)*x(1)*x(2)
        prec%d(2,2) = p(1)*x(1)**2 + p(1)*p(2) - 3*x(2)
        dfdx_prec => prec
      end if

      ! set dfdp
      if (present(dfdp)) then
        dfdp(1,1) = 0
        dfdp(1,2) = x(1)
        dfdp(2,1) = x(1)**2*x(2) + p(2)*x(2)
        dfdp(2,2) = p(1)*x(2)
      end if
    end subroutine
  end subroutine

  subroutine poly1D(x, p, f, dfdx, dfdp)
    real,              intent(in)  :: x
      !! argument
    real,              intent(in)  :: p(:)
      !! parameters
    real,              intent(out) :: f
      !! output function value
    real,    optional, intent(out) :: dfdx
      !! optional output derivative of f wrt x
    real,    optional, intent(out) :: dfdp(:)
      !! optional output derivatives of f wrt p

    f = (x + 2*p(1)) * (4*x**2 + p(2))
    if (present(dfdx)) dfdx = 4*x*(3*x + 4*p(1)) + p(2)
    if (present(dfdp)) then
      dfdp(1) = 2*(4*x**2 + p(2))
      dfdp(2) = x + 2*p(1)
    end if
  end subroutine

end module
