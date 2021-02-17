module test_newton_m
  use test_case_m
  use newton_m
  implicit none

  private
  public test_newton

contains

  subroutine test_newton()
    type(test_case) :: tc

    print "(A)", "test_newton"
    call tc%init("newton")

    block
      type(newton1D_opt) :: opt
      real :: x, dxdp(2), xe, dxedp(2)

      real,    parameter :: p(2)    = [3, 1]
      integer, parameter :: ipar(1) = [-2]

      ! expected values for f = (x - rp1*ip1) * (4*x**2 + rp2)
      xe       = p(1)*ipar(1)
      dxedp(1) = ipar(1)
      dxedp(2) = 0

      call opt%init(atol = 1e-16, rtol = 1e-14, xmin = -10*abs(xe), xmax = 10*abs(xe), log = .true., msg = "poly1D: ")
      call newton1D(poly1D, p, opt, -1.0, x, dxdp = dxdp, ipar = ipar)
      print *

      call tc%assert_eq(   xe,    x, opt%atol, "newton1D: x" )
      call tc%assert_eq(dxedp, dxdp, opt%atol, "newton1D: dxdp")
    end block

    block
      type(newton_opt) :: opt
      real :: x(2), dxdp(2,2), xe(2), dxedp(2,2)

      call opt%init(2, atol = 1e-16 * [1.0, 1.0], rtol = 1e-14 * [1.0, 1.0], log = .true., msg = "poly2D: ")
      call newton(poly2D, [1.0, 4.0], opt, [0.5, 6.0], x, dxdp = dxdp)
      print *

      ! expected values
      xe(1) = 1.0
      xe(2) = 5.0
      dxedp(1,1) = 1.0
      dxedp(1,2) = 0.0
      dxedp(2,1) = 7.0
      dxedp(2,2) = 1.0

      call tc%assert_eq(   xe,    x, opt%atol(1), "newton: x")
      call tc%assert_eq(dxedp, dxdp, opt%atol(1), "newton: dxdp")
    end block

    call tc%finish()
  end subroutine

  subroutine poly1D(x, p, f, ipar, dfdx, dfdp)
    real,              intent(in)  :: x
      !! argument
    real,              intent(in)  :: p(:)
      !! parameters
    real,              intent(out) :: f
      !! output function value
    integer, optional, intent(in)  :: ipar(:)
      !! integer paramters
    real,    optional, intent(out) :: dfdx
      !! optional output derivative of f wrt x
    real,    optional, intent(out) :: dfdp(:)
      !! optional output derivatives of f wrt p

    f = (x - p(1)*ipar(1)) * (4*x**2 + p(2))
    if (present(dfdx)) dfdx = 4*x*(3*x - 2*p(1)*ipar(1)) + p(2)
    if (present(dfdp)) then
      dfdp(1) = - ipar(1)*(4*x**2 + p(2))
      dfdp(2) = x - p(1)*ipar(1)
    end if
  end subroutine

  subroutine poly2D(x, p, f, dfdx, dfdp)
    real,                        intent(in) :: x(:)
      !! arguments
    real,                        intent(in) :: p(:)
      !! parameters
    real,                        intent(out) :: f(:)
      !! output function values
    class(matrix_real), pointer, intent(out) :: dfdx
      !! output pointer to jacobian of f wrt x
    real, optional,              intent(out) :: dfdp(:,:)
      !! optional output jacobian of f wrt p

    type(dense_real), target, save :: jaco

    f(1) = x(1)**3 + p(2)*x(1) - x(2)
    f(2) = p(1)*x(1)**2*x(2) + p(1)*p(2)*x(2) - x(2)**2

    if (.not. allocated(jaco%d)) then
      call jaco%init(2)
    end if

    jaco%d(1,1) = 3*x(1)**2 + p(2)
    jaco%d(1,2) = - 1.0
    jaco%d(2,1) = 2*p(1)*x(1)*x(2)
    jaco%d(2,2) = p(1)*x(1)**2 + p(1)*p(2) - 2*x(2)
    dfdx => jaco

    if (present(dfdp)) then
      dfdp(1,1) = 0.0
      dfdp(1,2) = x(1)
      dfdp(2,1) = x(1)**2*x(2) + p(2)*x(2)
      dfdp(2,2) = p(1)*x(2)
    end if
  end subroutine

end module
