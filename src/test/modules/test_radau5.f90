#include "../../util/macro.f90.inc"

module test_radau5_m

  use test_case_m
  use radau5_m
  use math_m

  implicit none

  private
  public test_radau5

contains

  subroutine test_radau5()
    type(test_case) :: tc

    print "(A)", "test_radau5"
    call tc%init("radau5")

    ! exponential
    block
      real, parameter   :: tol = 1e-10
      real              :: x0, x1, U0(1), P(1)
      integer           :: i, nsmp
      real, allocatable :: xsmp(:), Usmp(:,:), dUsmpdU0(:,:,:), dUsmpdP(:,:,:)
      type(ode_options) :: opt
      type(ode_result)  :: res

      call opt%init(1)

      x0    = 0.5
      x1    = 2.0
      U0(1) = 1.5
      P( 1) = 0.1

      ! samples
      nsmp = 5
      allocate (xsmp(nsmp), source = linspace(x0, x1, nsmp))

      call radau5(test_exp, x0, x1, xsmp, U0, P, opt, res)

      ! expected sample values
      allocate (Usmp(    1  ,nsmp))
      allocate (dUsmpdU0(1,1,nsmp))
      allocate (dUsmpdP( 1,1,nsmp))
      do i = 1, nsmp
        Usmp(    1  ,i) = U0(1) * exp(P(1) * (xsmp(i) - x0))
        dUsmpdU0(1,1,i) = exp(P(1) * (xsmp(i) - x0))
        dUsmpdP( 1,1,i) = (xsmp(i) - x0) * U0(1) * exp(P(1) * (xsmp(i) - x0))
      end do

      call tc%assert_eq(Usmp,     res%Usmp,     tol, "exp: Usmp"    )
      call tc%assert_eq(dUsmpdU0, res%dUsmpdU0, tol, "exp: dUsmpdU0")
      call tc%assert_eq(dUsmpdP,  res%dUsmpdP,  tol, "exp: dUsmpdP" )
    end block

    ! sine
    block
      real, parameter   :: tol = 1e-10
      real              :: x0, x1, U0(2), P(1)
      real              :: A, B, dAdU0(2), dAdP, dBdU0(2), dBdP
      integer           :: i, nsmp
      real, allocatable :: xsmp(:), Usmp(:,:), dUsmpdU0(:,:,:), dUsmpdP(:,:,:)
      type(ode_options) :: opt
      type(ode_result)  :: res

      call opt%init(2)

      x0    = 0.0
      x1    = -2.5 ! backwards
      U0(1) = 0.2
      U0(2) = -0.7
      P( 1) = 2.1

      ! samples
      nsmp = 7
      allocate (xsmp(nsmp), source = linspace(x1, x0, nsmp))

      call radau5(test_sin, x0, x1, xsmp, U0, P, opt, res)

      ! expected values ( U(x) = A * sin(P(1) * x) + B * cos(P(1) * x) )
      allocate (Usmp(    2  ,nsmp))
      allocate (dUsmpdU0(2,2,nsmp))
      allocate (dUsmpdP( 2,1,nsmp))
      A             = U0(1) * sin(P(1) * x0) + U0(2) * cos(P(1) * x0) / P(1)
      dAdU0(1)      = sin(P(1) * x0)
      dAdU0(2)      = cos(P(1) * x0) / P(1)
      dAdP          = U0(1) * cos(P(1) * x0) * x0 - U0(2) * sin(P(1) * x0) * x0 / P(1) - U0(2) * cos(P(1) * x0) / P(1)**2
      B             = U0(1) * cos(P(1) * x0) - U0(2) * sin(P(1) * x0) / P(1)
      dBdU0(1)      =  cos(P(1) * x0)
      dBdU0(2)      = -sin(P(1) * x0) / P(1)
      dBdP          = -U0(1) * sin(P(1) * x0) * x0 - U0(2) * cos(P(1) * x0) * x0 / P(1) + U0(2) * sin(P(1) * x0) / P(1)**2
      do i = 1, nsmp
        Usmp(    1  ,i) = A        * sin(P(1) * xsmp(i)) + B        * cos(P(1) * xsmp(i))
        Usmp(    2  ,i) = A * P(1) * cos(P(1) * xsmp(i)) - B * P(1) * sin(P(1) * xsmp(i))
        dUsmpdU0(1,:,i) = dAdU0 *        sin(P(1) * xsmp(i)) + dBdU0        * cos(P(1) * xsmp(i))
        dUsmpdU0(2,:,i) = dAdU0 * P(1) * cos(P(1) * xsmp(i)) - dBdU0 * P(1) * sin(P(1) * xsmp(i))
        dUsmpdP( 1,1,i) = dAdP * sin(P(1) * xsmp(i)) + A * cos(P(1) * xsmp(i)) * xsmp(i) + dBdP * cos(P(1) * xsmp(i)) - B * sin(P(1) * xsmp(i)) * xsmp(i)
        dUsmpdP( 2,1,i) = dAdP * P(1) * cos(P(1) * xsmp(i)) + A * cos(P(1) * xsmp(i)) - A * P(1) * sin(P(1) * xsmp(i)) * xsmp(i)  &
                      - dBdP * P(1) * sin(P(1) * xsmp(i)) - B * sin(P(1) * xsmp(i)) - B * P(1) * cos(P(1) * xsmp(i)) * xsmp(i)
      end do

      call tc%assert_eq(Usmp,     res%Usmp,     tol, "sin: Usmp"    )
      call tc%assert_eq(dUsmpdU0, res%dUsmpdU0, tol, "sin: UsmpdU0" )
      call tc%assert_eq(dUsmpdP,  res%dUsmpdP,  tol, "sin: UsmpdP"  )
    end block

    call tc%finish()
  end subroutine

  subroutine test_exp(x, U, P, f, dfdU, dfdP)
    real,              intent(in)  :: x
      !! x coordinate
    real,              intent(in)  :: U(:)
      !! state (1)
    real,              intent(in)  :: P(:)
      !! parameters (1)
    real, optional,    intent(out) :: f(:)
      !! output dU/dx (1)
    real, optional,    intent(out) :: dfdU(:,:)
      !! output derivatives of f wrt U (1,1)
    real, optional,    intent(out) :: dfdP(:,:)
      !! output derivatives of f wrt P (1,1)

    IGNORE(x)

    if (present(f   )) f(1) = P(1) * U(1)
    if (present(dfdU)) dfdU(1,1) = P(1)
    if (present(dfdP)) dfdP(1,1) = U(1)
  end subroutine

  subroutine test_sin(x, U, P, f, dfdU, dfdP)
    real,              intent(in)  :: x
      !! x coordinate
    real,              intent(in)  :: U(:)
      !! state (2)
    real,              intent(in)  :: P(:)
      !! parameters (1)
    real, optional,    intent(out) :: f(:)
      !! output dU/dx (2)
    real, optional,    intent(out) :: dfdU(:,:)
      !! output derivatives of f wrt U (2,2)
    real, optional,    intent(out) :: dfdP(:,:)
      !! output derivatives of f wrt P (2,1)

    IGNORE(x)

    if (present(f   )) then
      f(1) =             U(2)
      f(2) = - P(1)**2 * U(1)
    end if
    if (present(dfdU)) then
      dfdU(1,1) = 0
      dfdU(1,2) = 1
      dfdU(2,1) = - P(1)**2
      dfdU(2,2) = 0
    end if
    if (present(dfdP)) then
      dfdP(1,1) = 0
      dfdP(2,1) = - 2 * P(1) * U(1)
    end if
  end subroutine

end module
