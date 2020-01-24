#include "../../util/macro.f90.inc"

module test_radau5_m
  use test_case_m
  use radau5_m
  use util_m
  implicit none

contains

  subroutine test_radau5()
    type(test_case) :: tc

    print "(1A)", "test_radau5"
    call tc%init("radau5")

    ! exponential
    block
      real, parameter   :: tol = 1e-10
      real              :: x0, x1, U0(1), P(1)
      real              :: U1(1), dU1dU0(1,1), dU1dP(1,1), UA(1), dUAdU0(1,1), dUAdP(1,1)
      integer           :: i, ns
      real, allocatable :: xs(:), Us(:,:), dUsdU0(:,:,:), dUsdP(:,:,:)
      type(ode_options) :: opt
      type(ode_result)  :: res

      call opt%init(1)

      x0    = 0.5
      x1    = 2.0
      U0(1) = 1.5
      P( 1) = 0.1

      ! samples
      ns = 5
      allocate (xs(ns), source = linspace(x0, x1, ns))

      call radau5(test_exp, x0, x1, U0, P, opt, res, xs = xs)

      ! expected values (U(x) = U0(1) * exp(P(1) * (x - x0))
      U1(    1  ) = U0(1) * exp(P(1) * (x1 - x0))
      dU1dU0(1,1) = exp(P(1) * (x1 - x0))
      dU1dP( 1,1) = (x1 - x0) * U0(1) * exp(P(1) * (x1 - x0))
      UA(    1  ) = U0(1) * (exp(P(1) * (x1 - x0)) - 1.0) / (P(1) * (x1 - x0))
      dUAdU0(1,1) = (exp(P(1) * (x1 - x0)) - 1.0) / (P(1) * (x1 - x0))
      dUAdP( 1,1) = U0(1) * (exp(P(1) * (x1 - x0)) / (P(1)) - (exp(P(1) * (x1 - x0)) - 1.0) / (P(1)**2 * (x1 - x0)))

      ! expected sample values
      allocate (Us(    1  ,ns))
      allocate (dUsdU0(1,1,ns))
      allocate (dUsdP( 1,1,ns))
      do i = 1, ns
        Us(    1  ,i) = U0(1) * exp(P(1) * (xs(i) - x0))
        dUsdU0(1,1,i) = exp(P(1) * (xs(i) - x0))
        dUsdP( 1,1,i) = (xs(i) - x0) * U0(1) * exp(P(1) * (xs(i) - x0))
      end do

      call tc%assert_eq(Us,     res%Us,     tol, "exp: Us"    )
      call tc%assert_eq(dUsdU0, res%dUsdU0, tol, "exp: dUsdU0")
      call tc%assert_eq(dUsdP,  res%dUsdP,  tol, "exp: dUsdP" )
      call tc%assert_eq(U1,     res%U1,     tol, "exp: U1"    )
      call tc%assert_eq(dU1dU0, res%dU1dU0, tol, "exp: dU1dU0")
      call tc%assert_eq(dU1dP,  res%dU1dP,  tol, "exp: dU1dP" )
      call tc%assert_eq(UA,     res%UA,     tol, "exp: UA"    )
      call tc%assert_eq(dUAdU0, res%dUAdU0, tol, "exp: dUAdU0")
      call tc%assert_eq(dUAdP,  res%dUAdP,  tol, "exp: dUAdP" )
    end block

    ! sine
    block
      real, parameter   :: tol = 1e-10
      real              :: x0, x1, U0(2), P(1)
      real              :: A, B, dAdU0(2), dAdP, dBdU0(2), dBdP
      real              :: U1(2), dU1dU0(2,2), dU1dP(2,1), UA(2), dUAdU0(2,2), dUAdP(2,1)
      integer           :: i, ns
      real, allocatable :: xs(:), Us(:,:), dUsdU0(:,:,:), dUsdP(:,:,:)
      type(ode_options) :: opt
      type(ode_result)  :: res

      call opt%init(2)

      x0    = 0.0
      x1    = -2.5 ! backwards
      U0(1) = 0.2
      U0(2) = -0.7
      P( 1) = 2.1

      ! samples
      ns = 7
      allocate (xs(ns), source = linspace(x0, x1, ns))

      call radau5(test_sin, x0, x1, U0, P, opt, res, xs = xs)

      ! expected values ( U(x) = A * sin(P(1) * x) + B * cos(P(1) * x) )
      A             = U0(1) * sin(P(1) * x0) + U0(2) * cos(P(1) * x0) / P(1)
      dAdU0(1)      = sin(P(1) * x0)
      dAdU0(2)      = cos(P(1) * x0) / P(1)
      dAdP          = U0(1) * cos(P(1) * x0) * x0 - U0(2) * sin(P(1) * x0) * x0 / P(1) - U0(2) * cos(P(1) * x0) / P(1)**2
      B             = U0(1) * cos(P(1) * x0) - U0(2) * sin(P(1) * x0) / P(1)
      dBdU0(1)      =  cos(P(1) * x0)
      dBdU0(2)      = -sin(P(1) * x0) / P(1)
      dBdP          = -U0(1) * sin(P(1) * x0) * x0 - U0(2) * cos(P(1) * x0) * x0 / P(1) + U0(2) * sin(P(1) * x0) / P(1)**2

      U1(1)       = A        * sin(P(1) * x1) + B        * cos(P(1) * x1)
      U1(2)       = A * P(1) * cos(P(1) * x1) - B * P(1) * sin(P(1) * x1)
      dU1dU0(1,:) = dAdU0 *        sin(P(1) * x1) + dBdU0        * cos(P(1) * x1)
      dU1dU0(2,:) = dAdU0 * P(1) * cos(P(1) * x1) - dBdU0 * P(1) * sin(P(1) * x1)
      dU1dP(1,1)  = dAdP * sin(P(1) * x1) + A * cos(P(1) * x1) * x1 + dBdP * cos(P(1) * x1) - B * sin(P(1) * x1) * x1
      dU1dP(2,1)  = dAdP * P(1) * cos(P(1) * x1) + A * cos(P(1) * x1) - A * P(1) * sin(P(1) * x1) * x1  &
                  - dBdP * P(1) * sin(P(1) * x1) - B * sin(P(1) * x1) - B * P(1) * cos(P(1) * x1) * x1
      UA(1)       = -(    A*(cos(P(1)*x0) - cos(P(1)*x1)) -     B*( sin(P(1)*x0)    - sin(P(1)*x1)   ))/(P(1) * (x0 - x1))
      UA(2)       =  (    B*(cos(P(1)*x0) - cos(P(1)*x1)) +     A*( sin(P(1)*x0)    - sin(P(1)*x1)   ))/(x0 - x1)
      dUAdU0(1,:) = -(dAdU0*(cos(P(1)*x0) - cos(P(1)*x1)) - dBdU0*( sin(P(1)*x0)    - sin(P(1)*x1)   ))/(P(1) * (x0 - x1))
      dUAdU0(2,:) =  (dBdU0*(cos(P(1)*x0) - cos(P(1)*x1)) + dAdU0*( sin(P(1)*x0)    - sin(P(1)*x1)   ))/(x0 - x1)
      dUAdP(1,1)  = -( dAdP*(cos(P(1)*x0) - cos(P(1)*x1)) +     A*(-sin(P(1)*x0)*x0 + sin(P(1)*x1)*x1)                     &
                    -  dBdP*(sin(P(1)*x0) - sin(P(1)*x1)) -     B*( cos(P(1)*x0)*x0 - cos(P(1)*x1)*x1))/(P(1) * (x0 - x1)) &
                    +(    A*(cos(P(1)*x0) - cos(P(1)*x1)) -     B*( sin(P(1)*x0)    - sin(P(1)*x1)   ))/(P(1)**2 * (x0 - x1))
      dUAdP(2,1)  =  ( dBdP*(cos(P(1)*x0) - cos(P(1)*x1)) +     B*(-sin(P(1)*x0)*x0 + sin(P(1)*x1)*x1) &
                     + dAdP*(sin(P(1)*x0) - sin(P(1)*x1)) +     A*( cos(P(1)*x0)*x0 - cos(P(1)*x1)*x1))/(x0 - x1)

      ! expected sample values
      allocate (Us(    2  ,ns))
      allocate (dUsdU0(2,2,ns))
      allocate (dUsdP( 2,1,ns))
      do i = 1, ns
        Us(    1  ,i) = A        * sin(P(1) * xs(i)) + B        * cos(P(1) * xs(i))
        Us(    2  ,i) = A * P(1) * cos(P(1) * xs(i)) - B * P(1) * sin(P(1) * xs(i))
        dUsdU0(1,:,i) = dAdU0 *        sin(P(1) * xs(i)) + dBdU0        * cos(P(1) * xs(i))
        dUsdU0(2,:,i) = dAdU0 * P(1) * cos(P(1) * xs(i)) - dBdU0 * P(1) * sin(P(1) * xs(i))
        dUsdP( 1,1,i) = dAdP * sin(P(1) * xs(i)) + A * cos(P(1) * xs(i)) * xs(i) + dBdP * cos(P(1) * xs(i)) - B * sin(P(1) * xs(i)) * xs(i)
        dUsdP( 2,1,i) = dAdP * P(1) * cos(P(1) * xs(i)) + A * cos(P(1) * xs(i)) - A * P(1) * sin(P(1) * xs(i)) * xs(i)  &
                      - dBdP * P(1) * sin(P(1) * xs(i)) - B * sin(P(1) * xs(i)) - B * P(1) * cos(P(1) * xs(i)) * xs(i)
      end do

      call tc%assert_eq(Us,     res%Us,     tol, "sin: Us"    )
      call tc%assert_eq(dUsdU0, res%dUsdU0, tol, "sin: UsdU0" )
      call tc%assert_eq(dUsdP,  res%dUsdP,  tol, "sin: UsdP"  )
      call tc%assert_eq(U1,     res%U1,     tol, "sin: U1"    )
      call tc%assert_eq(dU1dU0, res%dU1dU0, tol, "sin: dU1dU0")
      call tc%assert_eq(dU1dP,  res%dU1dP,  tol, "sin: dU1dP" )
      call tc%assert_eq(UA,     res%UA,     tol, "sin: UA"    )
      call tc%assert_eq(dUAdU0, res%dUAdU0, tol, "exp: dUAdU0")
      call tc%assert_eq(dUAdP,  res%dUAdP,  tol, "exp: dUAdP" )
    end block

    call tc%finish()
  end subroutine test_radau5

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
