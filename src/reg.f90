module reg_m

  use error_m,         only: program_error
  use distributions_m, only: inv_fermi_dirac_integral_1h
  use fukushima_m,     only: fd1h, fdm1h, fdm3h, fdm5h, fdm7h, fdm9h, dfdm9h
  use newton_m,        only: newton1D_opt, newton1D
  use ieee_arithmetic, only: ieee_is_finite

  implicit none

  private
  public fermi_dirac_integral_1h_reg, inv_fermi_dirac_integral_1h_reg

  real, parameter :: A = 1e-15
  real, parameter :: B = 0.02

contains

  subroutine fermi_dirac_integral_1h_reg(eta, F, dF1, dF2, dF3, dF4, dF5, dF6)
    !! regularized fermi-dirac integral for j = 1/2 and derivatives
    real,           intent(in)  :: eta
      !! argument
    real, optional, intent(out) :: F
      !! output distribution function
    real, optional, intent(out) :: dF1
      !! output first derivative of F wrt eta
    real, optional, intent(out) :: dF2
      !! output second derivative of F wrt eta
    real, optional, intent(out) :: dF3
      !! output third derivative of F wrt eta
    real, optional, intent(out) :: dF4
      !! output fourth derivative of F wrt eta
    real, optional, intent(out) :: dF5
      !! output fifth derivative of F wrt eta
    real, optional, intent(out) :: dF6
      !! output sixth derivative of F wrt eta

    real, parameter :: G0 = 1.0 / gamma(1.5)
    real, parameter :: G1 = 1.0 / gamma(0.5)
    real, parameter :: G2 = 1.0 / gamma(-0.5)
    real, parameter :: G3 = 1.0 / gamma(-1.5)
    real, parameter :: G4 = 1.0 / gamma(-2.5)
    real, parameter :: G5 = 1.0 / gamma(-3.5)

    if (present(f  )) F   = (fd1h( eta) + A*     fd1h(B*eta)) * G0
    if (present(dF1)) dF1 = (fdm1h(eta) + A*B*   fdm1h(B*eta)) * G1
    if (present(dF2)) dF2 = (fdm3h(eta) + A*B**2*fdm3h(B*eta)) * G2
    if (present(dF3)) dF3 = (fdm5h(eta) + A*B**3*fdm5h(B*eta)) * G3
    if (present(dF4)) dF4 = (fdm7h(eta) + A*B**4*fdm7h(B*eta)) * G4
    if (present(dF5)) dF5 = (fdm9h(eta) + A*B**5*fdm9h(B*eta)) * G5
    if (present(dF6)) dF6 = (dfdm9h(eta) + A*B**6*dfdm9h(B*eta)) * G5
  end subroutine

  subroutine inv_fermi_dirac_integral_1h_reg(Fd, eta, deta)
    !! inverse of regularized fermi-dirac integral for j = 1/2 and derivative
    real,           intent(in)  :: Fd
      !! regularized fermi-dirac integral value (> 0)
    real, optional, intent(out) :: eta
      !! output argument of inverse fermi-dirac integral
    real, optional, intent(out) :: deta
      !! output derivative of eta wrt Fd

    real :: p(0), x0, dum, tmp1, tmp2
    type(newton1D_opt) :: nopt

    if (Fd<6.6e-11) then
      x0 = 1/B * log(Fd/A)
    else
      call inv_fermi_dirac_integral_1h(Fd, x0, dum)
    end if

    call nopt%init()
    call newton1D(fun, p, nopt, x0, eta)

    call fermi_dirac_integral_1h_reg(eta, tmp1, tmp2)
    if(present(deta)) deta = 1/tmp2
    if(.not. (ieee_is_finite(eta) .and. ieee_is_finite(deta))) call program_error("not finite eta in inv_fd")

  contains

    subroutine fun(x, p, f, dfdx, dfdp)
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

      real :: val

      call fermi_dirac_integral_1h_reg(x, val, dfdx)
      f = val - Fd
    end subroutine

  end subroutine

end module
