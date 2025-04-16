module fermi_m

  use distributions_m, only: fermi_dirac_integral_m1h, fermi_dirac_integral_1h, inv_fermi_dirac_integral_1h
  use error_m,         only: program_error
  use fukushima_m,     only: fd1h, fdm1h, fdm3h, fdm5h, fdm7h, fdm9h, dfdm9h
  use ieee_arithmetic, only: ieee_is_finite
  use newton_m,        only: newton1D, newton1D_opt

  implicit none

  real, parameter :: GAM     = sqrt(0.125)
  real, parameter :: ETA_SMALL = -16.0
  real, parameter :: ETA_TINY  = -36.0
  real, parameter :: ETA_MICRO = -50.0
  real, parameter :: ETA_REF   = -5000.0
  real, parameter :: F_SMALL   = 1.0 / (exp(-ETA_SMALL) + GAM)
  real, parameter :: F_TINY    = exp(ETA_TINY)
  real, parameter :: F_MICRO   = 1.0 / (exp(-ETA_MICRO) + GAM)
  real, parameter :: F_REF     = 1e-50
  real, parameter :: BETA      = log(F_REF / F_MICRO) / (ETA_REF - ETA_MICRO)

  real, parameter :: A = 1e-12 !1e-15!1e-8
  real, parameter :: B = 0.001 !0.004!0.003

contains

  subroutine fermi_dirac_integral_1h_reg_old(eta, f, dfdeta)
    !! fermi-dirac integral for j = 1/2 with regularization to avoid underflow
    real, intent(in)  :: eta
      !! argument
    real, intent(out) :: f
      !! value of fermi-dirac integral
    real, intent(out) :: dfdeta
      !! output derivative of f wrt eta

    real :: e

    if (eta < ETA_MICRO) then
      f      = F_MICRO * exp(BETA * (eta - ETA_MICRO))
      dfdeta = f * BETA
    elseif (eta > ETA_SMALL) then
      call fermi_dirac_integral_1h(eta, f, dfdeta)
    elseif (eta > ETA_TINY) then
      e      = exp(-eta)
      f      = 1.0 / (e + GAM)
      dfdeta = e / (e + GAM)**2
    else
      f      = exp(eta)
      dfdeta = f
    end if
  end subroutine

  subroutine inv_fermi_dirac_integral_1h_reg_old(f, eta, detadf)
    !! inverse fermi-dirac integral for j = 1/2 with regularization to avoid underflow
    real, intent(in)  :: f
      !! value of fermi-dirac integral
    real, intent(out) :: eta
      !! argument
    real, intent(out) :: detadf
      !! output derivative of eta wrt f

    if (f < F_MICRO) then
      eta    = ETA_MICRO + log(f / F_MICRO) / BETA
      detadf = 1.0 / (BETA * f)
    elseif (f > F_SMALL) then
      call inv_fermi_dirac_integral_1h(f, eta, detadf)
    elseif (f > F_TINY) then
      eta    = - log(1 / f - GAM)
      detadf = 1.0 / (f * (1 - GAM * f))
    else
      eta    = log(f)
      detadf = 1 / f
    end if
  end subroutine

  subroutine fermi_dirac_integral_m1h_reg(eta, f, dfdeta)
    !! derivative of fermi_dirac_integral_1h_reg
    real, intent(in)  :: eta
      !! argument
    real, intent(out) :: f
      !! value of fermi-dirac integral
    real, intent(out) :: dfdeta
      !! output derivative of y wrt x

    real :: e

    if (f < F_MICRO) then
      f      = F_MICRO * exp(BETA * (eta - ETA_MICRO)) * BETA
      dfdeta = f * BETA
    elseif (eta > ETA_SMALL) then
      call fermi_dirac_integral_m1h(eta, f, dfdeta)
    elseif (eta > ETA_TINY) then
      e      = exp(-eta)
      f      = e / (e + GAM)**2
      dfdeta = e * (e - GAM) / (e + GAM)**3
    else
      f      = exp(eta)
      dfdeta = f
    end if
  end subroutine

  subroutine fermi_dirac_generation_reg(eta, g, dgdeta)
    !! g ~ F_12 * exp(-eta) with regularization (~1 for small eta instead of large values)
    real, intent(in)  :: eta
      !! argument
    real, intent(out) :: g
      !! output value of generation factor
    real, intent(out) :: dgdeta
      !! output derivative of g wrt eta

    real :: e, f, dfdeta

    if (eta > ETA_SMALL) then
      call fermi_dirac_integral_1h(eta, f, dfdeta)
      e      = exp(-eta)
      g      = f * e
      dgdeta = dfdeta * e - g
    elseif (eta > ETA_TINY) then
      e      = GAM * exp(eta)
      g      = 1.0 / (1.0 + e)
      dgdeta = - e * g**2
    else
      g      = 1.0
      dgdeta = 0
    end if
  end subroutine

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
