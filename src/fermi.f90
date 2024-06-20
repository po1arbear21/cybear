module fermi_m

  use distributions_m, only: fermi_dirac_integral_m1h, fermi_dirac_integral_1h, inv_fermi_dirac_integral_1h

  implicit none

  real, parameter :: GAMMA     = sqrt(0.125)
  real, parameter :: ETA_SMALL = -16.0
  real, parameter :: ETA_TINY  = -36.0
  real, parameter :: ETA_MICRO = -25.0
  real, parameter :: ETA_REF   = -5000.0
  real, parameter :: F_SMALL   = 1.0 / (exp(-ETA_SMALL) + GAMMA)
  real, parameter :: F_TINY    = exp(ETA_TINY)
  real, parameter :: F_MICRO   = 1.0 / (exp(-ETA_MICRO) + GAMMA)
  real, parameter :: F_REF     = 1e-12
  real, parameter :: BETA      = log(F_REF / F_MICRO) / (ETA_REF - ETA_MICRO)

contains

  subroutine fermi_dirac_integral_1h_reg(eta, f, dfdeta)
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
      f      = 1.0 / (e + GAMMA)
      dfdeta = e / (e + GAMMA)**2
    else
      f      = exp(eta)
      dfdeta = f
    end if
  end subroutine

  subroutine inv_fermi_dirac_integral_1h_reg(f, eta, detadf)
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
      eta    = - log(1 / f - GAMMA)
      detadf = 1.0 / (f * (1 - GAMMA * f))
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
      f      = e / (e + GAMMA)**2
      dfdeta = e * (e - GAMMA) / (e + GAMMA)**3
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
      e      = GAMMA * exp(eta)
      g      = 1.0 / (1.0 + e)
      dgdeta = - e * g**2
    else
      g      = 1.0
      dgdeta = 0
    end if
  end subroutine

end module
