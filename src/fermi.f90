module fermi_m

  use distributions_m, only: fermi_dirac_integral_m1h, fermi_dirac_integral_1h, inv_fermi_dirac_integral_1h

  implicit none

  real, parameter :: GAMMA     = sqrt(0.125)
  real, parameter :: ETA_SMALL = -16.0
  real, parameter :: ETA_TINY  = -36.0
  real, parameter :: ETA_MICRO = -50.0
  real, parameter :: ETA_REF   = -5000.0
  real, parameter :: F_SMALL   = 1.0 / (exp(-ETA_SMALL) + GAMMA)
  real, parameter :: F_TINY    = exp(ETA_TINY)
  real, parameter :: F_MICRO   = exp(ETA_MICRO)
  real, parameter :: F_REF     = 1e-50

contains

  subroutine fermi_dirac_integral_m1h_reg(eta, f, dfdeta)
    !! derivative of fermi_dirac_integral_1h_reg
    real, intent(in)  :: eta
      !! argument
    real, intent(out) :: f
      !! value of fermi-dirac integral
    real, intent(out) :: dfdeta
      !! output derivative of y wrt x

    real :: e

    if (eta > ETA_SMALL) then
      call fermi_dirac_integral_m1h(eta, f, dfdeta)
    elseif (eta > ETA_TINY) then
      e      = exp(-eta)
      f      = e / (e + GAMMA)**2
      dfdeta = e * (e - GAMMA) / (e + GAMMA)**3
    elseif (eta > ETA_MICRO) then
      e      = exp(eta)
      f      = e
      dfdeta = e
    else
      e      = (ETA_MICRO - log(F_REF)) / (ETA_MICRO - ETA_REF)
      f      = F_REF * exp((eta - ETA_REF) * e) * e
      dfdeta = f * e
    end if
  end subroutine

  subroutine fermi_dirac_integral_1h_reg(eta, f, dfdeta)
    !! fermi-dirac integral for j = 1/2 with regularization to avoid underflow
    real, intent(in)  :: eta
      !! argument
    real, intent(out) :: f
      !! value of fermi-dirac integral
    real, intent(out) :: dfdeta
      !! output derivative of f wrt eta

    real :: e

    if (eta > ETA_SMALL) then
      call fermi_dirac_integral_1h(eta, f, dfdeta)
    elseif (eta > ETA_TINY) then
      e      = exp(-eta)
      f      = 1.0 / (e + GAMMA)
      dfdeta = e / (e + GAMMA)**2
    elseif (eta > ETA_MICRO) then
      e      = exp(eta)
      f      = e
      dfdeta = e
    else
      e      = (ETA_MICRO - log(F_REF)) / (ETA_MICRO - ETA_REF)
      f      = F_REF * exp((eta - ETA_REF) * e)
      dfdeta = f * e
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

    real :: e

    if (f > F_SMALL) then
      call inv_fermi_dirac_integral_1h(f, eta, detadf)
    elseif (f > F_TINY) then
      eta    = - log(1 / f - GAMMA)
      detadf = 1.0 / (f * (1 - GAMMA * f))
    elseif (f > F_MICRO) then
      eta    = log(f)
      detadf = 1 / f
    else
      e      = (ETA_MICRO - ETA_REF) / (ETA_MICRO - log(F_REF))
      eta    = ETA_REF + e * log(f / F_REF)
      detadf = e / f
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
