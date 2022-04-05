module distributions_m

  use, intrinsic :: ieee_arithmetic
  use math_m,     only: expm1, PI
#ifdef USE_QUADPACK
  use quadpack_m, only: quadpack_int
#endif

  implicit none

  private
  public bose_einstein,        d_bose_einstein
  public fermi_dirac,          d_fermi_dirac
  public maxwell_boltzmann,    d_maxwell_boltzmann
  public fermi_dirac_integral_approx

#ifdef USE_QUADPACK
  public fermi_dirac_integral
#endif

contains

  real function maxwell_boltzmann(E) result(f)
    !! maxwell-boltzmann distribution
    real, intent(in) :: E
      !! energy

    f = exp(-E)
  end function

  real function d_maxwell_boltzmann(E) result(dfdE)
    !! derivative of maxwell-boltzmann distribution
    real, intent(in) :: E
      !! energy

    dfdE = - exp(-E)
  end function

  real function bose_einstein(E) result(f)
    !! bose-einstein distribution
    real, intent(in) :: E
      !! energy

    f = 1 / expm1(E)
  end function

  real function d_bose_einstein(E) result(dfdE)
    !! derivative of bose-einstein distribution
    real, intent(in) :: E
      !! energy

    dfdE = - exp(E) / expm1(E)**2
  end function

  real function fermi_dirac(E) result(f)
    !! fermi-dirac distribution
    real, intent(in) :: E
      !! energy

    f = 1 / (exp(E) + 1)
  end function

  real function d_fermi_dirac(E) result(dfdE)
    !! derivative of fermi-dirac distribution
    real, intent(in) :: E
      !! energy

    dfdE = - exp(E) / (exp(E) + 1)**2
  end function

  subroutine fermi_dirac_integral_approx(x, y, dydx)
    !! Approximation of Fermi-Dirac integral F_/2
    !!  wiki: https://de.wikipedia.org/wiki/Fermi-Dirac-Integral
    real,           intent(in)  :: x
    real,           intent(out) :: y
    real, optional, intent(out) :: dydx
      !! derivative

    if (x < 1.3) then
      y = 1/(exp(-x) + 0.27)
      if (present(dydx)) dydx = exp(-x)/(exp(-x) + 0.27)**2
    else
      y = 4/(3*sqrt(PI)) * (x*x + PI*PI/6)**0.75
      if (present(dydx)) dydx =  2*x/sqrt(PI) / (x*x + PI*PI/6)**0.25
    end if
  end subroutine


#ifdef USE_QUADPACK
  recursive subroutine fermi_dirac_integral(x, y, j, atol, rtol, dydx)
    !! Fermi-Dirac integral
    !!  wiki: https://de.wikipedia.org/wiki/Fermi-Dirac-Integral
    !! y = y(x) = int_0^\infty \frac{t}{1+\exp(t-x)} \dd{t}

    use quadpack_m, only: quadpack_int

    real,           intent(in)  :: x
    real,           intent(out) :: y
    real, optional, intent(in)  :: j
      !! index j. default: 1/2
    real, optional, intent(in)  :: atol
      !! optional absolute tolerance. default: 1e-10
    real, optional, intent(in)  :: rtol
      !! optional relative tolerance. default: 1e-6
    real, optional, intent(out) :: dydx
      !! derivative

    real :: atol_, rtol_, j_

    j_ = 0.5
    if (present(j)) j_ = j
    atol_ = 1e-10
    if (present(atol)) atol_ = atol
    rtol_ = 1e-6
    if (present(rtol)) rtol_ = rtol

    call quadpack_int(f, 0.0, ieee_value(0.0, ieee_positive_inf), atol_, rtol_, y)
    y = y / gamma(j_+1)

    ! derivative: d(F_j(x))/dx = F_{j-1}(x)
    if (present(dydx)) call fermi_dirac_integral(x, dydx, j=j_-1, atol=atol, rtol=rtol)

  contains

    real function f(t) result(y)
      !! integrand: f(t) = f(t; x) = \frac{t^j}{1+\exp(t-x)}
      real, intent(in) :: t
        !! integration variable

      y = t**j_ / (1 + exp(t-x))
    end function

  end subroutine
#endif
end module


