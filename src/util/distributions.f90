module distributions_m
  implicit none

  private
  public distr_feq
  public bose_einstein,     d_bose_einstein
  public fermi_dirac,       d_fermi_dirac
  public maxwell_boltzmann, d_maxwell_boltzmann

  abstract interface
    real function distr_feq(E)
      !! to easily switch between different distributions by procedure pointers
      real, intent(in) :: E
        !! energy
    end function
  end interface

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

    f = 1 / (exp(E) - 1)
  end function

  real function d_bose_einstein(E) result(dfdE)
    !! derivative of bose-einstein distribution
    real, intent(in) :: E
      !! energy

    dfdE = - exp(E) / (exp(E) - 1)**2
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

end module
