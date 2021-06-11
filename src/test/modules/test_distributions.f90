module test_distributions_m

  use distributions_m, only: bose_einstein, d_bose_einstein, fermi_dirac, d_fermi_dirac, maxwell_boltzmann, d_maxwell_boltzmann
  use test_case_m,     only: test_case

  implicit none

  private
  public test_distributions

contains

  subroutine test_distributions()
    type(test_case) :: tc

    call tc%init("distributions")

    ! BE
    call tc%assert_eq( 0.99999950000008333e6,   bose_einstein(1e-6), 1e-9,         "BE")
    call tc%assert_eq( 0.15651764274966565,     bose_einstein(2.0),  epsilon(1.0), "BE")
    call tc%assert_eq(-0.18101541524157762,   d_bose_einstein(2.0),  epsilon(1.0), "BE")

    ! FD
    call tc%assert_eq( 0.5,                   fermi_dirac(0.0), epsilon(1.0), "FD")
    call tc%assert_eq( 0.11920292202211755,   fermi_dirac(2.0), epsilon(1.0), "FD")
    call tc%assert_eq(-0.10499358540350651, d_fermi_dirac(2.0), epsilon(1.0), "FD")

    ! MB
    call tc%assert_eq( 1.0,         maxwell_boltzmann(0.0), epsilon(1.0), "MB")
    call tc%assert_eq( exp(-2.0),   maxwell_boltzmann(2.0), epsilon(1.0), "MB")
    call tc%assert_eq(-1.0,       d_maxwell_boltzmann(0.0), epsilon(1.0), "MB")
    call tc%assert_eq(-exp(-2.0), d_maxwell_boltzmann(2.0), epsilon(1.0), "MB")

    call tc%finish()
  end subroutine

end module
