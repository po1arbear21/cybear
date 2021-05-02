module test_distributions_m

  use distributions_m
  use test_case_m, only: test_case

  implicit none

  private
  public test_distributions

contains

  subroutine test_distributions()
    type(test_case) :: tc

    call tc%init("distributions")

    ! BE
    call tc%assert_eq( 0.15651764274966565,   bose_einstein(2.0), 1e-15, "BE")
    call tc%assert_eq(-0.18101541524157762, d_bose_einstein(2.0), 1e-15, "BE")

    ! FD
    call tc%assert_eq( 0.5,                   fermi_dirac(0.0), epsilon(1.0), "FD")
    call tc%assert_eq( 0.11920292202211755,   fermi_dirac(2.0), 1e-15,        "FD")
    call tc%assert_eq(-0.10499358540350651, d_fermi_dirac(2.0), 1e-15,        "FD")

    ! MB
    call tc%assert_eq( 1.0,         maxwell_boltzmann(0.0), epsilon(1.0), "MB")
    call tc%assert_eq( exp(-2.0),   maxwell_boltzmann(2.0), epsilon(1.0), "MB")
    call tc%assert_eq(-1.0,       d_maxwell_boltzmann(0.0), epsilon(1.0), "MB")
    call tc%assert_eq(-exp(-2.0), d_maxwell_boltzmann(2.0), epsilon(1.0), "MB")

    ! usage of pointer
    block
      real, parameter               :: e = 1.234
      procedure(distr_feq), pointer :: feq

      feq => bose_einstein
      call tc%assert_eq(bose_einstein(e), feq(e), 0.0, "fptr: BE")

      feq => fermi_dirac
      call tc%assert_eq(fermi_dirac(e),   feq(e), 0.0, "fptr: FD")
    end block

    call tc%finish()
  end subroutine

end module
