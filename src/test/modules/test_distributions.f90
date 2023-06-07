m4_include(../../util/macro.f90.inc)

module test_distributions_m

  use distributions_m, only: bose_einstein, d_bose_einstein, fermi_dirac, d_fermi_dirac, maxwell_boltzmann, d_maxwell_boltzmann
  use distributions_m, only: fermi_dirac_integral_m1h, fermi_dirac_integral_1h, inv_fermi_dirac_integral_1h, fermi_dirac_integral_3h
  use test_case_m,     only: test_case
  use util_m,          only: int2str

  implicit none

  private
  public test_distributions

contains

  subroutine test_distributions()
    type(test_case) :: tc

    call tc%init("distributions")

    ! MB
    call tc%assert_eq( 1.0,         maxwell_boltzmann(0.0), epsilon(1.0), "MB")
    call tc%assert_eq( exp(-2.0),   maxwell_boltzmann(2.0), epsilon(1.0), "MB")
    call tc%assert_eq(-1.0,       d_maxwell_boltzmann(0.0), epsilon(1.0), "MB")
    call tc%assert_eq(-exp(-2.0), d_maxwell_boltzmann(2.0), epsilon(1.0), "MB")

    ! BE
    call tc%assert_eq( 0.99999950000008333e6,   bose_einstein(1e-6), 1e-9,         "BE")
    call tc%assert_eq( 0.15651764274966565,     bose_einstein(2.0),  epsilon(1.0), "BE")
    call tc%assert_eq(-0.18101541524157762,   d_bose_einstein(2.0),  epsilon(1.0), "BE")

    ! FD
    call tc%assert_eq( 0.5,                   fermi_dirac(0.0), epsilon(1.0), "FD")
    call tc%assert_eq( 0.11920292202211755,   fermi_dirac(2.0), epsilon(1.0), "FD")
    call tc%assert_eq(-0.10499358540350651, d_fermi_dirac(2.0), epsilon(1.0), "FD")

    ! FD integral
    m4_divert(m4_ifdef({m4_quadpack},0,-1))
    block
      use distributions_m, only: fermi_dirac_integral
      real :: res

      call fermi_dirac_integral(0.0, res)
      call tc%assert_eq(0.76514702462540794, res, 1e-8, "FD intg")
    end block
    m4_divert(0)

    ! FD integral -1/2, 1/2 and 3/2
    block
      real, parameter :: x(16) = [-10.0, -3.0, -1.5, -1.0, 0.2, 1.5, 2.5, 3.5, 6.0, 8.0, 14.0, 17.0, 23.0, 36.0, 42.0, 53.0]
      real, parameter :: Fm1h_exp(16) = [ &
        4.5398472360805495E-05, 4.8102635332204082E-02, 1.9330537025505338E-01, 2.9402761761145122E-01, &
        6.8317008347614267E-01, 1.2493233478527122E+00, 1.6663128834358317E+00, 2.0261936991004601E+00, &
        2.7271485514413181E+00, 3.1691247564959307E+00, 4.2129336531051786E+00, 4.6457008235030807E+00, &
        5.4072740945884877E+00, 6.7681194771811626E+00, 7.3110237937808748E+00, 8.2135198499452619E+00 ]
      real, parameter :: F1h_exp(16) = [ &
        4.5399201052641328E-05, 4.8933705696495779E-02, 2.0739818703202982E-01, 3.2779515926071155E-01, &
        8.9388099125744952E-01, 2.1448608775831140E+00, 3.6069753779503733E+00, 5.4580443887454599E+00, &
        1.1446599782570889E+01, 1.7355020120423224E+01, 3.9654596251512965E+01, 5.2953282656484951E+01, &
        8.3170418897886679E+01, 1.6264137964645799E+02, 2.0489979057766496E+02, 2.9038111053946135E+02 ]
      real, parameter :: F3h_exp(16) = [ &
        4.5399565404561763E-05, 4.9356612790684162E-02, 2.1497282584741158E-01, 3.4667479479905742E-01, &
        1.0328417428619416E+00, 2.9277494127932173E+00, 5.7688793122105166E+00, 1.0271411150617323E+01, &
        3.1038703861050414E+01, 5.9693209004975021E+01, 2.2760479848966514E+02, 3.6619220404853726E+02, &
        7.7228128412821934E+02, 2.3509412157009575E+03, 3.4519365010156854E+03, 6.1668874670831883E+03 ]

      integer :: i
      real    :: Fm1h(16), F1h(16), IF1h(16), F3h(16), dFm1h(16), dF1h(16), dIF1h(16), dF3h(16)

      do i = 1, 16
        call fermi_dirac_integral_m1h(x(i), Fm1h(i), dFm1h(i))
        call fermi_dirac_integral_1h(x(i), F1h(i), dF1h(i))
        call inv_fermi_dirac_integral_1h(F1h_exp(i), IF1h(i), dIF1h(i))
        call fermi_dirac_integral_3h(x(i), F3h(i), dF3h(i))

        call tc%assert_eq(Fm1h_exp(i), Fm1h(i), 5e-16 * abs(Fm1h_exp(i)), "FD integral -1/2: "     // int2str(i))
        call tc%assert_eq( F1h_exp(i),  F1h(i), 5e-16 * abs( F1h_exp(i)), "FD integral 1/2: "      // int2str(i))
        call tc%assert_eq(       x(i), IF1h(i), 5e-15 * abs(       x(i)), "Inv. FD integral 1/2: " // int2str(i))
        call tc%assert_eq( F3h_exp(i),  F3h(i), 5e-16 * abs( F3h_exp(i)), "FD integral 3/2: "      // int2str(i))
      end do
    end block

    call tc%finish()
  end subroutine

end module
