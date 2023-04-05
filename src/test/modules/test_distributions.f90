m4_include(../../util/macro.f90.inc)

module test_distributions_m

  use distributions_m, only: bose_einstein, d_bose_einstein, fermi_dirac, d_fermi_dirac, maxwell_boltzmann, d_maxwell_boltzmann
  use distributions_m, only: fermi_dirac_integral_1h, inv_fermi_dirac_integral_1h, fermi_dirac_integral_m1h
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

    ! FD integral 1/2 and -1/2
    block
      real, parameter :: xx(16) = [-10.0, -3.0, -1.5, -1.0, 0.2, 1.5, 2.5, 3.5, 6.0, 8.0, 14.0, 17.0, 23.0, 36.0, 42.0, 53.0]
      real, parameter :: yy(16) = [ &
        4.5399201052641328E-05, 4.8933705696495779E-02, 2.0739818703202982E-01, 3.2779515926071155E-01, &
        8.9388099125744952E-01, 2.1448608775831140E+00, 3.6069753779503733E+00, 5.4580443887454599E+00, &
        1.1446599782570889E+01, 1.7355020120423224E+01, 3.9654596251512965E+01, 5.2953282656484951E+01, &
        8.3170418897886679E+01, 1.6264137964645799E+02, 2.0489979057766496E+02, 2.9038111053946135E+02 ]
      real, parameter :: zz(16) = [ &
        4.5398472360805495E-05, 4.8102635332204082E-02, 1.9330537025505338E-01, 2.9402761761145122E-01, &
        6.8317008347614267E-01, 1.2493233478527122E+00, 1.6663128834358317E+00, 2.0261936991004601E+00, &
        2.7271485514413181E+00, 3.1691247564959307E+00, 4.2129336531051786E+00, 4.6457008235030807E+00, &
        5.4072740945884877E+00, 6.7681194771811626E+00, 7.3110237937808748E+00, 8.2135198499452619E+00 ]

      integer :: i
      real    :: x(16), y(16), z(16), dydx(16), dxdy(16), dzdx(16)

      do i = 1, 16
        call fermi_dirac_integral_1h(xx(i), y(i), dydx(i))
        call fermi_dirac_integral_m1h(xx(i), z(i), dzdx(i))
        call inv_fermi_dirac_integral_1h(yy(i), x(i), dxdy(i))

        call tc%assert_eq(yy(i), y(i), 5e-16 * abs(yy(i)), "FD integral 1/2: "//int2str(i))
        call tc%assert_eq(zz(i), z(i), 5e-16 * abs(zz(i)), "FD integral -1/2: "//int2str(i))
        call tc%assert_eq(xx(i), x(i), 5e-15 * abs(xx(i)), "Inverse FD integral 1/2: "//int2str(i))
      end do
    end block

    call tc%finish()
  end subroutine

end module
