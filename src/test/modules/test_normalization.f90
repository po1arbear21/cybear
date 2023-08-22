m4_include(../../util/macro.f90.inc)

module test_normalization_m

  use math_m,          only: PI
  use normalization_m, only: denorm, destruct_normconst, init_normconst, norm, get_norm_value
  use test_case_m,     only: test_case
  m4_ifdef({m4_intel},{use ifport})

  implicit none

  private
  public test_normalization

contains

  subroutine test_normalization()
    type(test_case) :: tc

    call tc%init("normalization")
    call test_norm_interface()
    call test_norm_compare()
    call tc%finish()
  contains

    subroutine test_norm_interface()
      real, parameter :: T_K = 300.0

      call init_normconst(T_K)

      ! check if norming and denorming is exchanged
      block
        real :: ener_eV, ener, ener_exp

        ener_eV  = 5.0                                ! 5eV
        ener_exp = 5.0 / (8.617333262145e-5 * 300.0)  ! = ener_eV / (kB_eV_per_K * T_K) = approx 2e3

        ener = norm(ener_eV, "eV")

        call tc%assert_eq(ener_exp, ener, 1e-13, 1e-15, "norming 5eV")
      end block

      ! check norm*denorm = identity
      block
        character(:), allocatable         :: unit
        character(*), parameter           :: unit_cm = "cm", unit_eV = "eV", unit_A_cm = "A/cm"
        integer                           :: i, j, n, i_unit
        real                              :: val_scal, val_scal_exp, nval_scal
        real, dimension(:),   allocatable :: val_arr, val_arr_exp, nval_arr
        real, dimension(:,:), allocatable :: val_mat, val_mat_exp, nval_mat

        n = 4
        allocate(val_arr(n),   val_arr_exp(n),   nval_arr(n))
        allocate(val_mat(n,n), val_mat_exp(n,n), nval_mat(n,n))

        allocate (character(0) :: unit)      ! remove gfortran warning

        do i_unit = 1, 3
          select case (i_unit)
            case (1)
              unit = unit_A_cm
            case (2)
              unit = unit_cm
            case (3)
              unit = unit_eV
          end select

          ! scalar
          val_scal_exp = rand()
          nval_scal    = norm(val_scal_exp, unit)
          val_scal     = denorm(nval_scal, unit)
          call tc%assert_eq(val_scal_exp, val_scal, 1e-14, 1e-16, "denorm(norm(scalar))")

          ! array
          do i = 1, n
            val_arr_exp(i) = rand()
          end do
          nval_arr = norm(val_arr_exp, unit)
          val_arr  = denorm(nval_arr, unit)
          call tc%assert_eq(val_arr_exp, val_arr, 1e-14, 1e-16, "denorm(norm(array))")

          ! matrix
          do i = 1, n
            do j = 1, n
              val_mat_exp(i,j) = rand()
            end do
          end do
          nval_mat = norm(val_mat_exp, unit)
          val_mat  = denorm(nval_mat, unit)
          call tc%assert_eq(val_mat_exp, val_mat, 1e-14, 1e-16, "denorm(norm(matrix))")
        end do
      end block

      ! check deg
      block
        real :: deg

        deg = norm(45.0, "deg")
        call tc%assert_eq(PI/4.0, deg, 1e-14, 1e-16, "norming degrees")
      end block

      call destruct_normconst()
    end subroutine

    subroutine test_norm_compare()
      real, parameter :: T_K = 300.0

      real, parameter :: EC = 1.602176634e-19
      real, parameter :: EM = 9.1093837015e-31
      real, parameter :: PLANCK = 6.62607015e-34 / (2*PI*EC)
      real, parameter :: BOLTZ = 1.380649e-23/EC
      real, parameter :: EPS0 = 8.8541878128e-12
      real, parameter :: MU0 = 1.25663706212e-6
      real, parameter :: TERA = 1e-12
      real, parameter :: GIGA = 1e-9
      real, parameter :: MEGA = 1e-6
      real, parameter :: KILO = 1e-3
      real, parameter :: CENTI = 1e2
      real, parameter :: MILLI = 1e3
      real, parameter :: MICRO = 1e6
      real, parameter :: NANO  = 1e9
      real, parameter :: PICO  = 1e12
      real, parameter :: FEMTO = 1e15

      ! local variables
      real :: ampere, coulomb, diel, farad, henry, hertz, kelvin, kilogram, meter, ohm, permby, second, volt, watt

      call init_normconst(T_K)

      ! basic units
      volt     = BOLTZ * T_K
      meter    = PLANCK / sqrt(EM / EC * volt)
      second   = PLANCK / volt
      hertz    = 1.0 / second
      kilogram = EM
      ampere   = EC / second
      coulomb  = ampere * second
      farad    = coulomb / volt
      ohm      = volt / ampere
      henry    = ohm * second
      watt     = volt * ampere
      kelvin   = T_K
      diel     = sqrt(EM * EC / volt) / (PLANCK * EPS0)
      permby   = henry / meter / MU0

      call tc%assert_eq(1.0, get_norm_value("((((1))))"), 0.0, 0.0, "test parentheses")

      call tc%assert_eq(1.0, get_norm_value("1"), 1e-14, 1e-16, "compare to old normalization: 1")
      call tc%assert_eq(180 / PI, get_norm_value("deg"), 1e-14, 1e-16, "compare to old normalization: deg")

      call tc%assert_eq(       meter,     get_norm_value("m"   ), 1e-14, 1e-16, "compare to old normalization: m"   )
      call tc%assert_eq(       meter**2,  get_norm_value("m^2" ), 1e-14, 1e-16, "compare to old normalization: m^2" )
      call tc%assert_eq(       meter**3 , get_norm_value("m^3" ), 1e-14, 1e-16, "compare to old normalization: m^3" )
      call tc%assert_eq( CENTI*meter    , get_norm_value("cm"  ), 1e-14, 1e-16, "compare to old normalization: cm"  )
      call tc%assert_eq((CENTI*meter)**2, get_norm_value("cm^2"), 1e-14, 1e-16, "compare to old normalization: cm^2")
      call tc%assert_eq((CENTI*meter)**3, get_norm_value("cm^3"), 1e-14, 1e-16, "compare to old normalization: cm^3")
      call tc%assert_eq( MILLI*meter    , get_norm_value("mm"  ), 1e-14, 1e-16, "compare to old normalization: mm"  )
      call tc%assert_eq((MILLI*meter)**2, get_norm_value("mm^2"), 1e-14, 1e-16, "compare to old normalization: mm^2")
      call tc%assert_eq((MILLI*meter)**3, get_norm_value("mm^3"), 1e-14, 1e-16, "compare to old normalization: mm^3")
      call tc%assert_eq( MICRO*meter    , get_norm_value("um"  ), 1e-14, 1e-16, "compare to old normalization: um"  )
      call tc%assert_eq((MICRO*meter)**2, get_norm_value("um^2"), 1e-14, 1e-16, "compare to old normalization: um^2")
      call tc%assert_eq((MICRO*meter)**3, get_norm_value("um^3"), 1e-14, 1e-16, "compare to old normalization: um^3")
      call tc%assert_eq( NANO *meter    , get_norm_value("nm"  ), 1e-14, 1e-16, "compare to old normalization: nm"  )
      call tc%assert_eq((NANO *meter)**2, get_norm_value("nm^2"), 1e-14, 1e-16, "compare to old normalization: nm^2")
      call tc%assert_eq((NANO *meter)**3, get_norm_value("nm^3"), 1e-14, 1e-16, "compare to old normalization: nm^3")
      call tc%assert_eq( PICO *meter    , get_norm_value("pm"  ), 1e-14, 1e-16, "compare to old normalization: pm"  )
      call tc%assert_eq((PICO *meter)**2, get_norm_value("pm^2"), 1e-14, 1e-16, "compare to old normalization: pm^2")
      call tc%assert_eq((PICO *meter)**3, get_norm_value("pm^3"), 1e-14, 1e-16, "compare to old normalization: pm^3")
      call tc%assert_eq( FEMTO*meter    , get_norm_value("fm"  ), 1e-14, 1e-16, "compare to old normalization: fm"  )
      call tc%assert_eq((FEMTO*meter)**2, get_norm_value("fm^2"), 1e-14, 1e-16, "compare to old normalization: fm^2")
      call tc%assert_eq((FEMTO*meter)**3, get_norm_value("fm^3"), 1e-14, 1e-16, "compare to old normalization: fm^3")

      call tc%assert_eq(1 /        meter    , get_norm_value("1/m"   ), 1e-14, 1e-16, "compare to old normalization: 1/m"   )
      call tc%assert_eq(1 /        meter**2 , get_norm_value("1/m^2" ), 1e-14, 1e-16, "compare to old normalization: 1/m^2" )
      call tc%assert_eq(1 /        meter**3 , get_norm_value("1/m^3" ), 1e-14, 1e-16, "compare to old normalization: 1/m^3" )
      call tc%assert_eq(1 / (CENTI*meter)   , get_norm_value("1/cm"  ), 1e-14, 1e-16, "compare to old normalization: 1/cm"  )
      call tc%assert_eq(1 / (CENTI*meter)**2, get_norm_value("1/cm^2"), 1e-14, 1e-16, "compare to old normalization: 1/cm^2")
      call tc%assert_eq(1 / (CENTI*meter)**3, get_norm_value("1/cm^3"), 1e-14, 1e-16, "compare to old normalization: 1/cm^3")
      call tc%assert_eq(1 / (MILLI*meter)   , get_norm_value("1/mm"  ), 1e-14, 1e-16, "compare to old normalization: 1/mm"  )
      call tc%assert_eq(1 / (MILLI*meter)**2, get_norm_value("1/mm^2"), 1e-14, 1e-16, "compare to old normalization: 1/mm^2")
      call tc%assert_eq(1 / (MILLI*meter)**3, get_norm_value("1/mm^3"), 1e-14, 1e-16, "compare to old normalization: 1/mm^3")
      call tc%assert_eq(1 / (MICRO*meter)   , get_norm_value("1/um"  ), 1e-14, 1e-16, "compare to old normalization: 1/um"  )
      call tc%assert_eq(1 / (MICRO*meter)**2, get_norm_value("1/um^2"), 1e-14, 1e-16, "compare to old normalization: 1/um^2")
      call tc%assert_eq(1 / (MICRO*meter)**3, get_norm_value("1/um^3"), 1e-14, 1e-16, "compare to old normalization: 1/um^3")
      call tc%assert_eq(1 / (NANO *meter)   , get_norm_value("1/nm"  ), 1e-14, 1e-16, "compare to old normalization: 1/nm"  )
      call tc%assert_eq(1 / (NANO *meter)**2, get_norm_value("1/nm^2"), 1e-14, 1e-16, "compare to old normalization: 1/nm^2")
      call tc%assert_eq(1 / (NANO *meter)**3, get_norm_value("1/nm^3"), 1e-14, 1e-16, "compare to old normalization: 1/nm^3")
      call tc%assert_eq(1 / (PICO *meter)   , get_norm_value("1/pm"  ), 1e-14, 1e-16, "compare to old normalization: 1/pm"  )
      call tc%assert_eq(1 / (PICO *meter)**2, get_norm_value("1/pm^2"), 1e-14, 1e-16, "compare to old normalization: 1/pm^2")
      call tc%assert_eq(1 / (PICO *meter)**3, get_norm_value("1/pm^3"), 1e-14, 1e-16, "compare to old normalization: 1/pm^3")
      call tc%assert_eq(1 / (FEMTO*meter)   , get_norm_value("1/fm"  ), 1e-14, 1e-16, "compare to old normalization: 1/fm"  )
      call tc%assert_eq(1 / (FEMTO*meter)**2, get_norm_value("1/fm^2"), 1e-14, 1e-16, "compare to old normalization: 1/fm^2")
      call tc%assert_eq(1 / (FEMTO*meter)**3, get_norm_value("1/fm^3"), 1e-14, 1e-16, "compare to old normalization: 1/fm^3")

      call tc%assert_eq(      second, get_norm_value("s" ), 1e-14, 1e-16, "compare to old normalization: s" )
      call tc%assert_eq(MILLI*second, get_norm_value("ms"), 1e-14, 1e-16, "compare to old normalization: ms")
      call tc%assert_eq(MICRO*second, get_norm_value("us"), 1e-14, 1e-16, "compare to old normalization: us")
      call tc%assert_eq(NANO *second, get_norm_value("ns"), 1e-14, 1e-16, "compare to old normalization: ns")
      call tc%assert_eq(PICO *second, get_norm_value("ps"), 1e-14, 1e-16, "compare to old normalization: ps")
      call tc%assert_eq(FEMTO*second, get_norm_value("fs"), 1e-14, 1e-16, "compare to old normalization: fs")

      call tc%assert_eq( 1 /        second , get_norm_value("1/s "), 1e-14, 1e-16, "compare to old normalization: 1/s")
      call tc%assert_eq( 1 / (MILLI*second), get_norm_value("1/ms"), 1e-14, 1e-16, "compare to old normalization: 1/ms")
      call tc%assert_eq( 1 / (MICRO*second), get_norm_value("1/us"), 1e-14, 1e-16, "compare to old normalization: 1/us")
      call tc%assert_eq( 1 / (NANO *second), get_norm_value("1/ns"), 1e-14, 1e-16, "compare to old normalization: 1/ns")
      call tc%assert_eq( 1 / (PICO *second), get_norm_value("1/ps"), 1e-14, 1e-16, "compare to old normalization: 1/ps")
      call tc%assert_eq( 1 / (FEMTO*second), get_norm_value("1/fs"), 1e-14, 1e-16, "compare to old normalization: 1/fs")

      call tc%assert_eq(      hertz,get_norm_value("Hz" ), 1e-14, 1e-16, "compare to old normalization: Hz" )
      call tc%assert_eq( KILO*hertz,get_norm_value("kHz"), 1e-14, 1e-16, "compare to old normalization: kHz")
      call tc%assert_eq( MEGA*hertz,get_norm_value("MHz"), 1e-14, 1e-16, "compare to old normalization: MHz")
      call tc%assert_eq( GIGA*hertz,get_norm_value("GHz"), 1e-14, 1e-16, "compare to old normalization: GHz")
      call tc%assert_eq( TERA*hertz,get_norm_value("THz"), 1e-14, 1e-16, "compare to old normalization: THz")

      call tc%assert_eq(       meter     / second   , get_norm_value("m/s"     ), 1e-14, 1e-16, "compare to old normalization: m/s"     )
      call tc%assert_eq((CENTI*meter)    / second   , get_norm_value("cm/s"    ), 1e-14, 1e-16, "compare to old normalization: cm/s"    )
      call tc%assert_eq(       meter**2  / second**2, get_norm_value("m^2/s^2" ), 1e-14, 1e-16, "compare to old normalization: m^2/s^2" )
      call tc%assert_eq((CENTI*meter)**2 / second**2, get_norm_value("cm^2/s^2"), 1e-14, 1e-16, "compare to old normalization: cm^2/s^2")

      call tc%assert_eq((CENTI*meter)**2 / second, get_norm_value("cm^2/s"), 1e-14, 1e-16, "compare to old normalization: cm^2/s")

      call tc%assert_eq(1 / (CENTI*meter)**2 / second, get_norm_value("1/cm^2/s"), 1e-14, 1e-16, "compare to old normalization: 1/cm^2/s")
      call tc%assert_eq(1 / (MICRO*meter)**2 / second, get_norm_value("1/um^2/s"), 1e-14, 1e-16, "compare to old normalization: 1/um^2/s")
      call tc%assert_eq(1 / (CENTI*meter)**3 / second, get_norm_value("1/cm^3/s"), 1e-14, 1e-16, "compare to old normalization: 1/cm^3/s")
      call tc%assert_eq(1 /        meter **3 / second, get_norm_value("1/m^3/s" ), 1e-14, 1e-16, "compare to old normalization: 1/m^3/s" )

      call tc%assert_eq( kilogram                         , get_norm_value("kg"    ), 1e-14, 1e-16, "compare to old normalization: kg"    )
      call tc%assert_eq( 1 / kilogram                     , get_norm_value("1/kg"  ), 1e-14, 1e-16, "compare to old normalization: 1/kg"  )
      call tc%assert_eq( kilogram       /        meter **3, get_norm_value("kg/m^3"), 1e-14, 1e-16, "compare to old normalization: kg/m^3")
      call tc%assert_eq((kilogram/KILO) / (CENTI*meter)**3, get_norm_value("g/cm^3"), 1e-14, 1e-16, "compare to old normalization: g/cm^3")

      call tc%assert_eq( KILO*volt, get_norm_value("kV"), 1e-14, 1e-16, "compare to old normalization: kV")
      call tc%assert_eq(      volt, get_norm_value("V" ), 1e-14, 1e-16, "compare to old normalization: V" )
      call tc%assert_eq(MILLI*volt, get_norm_value("mV"), 1e-14, 1e-16, "compare to old normalization: mV")
      call tc%assert_eq(MICRO*volt, get_norm_value("uV"), 1e-14, 1e-16, "compare to old normalization: uV")
      call tc%assert_eq( NANO*volt, get_norm_value("nV"), 1e-14, 1e-16, "compare to old normalization: nV")

      call tc%assert_eq(      volt  /        meter , get_norm_value("V/m"  ), 1e-14, 1e-16, "compare to old normalization: V/m"  )
      call tc%assert_eq(      volt  / (CENTI*meter), get_norm_value("V/cm" ), 1e-14, 1e-16, "compare to old normalization: V/cm" )
      call tc%assert_eq((KILO*volt) / (CENTI*meter), get_norm_value("kV/cm"), 1e-14, 1e-16, "compare to old normalization: kV/cm")
      call tc%assert_eq(1 /   volt                 , get_norm_value("1/V"  ), 1e-14, 1e-16, "compare to old normalization: 1/V"  )

      call tc%assert_eq(      volt                              , get_norm_value("eV"       ), 1e-14, 1e-16, "compare to old normalization: eV"       )
      call tc%assert_eq(MILLI*volt                              , get_norm_value("meV"      ), 1e-14, 1e-16, "compare to old normalization: meV"      )
      call tc%assert_eq(      volt /  (CENTI*meter)             , get_norm_value("eV/cm"    ), 1e-14, 1e-16, "compare to old normalization: eV/cm"    )
      call tc%assert_eq(      volt /  (CENTI*meter)**3          , get_norm_value("eV/cm^3"  ), 1e-14, 1e-16, "compare to old normalization: eV/cm^3"  )
      call tc%assert_eq(      volt / ((CENTI*meter)**2 * second), get_norm_value("eV/cm^2/s"), 1e-14, 1e-16, "compare to old normalization: eV/cm^2/s")
      call tc%assert_eq(1 /   volt                              , get_norm_value("1/eV"     ), 1e-14, 1e-16, "compare to old normalization: 1/eV"     )
      call tc%assert_eq(1 / (CENTI*meter)**3 / volt             , get_norm_value("1/cm^3/eV"), 1e-14, 1e-16, "compare to old normalization: 1/cm^3/eV")

      call tc%assert_eq((CENTI*meter)**2 / volt / second, get_norm_value("cm^2/V/s"), 1e-14, 1e-16, "compare to old normalization: cm^2/V/s")
      call tc%assert_eq(                   volt / second, get_norm_value("V/s"     ), 1e-14, 1e-16, "compare to old normalization: V/s"     )
      call tc%assert_eq(                 second / volt  , get_norm_value("s/V"     ), 1e-14, 1e-16, "compare to old normalization: s/V"     )

      call tc%assert_eq((volt/(CENTI*meter))**(-2.0/3.0) * kelvin * (CENTI*meter) / second, get_norm_value("(V/cm)^(-2/3)*K*cm/s"), 1e-14, 1e-16, "CURSED: (V/cm)^(-2/3)*K*cm/s")

      ! Adjusted "As" -> "A*s"
      call tc%assert_eq(       ampere                    , get_norm_value("A"      ), 1e-14, 1e-16, "compare to old normalization: A"      )
      call tc%assert_eq(       ampere  * second          , get_norm_value("A*s"    ), 1e-14, 1e-16, "compare to old normalization: As"     )
      call tc%assert_eq(       ampere  / meter           , get_norm_value("A/m"    ), 1e-14, 1e-16, "compare to old normalization: A/m"    )
      call tc%assert_eq(       ampere  / (CENTI*meter)   , get_norm_value("A/cm"   ), 1e-14, 1e-16, "compare to old normalization: A/cm"   )
      call tc%assert_eq(       ampere  / (MILLI*meter)   , get_norm_value("A/mm"   ), 1e-14, 1e-16, "compare to old normalization: A/mm"   )
      call tc%assert_eq(       ampere  / (MICRO*meter)   , get_norm_value("A/um"   ), 1e-14, 1e-16, "compare to old normalization: A/um"   )
      call tc%assert_eq(       ampere  / (NANO *meter)   , get_norm_value("A/nm"   ), 1e-14, 1e-16, "compare to old normalization: A/nm"   )
      call tc%assert_eq(       ampere  / (PICO *meter)   , get_norm_value("A/pm"   ), 1e-14, 1e-16, "compare to old normalization: A/pm"   )
      call tc%assert_eq(       ampere  / (FEMTO*meter)   , get_norm_value("A/fm"   ), 1e-14, 1e-16, "compare to old normalization: A/fm"   )
      call tc%assert_eq(       ampere  / (CENTI*meter)**2, get_norm_value("A/cm^2" ), 1e-14, 1e-16, "compare to old normalization: A/cm^2" )
      call tc%assert_eq(       ampere  / (MICRO*meter)**2, get_norm_value("A/um^2" ), 1e-14, 1e-16, "compare to old normalization: A/um^2" )
      call tc%assert_eq((MILLI*ampere) / (MICRO*meter)**2, get_norm_value("mA/um^2"), 1e-14, 1e-16, "compare to old normalization: mA/um^2")

      call tc%assert_eq(ampere / volt / (      meter), get_norm_value("A/V/m"  ), 1e-14, 1e-16, "compare to old normalization: A/V/m" )
      call tc%assert_eq(ampere / volt / (CENTI*meter), get_norm_value("A/V/cm" ), 1e-14, 1e-16, "compare to old normalization: A/V/cm")
      call tc%assert_eq(ampere / volt / (MILLI*meter), get_norm_value("A/V/mm" ), 1e-14, 1e-16, "compare to old normalization: A/V/mm")
      call tc%assert_eq(ampere / volt / (MICRO*meter), get_norm_value("A/V/um" ), 1e-14, 1e-16, "compare to old normalization: A/V/um")
      call tc%assert_eq(ampere / volt / (NANO *meter), get_norm_value("A/V/nm" ), 1e-14, 1e-16, "compare to old normalization: A/V/nm")
      call tc%assert_eq(ampere / volt / (PICO *meter), get_norm_value("A/V/pm" ), 1e-14, 1e-16, "compare to old normalization: A/V/pm")
      call tc%assert_eq(ampere / volt / (FEMTO*meter), get_norm_value("A/V/fm" ), 1e-14, 1e-16, "compare to old normalization: A/V/fm")

      call tc%assert_eq(ampere / volt / (      meter)**2, get_norm_value("A/V/m^2" ), 1e-14, 1e-16, "compare to old normalization: A/V/m^2" )
      call tc%assert_eq(ampere / volt / (CENTI*meter)**2, get_norm_value("A/V/cm^2"), 1e-14, 1e-16, "compare to old normalization: A/V/cm^2")
      call tc%assert_eq(ampere / volt / (MILLI*meter)**2, get_norm_value("A/V/mm^2"), 1e-14, 1e-16, "compare to old normalization: A/V/mm^2")
      call tc%assert_eq(ampere / volt / (MICRO*meter)**2, get_norm_value("A/V/um^2"), 1e-14, 1e-16, "compare to old normalization: A/V/um^2")
      call tc%assert_eq(ampere / volt / (NANO *meter)**2, get_norm_value("A/V/nm^2"), 1e-14, 1e-16, "compare to old normalization: A/V/nm^2")
      call tc%assert_eq(ampere / volt / (PICO *meter)**2, get_norm_value("A/V/pm^2"), 1e-14, 1e-16, "compare to old normalization: A/V/pm^2")
      call tc%assert_eq(ampere / volt / (FEMTO*meter)**2, get_norm_value("A/V/fm^2"), 1e-14, 1e-16, "compare to old normalization: A/V/fm^2")

      call tc%assert_eq(coulomb                   , get_norm_value("C"     ), 1e-14, 1e-16, "compare to old normalization: C"     )
      call tc%assert_eq(coulomb /        meter    , get_norm_value("C/m"   ), 1e-14, 1e-16, "compare to old normalization: C/m"   )
      call tc%assert_eq(coulomb / (CENTI*meter)   , get_norm_value("C/cm"  ), 1e-14, 1e-16, "compare to old normalization: C/cm"  )
      call tc%assert_eq(coulomb / (MILLI*meter)   , get_norm_value("C/mm"  ), 1e-14, 1e-16, "compare to old normalization: C/mm"  )
      call tc%assert_eq(coulomb / (MICRO*meter)   , get_norm_value("C/um"  ), 1e-14, 1e-16, "compare to old normalization: C/um"  )
      call tc%assert_eq(coulomb / (NANO *meter)   , get_norm_value("C/nm"  ), 1e-14, 1e-16, "compare to old normalization: C/nm"  )
      call tc%assert_eq(coulomb / (PICO *meter)   , get_norm_value("C/pm"  ), 1e-14, 1e-16, "compare to old normalization: C/pm"  )
      call tc%assert_eq(coulomb / (FEMTO*meter)   , get_norm_value("C/fm"  ), 1e-14, 1e-16, "compare to old normalization: C/fm"  )
      call tc%assert_eq(coulomb / (      meter)**2, get_norm_value("C/m^2" ), 1e-14, 1e-16, "compare to old normalization: C/m^2" )
      call tc%assert_eq(coulomb / (CENTI*meter)**2, get_norm_value("C/cm^2"), 1e-14, 1e-16, "compare to old normalization: C/cm^2")
      call tc%assert_eq(coulomb / (MILLI*meter)**2, get_norm_value("C/mm^2"), 1e-14, 1e-16, "compare to old normalization: C/mm^2")
      call tc%assert_eq(coulomb / (MICRO*meter)**2, get_norm_value("C/um^2"), 1e-14, 1e-16, "compare to old normalization: C/um^2")
      call tc%assert_eq(coulomb / (NANO *meter)**2, get_norm_value("C/nm^2"), 1e-14, 1e-16, "compare to old normalization: C/nm^2")
      call tc%assert_eq(coulomb / (PICO *meter)**2, get_norm_value("C/pm^2"), 1e-14, 1e-16, "compare to old normalization: C/pm^2")
      call tc%assert_eq(coulomb / (FEMTO*meter)**2, get_norm_value("C/fm^2"), 1e-14, 1e-16, "compare to old normalization: C/fm^2")
      call tc%assert_eq(coulomb / (      meter)**3, get_norm_value("C/m^3" ), 1e-14, 1e-16, "compare to old normalization: C/m^3" )
      call tc%assert_eq(coulomb / (CENTI*meter)**3, get_norm_value("C/cm^3"), 1e-14, 1e-16, "compare to old normalization: C/cm^3")
      call tc%assert_eq(coulomb / (MILLI*meter)**3, get_norm_value("C/mm^3"), 1e-14, 1e-16, "compare to old normalization: C/mm^3")
      call tc%assert_eq(coulomb / (MICRO*meter)**3, get_norm_value("C/um^3"), 1e-14, 1e-16, "compare to old normalization: C/um^3")
      call tc%assert_eq(coulomb / (NANO *meter)**3, get_norm_value("C/nm^3"), 1e-14, 1e-16, "compare to old normalization: C/nm^3")
      call tc%assert_eq(coulomb / (PICO *meter)**3, get_norm_value("C/pm^3"), 1e-14, 1e-16, "compare to old normalization: C/pm^3")
      call tc%assert_eq(coulomb / (FEMTO*meter)**3, get_norm_value("C/fm^3"), 1e-14, 1e-16, "compare to old normalization: C/fm^3")

      call tc%assert_eq(1 / ohm, get_norm_value("A/V" ), 1e-14, 1e-16, "compare to old normalization: A/V")

      call tc%assert_eq(MEGA *ohm, get_norm_value("MOhm"), 1e-14, 1e-16, "compare to old normalization: MOhm")
      call tc%assert_eq(KILO *ohm, get_norm_value("kOhm"), 1e-14, 1e-16, "compare to old normalization: kOhm")
      call tc%assert_eq(      ohm, get_norm_value("V/A" ), 1e-14, 1e-16, "compare to old normalization: V/A")
      call tc%assert_eq(      ohm, get_norm_value("Ohm" ), 1e-14, 1e-16, "compare to old normalization: Ohm")
      call tc%assert_eq(MILLI*ohm, get_norm_value("mOhm"), 1e-14, 1e-16, "compare to old normalization: mOhm")
      call tc%assert_eq(MICRO*ohm, get_norm_value("uOhm"), 1e-14, 1e-16, "compare to old normalization: uOhm")
      call tc%assert_eq(NANO *ohm, get_norm_value("nOhm"), 1e-14, 1e-16, "compare to old normalization: nOhm")
      call tc%assert_eq(PICO *ohm, get_norm_value("pOhm"), 1e-14, 1e-16, "compare to old normalization: pOhm")
      call tc%assert_eq(FEMTO*ohm, get_norm_value("fOhm"), 1e-14, 1e-16, "compare to old normalization: fOhm")

      call tc%assert_eq(      farad, get_norm_value("F" ), 1e-14, 1e-16, "compare to old normalization: F")
      call tc%assert_eq(MILLI*farad, get_norm_value("mF"), 1e-14, 1e-16, "compare to old normalization: mF")
      call tc%assert_eq(MICRO*farad, get_norm_value("uF"), 1e-14, 1e-16, "compare to old normalization: uF")
      call tc%assert_eq(NANO *farad, get_norm_value("nF"), 1e-14, 1e-16, "compare to old normalization: nF")
      call tc%assert_eq(PICO *farad, get_norm_value("pF"), 1e-14, 1e-16, "compare to old normalization: pF")
      call tc%assert_eq(FEMTO*farad, get_norm_value("fF"), 1e-14, 1e-16, "compare to old normalization: fF")

      call tc%assert_eq(      farad / meter, get_norm_value("F/m" ), 1e-14, 1e-16, "compare to old normalization: F/m" )
      call tc%assert_eq(MILLI*farad / meter, get_norm_value("mF/m"), 1e-14, 1e-16, "compare to old normalization: mF/m")
      call tc%assert_eq(MICRO*farad / meter, get_norm_value("uF/m"), 1e-14, 1e-16, "compare to old normalization: uF/m")
      call tc%assert_eq(NANO *farad / meter, get_norm_value("nF/m"), 1e-14, 1e-16, "compare to old normalization: nF/m")
      call tc%assert_eq(PICO *farad / meter, get_norm_value("pF/m"), 1e-14, 1e-16, "compare to old normalization: pF/m")
      call tc%assert_eq(FEMTO*farad / meter, get_norm_value("fF/m"), 1e-14, 1e-16, "compare to old normalization: fF/m")

      call tc%assert_eq(      henry, get_norm_value("H" ), 1e-14, 1e-16, "compare to old normalization: H" )
      call tc%assert_eq(MILLI*henry, get_norm_value("mH"), 1e-14, 1e-16, "compare to old normalization: mH")
      call tc%assert_eq(MICRO*henry, get_norm_value("uH"), 1e-14, 1e-16, "compare to old normalization: uH")
      call tc%assert_eq(NANO *henry, get_norm_value("nH"), 1e-14, 1e-16, "compare to old normalization: nH")
      call tc%assert_eq(PICO *henry, get_norm_value("pH"), 1e-14, 1e-16, "compare to old normalization: pH")
      call tc%assert_eq(FEMTO*henry, get_norm_value("fH"), 1e-14, 1e-16, "compare to old normalization: fH")

      call tc%assert_eq(      henry / meter, get_norm_value("H/m" ), 1e-14, 1e-16, "compare to old normalization: H/m" )
      call tc%assert_eq(MILLI*henry / meter, get_norm_value("mH/m"), 1e-14, 1e-16, "compare to old normalization: mH/m")
      call tc%assert_eq(MICRO*henry / meter, get_norm_value("uH/m"), 1e-14, 1e-16, "compare to old normalization: uH/m")
      call tc%assert_eq(NANO *henry / meter, get_norm_value("nH/m"), 1e-14, 1e-16, "compare to old normalization: nH/m")
      call tc%assert_eq(PICO *henry / meter, get_norm_value("pH/m"), 1e-14, 1e-16, "compare to old normalization: pH/m")
      call tc%assert_eq(FEMTO*henry / meter, get_norm_value("fH/m"), 1e-14, 1e-16, "compare to old normalization: fH/m")

      call tc%assert_eq(      watt, get_norm_value("W" ), 1e-14, 1e-16, "compare to old normalization: W")
      call tc%assert_eq(MILLI*watt, get_norm_value("mW"), 1e-14, 1e-16, "compare to old normalization: mW")
      call tc%assert_eq(MICRO*watt, get_norm_value("uW"), 1e-14, 1e-16, "compare to old normalization: uW")
      call tc%assert_eq(NANO *watt, get_norm_value("nW"), 1e-14, 1e-16, "compare to old normalization: nW")
      call tc%assert_eq(PICO *watt, get_norm_value("pW"), 1e-14, 1e-16, "compare to old normalization: pW")
      call tc%assert_eq(FEMTO*watt, get_norm_value("fW"), 1e-14, 1e-16, "compare to old normalization: fW")

      call tc%assert_eq(watt / (        meter), get_norm_value("W/m" ), 1e-14, 1e-16, "compare to old normalization: W/m")
      call tc%assert_eq(watt / (CENTI * meter), get_norm_value("W/cm"), 1e-14, 1e-16, "compare to old normalization: W/cm")
      call tc%assert_eq(watt / (MILLI * meter), get_norm_value("W/mm"), 1e-14, 1e-16, "compare to old normalization: W/mm")
      call tc%assert_eq(watt / (MICRO * meter), get_norm_value("W/um"), 1e-14, 1e-16, "compare to old normalization: W/um")
      call tc%assert_eq(watt / (NANO  * meter), get_norm_value("W/nm"), 1e-14, 1e-16, "compare to old normalization: W/nm")
      call tc%assert_eq(watt / (PICO  * meter), get_norm_value("W/pm"), 1e-14, 1e-16, "compare to old normalization: W/pm")
      call tc%assert_eq(watt / (FEMTO * meter), get_norm_value("W/fm"), 1e-14, 1e-16, "compare to old normalization: W/fm")

      call tc%assert_eq(ampere / watt, get_norm_value("A/W"), 1e-14, 1e-16, "compare to old normalization: A/W")

      call tc%assert_eq(diel, get_norm_value("eps0"), 1e-14, 1e-16, "compare to old normalization: eps0")

      call tc%assert_eq(permby, get_norm_value("mu0"), 1e-14, 1e-16, "compare to old normalization: mu0")

      call tc%assert_eq(kelvin, get_norm_value("K"), 1e-14, 1e-16, "compare to old normalization: K")

      call destruct_normconst()
    end subroutine
  end subroutine

end module
