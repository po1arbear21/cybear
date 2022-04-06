#include "macro.f90.inc"

module normalization_m

  use error_m,  only: program_error
  use map_m,    only: map_string_real, mapnode_string_real
  use math_m,   only: PI
  use string_m, only: new_string

  implicit none

  private
  public init_normconst, destruct_normconst
  public norm, denorm
  public normalization

  type normalization
    type(map_string_real) :: unit_const
      !! phyiscal unit tokens (e.g. "eV" or "kV/cm") -> corresponding normalization constant
  contains
    procedure :: init     => normalization_init
    procedure :: destruct => normalization_destruct
  end type

  type(normalization) :: normconst
    !! global normalization object

  interface norm
    module procedure :: norm_scalar_r
    module procedure :: norm_array_r
    module procedure :: norm_matrix_r
    module procedure :: norm_scalar_c
    module procedure :: norm_array_c
    module procedure :: norm_matrix_c
  end interface

  interface denorm
    module procedure :: denorm_scalar_r
    module procedure :: denorm_array_r
    module procedure :: denorm_matrix_r
    module procedure :: denorm_scalar_c
    module procedure :: denorm_array_c
    module procedure :: denorm_matrix_c
  end interface

contains

  subroutine init_normconst(T)
    !! initialize global normalization object

    real, intent(in) :: T
      !! temperature in Kelvin

    call normconst%init(T)
  end subroutine

  subroutine destruct_normconst()
    !! destruct global normalization object

    call normconst%destruct()
  end subroutine

#define T r
#define TT real
#define NORM
#include "normalization_imp.f90.inc"

#define T c
#define TT complex
#define NORM
#include "normalization_imp.f90.inc"

#define T r
#define TT real
#define DENORM
#include "normalization_imp.f90.inc"

#define T c
#define TT complex
#define DENORM
#include "normalization_imp.f90.inc"

  function get_norm_value(unit, n) result(val)
    !! usage:
    !!  in norm_X: normed_value = denormed_value * val

    character(*),                  intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
    real                                      :: val

    type(mapnode_string_real), pointer :: node

    ! search for unit
    if (present(n)) then
      node => n%unit_const%find(new_string(unit))
    else
      node => normconst%unit_const%find(new_string(unit))
    end if
    if (.not. associated(node)) then
      print "(A)", "Unit: "//trim(unit)
      call program_error("Unit not found in normalization object!")
    end if

    ! return value
    val = node%value
  end function

  subroutine normalization_init(this, T)
    !! initialize normalization constants

    class(normalization), intent(out) :: this
    real,                 intent(in)  :: T
      !! temperature

    ! constants
    real, parameter :: EC     = 1.602176634e-19               ! elementary charge         [ As   ]
      !! exact
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?e
    real, parameter :: EM     = 9.1093837015e-31              ! electron rest mass        [ kg   ]
      !! rel err: 3e-10
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?me
    real, parameter :: PLANCK = 6.62607015e-34 / (2*PI*EC)    ! reduced Planck's constant [ eVs  ]
      !! exact
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?h
    real, parameter :: BOLTZ  = 1.380649e-23/EC               ! Boltzmann's constant      [ eV/K ]
      !! exact
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?k
    real, parameter :: EPS0   = 8.8541878128e-12              ! vacuum permittivity       [ F/m  ]
      !! rel err: 1.5e-10
      !! nist link: https://physics.nist.gov/cgi-bin/cuu/Value?ep0

    ! metric prefixes: INVERSE values for normalization
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
    real :: meter, second, kilogram, volt, ampere, kelvin, diel

    ! basic units
    volt     = BOLTZ * T
    meter    = PLANCK / sqrt(EM / EC * volt)
    second   = PLANCK / volt
    kilogram = EM
    ampere   = EC / second
    kelvin   = T
    diel     = sqrt(EM * EC / volt) / (PLANCK * EPS0)

    ! initialize map
    call this%unit_const%init()

    call this%unit_const%insert(new_string("1"), 1.0)
    call this%unit_const%insert(new_string("deg"), 180 / PI)

    call this%unit_const%insert(new_string("m"   ),        meter     )
    call this%unit_const%insert(new_string("m^2" ),        meter**2  )
    call this%unit_const%insert(new_string("m^3" ),        meter**3  )
    call this%unit_const%insert(new_string("cm"  ),  CENTI*meter     )
    call this%unit_const%insert(new_string("cm^2"), (CENTI*meter)**2 )
    call this%unit_const%insert(new_string("cm^3"), (CENTI*meter)**3 )
    call this%unit_const%insert(new_string("mm"  ),  MILLI*meter     )
    call this%unit_const%insert(new_string("mm^2"), (MILLI*meter)**2 )
    call this%unit_const%insert(new_string("mm^3"), (MILLI*meter)**3 )
    call this%unit_const%insert(new_string("um"  ),  MICRO*meter     )
    call this%unit_const%insert(new_string("um^2"), (MICRO*meter)**2 )
    call this%unit_const%insert(new_string("um^3"), (MICRO*meter)**3 )
    call this%unit_const%insert(new_string("nm"  ),  NANO *meter     )
    call this%unit_const%insert(new_string("nm^2"), (NANO *meter)**2 )
    call this%unit_const%insert(new_string("nm^3"), (NANO *meter)**3 )
    call this%unit_const%insert(new_string("pm"  ),  PICO *meter     )
    call this%unit_const%insert(new_string("pm^2"), (PICO *meter)**2 )
    call this%unit_const%insert(new_string("pm^3"), (PICO *meter)**3 )
    call this%unit_const%insert(new_string("fm"  ),  FEMTO*meter     )
    call this%unit_const%insert(new_string("fm^2"), (FEMTO*meter)**2 )
    call this%unit_const%insert(new_string("fm^3"), (FEMTO*meter)**3 )

    call this%unit_const%insert(new_string("1/m"   ), 1 /        meter     )
    call this%unit_const%insert(new_string("1/m^2" ), 1 /        meter**2  )
    call this%unit_const%insert(new_string("1/m^3" ), 1 /        meter**3  )
    call this%unit_const%insert(new_string("1/cm"  ), 1 / (CENTI*meter)     )
    call this%unit_const%insert(new_string("1/cm^2"), 1 / (CENTI*meter)**2 )
    call this%unit_const%insert(new_string("1/cm^3"), 1 / (CENTI*meter)**3 )
    call this%unit_const%insert(new_string("1/mm"  ), 1 / (MILLI*meter)     )
    call this%unit_const%insert(new_string("1/mm^2"), 1 / (MILLI*meter)**2 )
    call this%unit_const%insert(new_string("1/mm^3"), 1 / (MILLI*meter)**3 )
    call this%unit_const%insert(new_string("1/um"  ), 1 / (MICRO*meter)     )
    call this%unit_const%insert(new_string("1/um^2"), 1 / (MICRO*meter)**2 )
    call this%unit_const%insert(new_string("1/um^3"), 1 / (MICRO*meter)**3 )
    call this%unit_const%insert(new_string("1/nm"  ), 1 / (NANO *meter)     )
    call this%unit_const%insert(new_string("1/nm^2"), 1 / (NANO *meter)**2 )
    call this%unit_const%insert(new_string("1/nm^3"), 1 / (NANO *meter)**3 )
    call this%unit_const%insert(new_string("1/pm"  ), 1 / (PICO *meter)     )
    call this%unit_const%insert(new_string("1/pm^2"), 1 / (PICO *meter)**2 )
    call this%unit_const%insert(new_string("1/pm^3"), 1 / (PICO *meter)**3 )
    call this%unit_const%insert(new_string("1/fm"  ), 1 / (FEMTO*meter)     )
    call this%unit_const%insert(new_string("1/fm^2"), 1 / (FEMTO*meter)**2 )
    call this%unit_const%insert(new_string("1/fm^3"), 1 / (FEMTO*meter)**3 )

    call this%unit_const%insert(new_string("s" ),       second )
    call this%unit_const%insert(new_string("ms"), MILLI*second )
    call this%unit_const%insert(new_string("us"), MICRO*second )
    call this%unit_const%insert(new_string("ns"), NANO *second )
    call this%unit_const%insert(new_string("ps"), PICO *second )
    call this%unit_const%insert(new_string("fs"), FEMTO*second )

    call this%unit_const%insert(new_string("1/s" ), 1 /        second  )
    call this%unit_const%insert(new_string("1/ms"), 1 / (MILLI*second) )
    call this%unit_const%insert(new_string("1/us"), 1 / (MICRO*second) )
    call this%unit_const%insert(new_string("1/ns"), 1 / (NANO *second) )
    call this%unit_const%insert(new_string("1/ps"), 1 / (PICO *second) )
    call this%unit_const%insert(new_string("1/fs"), 1 / (FEMTO*second) )

    associate (hertz => 1 / second)
      call this%unit_const%insert(new_string("Hz" ),      hertz )
      call this%unit_const%insert(new_string("kHz"), KILO*hertz )
      call this%unit_const%insert(new_string("MHz"), MEGA*hertz )
      call this%unit_const%insert(new_string("GHz"), GIGA*hertz )
      call this%unit_const%insert(new_string("THz"), TERA*hertz )
    end associate

    call this%unit_const%insert(new_string("m/s"    ),         meter     / second    )
    call this%unit_const%insert(new_string("cm/s"   ),  (CENTI*meter)    / second    )
    call this%unit_const%insert(new_string("m^2/s^2"),         meter**2  / second**2 )
    call this%unit_const%insert(new_string("cm^2/s^2"), (CENTI*meter)**2 / second**2 )

    call this%unit_const%insert(new_string("cm^2/s"), (CENTI*meter)**2 / second )

    call this%unit_const%insert(new_string("1/cm^2/s"), 1 / (CENTI*meter)**2 / second )
    call this%unit_const%insert(new_string("1/um^2/s"), 1 / (MICRO*meter)**2 / second )
    call this%unit_const%insert(new_string("1/m^3/s" ), 1 /        meter**3  / second )

    call this%unit_const%insert(new_string("kg"    ),  kilogram                          )
    call this%unit_const%insert(new_string("kg/m^3"),  kilogram       /        meter **3 )
    call this%unit_const%insert(new_string("g/cm^3"), (kilogram/KILO) / (CENTI*meter)**3 )

    call this%unit_const%insert(new_string("kV"   ),  KILO*volt)
    call this%unit_const%insert(new_string("V"    ),       volt)
    call this%unit_const%insert(new_string("mV"   ), MILLI*volt)
    call this%unit_const%insert(new_string("uV"   ), MICRO*volt)
    call this%unit_const%insert(new_string("nV"   ),  NANO*volt)

    call this%unit_const%insert(new_string("V/m"  ),       volt  /        meter  )
    call this%unit_const%insert(new_string("V/cm" ),       volt  / (CENTI*meter) )
    call this%unit_const%insert(new_string("kV/cm"), (KILO*volt) / (CENTI*meter) )
    call this%unit_const%insert(new_string("1/V"  ), 1 /   volt                  )

    call this%unit_const%insert(new_string("eV"       ),       volt                               )
    call this%unit_const%insert(new_string("meV"      ), MILLI*volt                               )
    call this%unit_const%insert(new_string("eV/cm"    ),       volt /  (CENTI*meter)              )
    call this%unit_const%insert(new_string("eV/cm^3"  ),       volt /  (CENTI*meter)**3           )
    call this%unit_const%insert(new_string("eV/cm^2/s"),       volt / ((CENTI*meter)**2 * second) )
    call this%unit_const%insert(new_string("1/eV"     ), 1 /   volt                               )
    call this%unit_const%insert(new_string("1/cm^3/eV"), 1 / (CENTI*meter)**3 / volt              )

    call this%unit_const%insert(new_string("cm^2/V/s"            ), (CENTI*meter)**2 / volt / second )
    call this%unit_const%insert(new_string("V/s"                 ),                    volt / second )
    call this%unit_const%insert(new_string("s/V"                 ),                  second / volt   )
    call this%unit_const%insert(new_string("(V/cm)^(-2/3)*K*cm/s"), &
      & (volt/(CENTI*meter))**(-2.0/3.0) * kelvin * (CENTI*meter) / second )

    call this%unit_const%insert(new_string("A"      ),        ampere                         )
    call this%unit_const%insert(new_string("As"     ),        ampere  * second               )
    call this%unit_const%insert(new_string("A/m"    ),        ampere  / meter                )
    call this%unit_const%insert(new_string("A/cm"   ),        ampere  / (CENTI*meter)        )
    call this%unit_const%insert(new_string("A/mm"   ),        ampere  / (MILLI*meter)        )
    call this%unit_const%insert(new_string("A/um"   ),        ampere  / (MICRO*meter)        )
    call this%unit_const%insert(new_string("A/nm"   ),        ampere  / (NANO *meter)        )
    call this%unit_const%insert(new_string("A/pm"   ),        ampere  / (PICO *meter)        )
    call this%unit_const%insert(new_string("A/fm"   ),        ampere  / (FEMTO*meter)        )
    call this%unit_const%insert(new_string("A/cm^2" ),        ampere  / (CENTI*meter)**2     )
    call this%unit_const%insert(new_string("A/um^2" ),        ampere  / (MICRO*meter)**2     )
    call this%unit_const%insert(new_string("mA/um^2"), (MILLI*ampere) / (MICRO*meter)**2     )

    call this%unit_const%insert(new_string("A/V/m"  ), ampere / volt / (      meter))
    call this%unit_const%insert(new_string("A/V/cm" ), ampere / volt / (CENTI*meter))
    call this%unit_const%insert(new_string("A/V/mm" ), ampere / volt / (MILLI*meter))
    call this%unit_const%insert(new_string("A/V/um" ), ampere / volt / (MICRO*meter))
    call this%unit_const%insert(new_string("A/V/nm" ), ampere / volt / (NANO *meter))
    call this%unit_const%insert(new_string("A/V/pm" ), ampere / volt / (PICO *meter))
    call this%unit_const%insert(new_string("A/V/fm" ), ampere / volt / (FEMTO*meter))

    call this%unit_const%insert(new_string("A/V/m^2" ), ampere / volt / (      meter)**2 )
    call this%unit_const%insert(new_string("A/V/cm^2"), ampere / volt / (CENTI*meter)**2 )
    call this%unit_const%insert(new_string("A/V/mm^2"), ampere / volt / (MILLI*meter)**2 )
    call this%unit_const%insert(new_string("A/V/um^2"), ampere / volt / (MICRO*meter)**2 )
    call this%unit_const%insert(new_string("A/V/nm^2"), ampere / volt / (NANO *meter)**2 )
    call this%unit_const%insert(new_string("A/V/pm^2"), ampere / volt / (PICO *meter)**2 )
    call this%unit_const%insert(new_string("A/V/fm^2"), ampere / volt / (FEMTO*meter)**2 )

    associate (coulomb => ampere * second)
      call this%unit_const%insert(new_string("C"     ), coulomb                    )
      call this%unit_const%insert(new_string("C/m"   ), coulomb /        meter     )
      call this%unit_const%insert(new_string("C/cm"  ), coulomb / (CENTI*meter)    )
      call this%unit_const%insert(new_string("C/mm"  ), coulomb / (MILLI*meter)    )
      call this%unit_const%insert(new_string("C/um"  ), coulomb / (MICRO*meter)    )
      call this%unit_const%insert(new_string("C/nm"  ), coulomb / (NANO *meter)    )
      call this%unit_const%insert(new_string("C/pm"  ), coulomb / (PICO *meter)    )
      call this%unit_const%insert(new_string("C/fm"  ), coulomb / (FEMTO*meter)    )
      call this%unit_const%insert(new_string("C/m^2" ), coulomb / (      meter)**2 )
      call this%unit_const%insert(new_string("C/cm^2"), coulomb / (CENTI*meter)**2 )
      call this%unit_const%insert(new_string("C/mm^2"), coulomb / (MILLI*meter)**2 )
      call this%unit_const%insert(new_string("C/um^2"), coulomb / (MICRO*meter)**2 )
      call this%unit_const%insert(new_string("C/nm^2"), coulomb / (NANO *meter)**2 )
      call this%unit_const%insert(new_string("C/pm^2"), coulomb / (PICO *meter)**2 )
      call this%unit_const%insert(new_string("C/fm^2"), coulomb / (FEMTO*meter)**2 )
      call this%unit_const%insert(new_string("C/m^3" ), coulomb / (      meter)**3 )
      call this%unit_const%insert(new_string("C/cm^3"), coulomb / (CENTI*meter)**3 )
      call this%unit_const%insert(new_string("C/mm^3"), coulomb / (MILLI*meter)**3 )
      call this%unit_const%insert(new_string("C/um^3"), coulomb / (MICRO*meter)**3 )
      call this%unit_const%insert(new_string("C/nm^3"), coulomb / (NANO *meter)**3 )
      call this%unit_const%insert(new_string("C/pm^3"), coulomb / (PICO *meter)**3 )
      call this%unit_const%insert(new_string("C/fm^3"), coulomb / (FEMTO*meter)**3 )
    end associate

    associate (ohm => volt / ampere)
      call this%unit_const%insert(new_string("A/V" ), 1 / ohm )

      call this%unit_const%insert(new_string("MOhm"), MEGA *ohm )
      call this%unit_const%insert(new_string("kOhm"), KILO *ohm )
      call this%unit_const%insert(new_string("V/A" ),       ohm )
      call this%unit_const%insert(new_string("Ohm" ),       ohm )
      call this%unit_const%insert(new_string("mOhm"), MILLI*ohm )
      call this%unit_const%insert(new_string("uOhm"), MICRO*ohm )
      call this%unit_const%insert(new_string("nOhm"), NANO *ohm )
      call this%unit_const%insert(new_string("pOhm"), PICO *ohm )
      call this%unit_const%insert(new_string("fOhm"), FEMTO*ohm )
    end associate

    associate (farad => ampere * second / volt)
      call this%unit_const%insert(new_string("F" ),       farad  )
      call this%unit_const%insert(new_string("mF"), MILLI*farad )
      call this%unit_const%insert(new_string("uF"), MICRO*farad )
      call this%unit_const%insert(new_string("nF"), NANO *farad )
      call this%unit_const%insert(new_string("pF"), PICO *farad )
      call this%unit_const%insert(new_string("fF"), FEMTO*farad )
    end associate

    associate (henry => volt * second / ampere)
      call this%unit_const%insert(new_string("H" ),       henry )
      call this%unit_const%insert(new_string("mH"), MILLI*henry )
      call this%unit_const%insert(new_string("uH"), MICRO*henry )
      call this%unit_const%insert(new_string("nH"), NANO *henry )
      call this%unit_const%insert(new_string("pH"), PICO *henry )
      call this%unit_const%insert(new_string("fH"), FEMTO*henry )
    end associate

    call this%unit_const%insert(new_string("eps0"), diel)

    call this%unit_const%insert(new_string("K"), kelvin)
  end subroutine

  subroutine normalization_destruct(this)
    !! destruct normalization constants (release memory)
    class(normalization), intent(inout) :: this

    ! destruct unit constants
    call this%unit_const%destruct()
  end subroutine

end module
