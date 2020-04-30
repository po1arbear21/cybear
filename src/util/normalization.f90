#include "assert.f90.inc"

module normalization_m
  use error_m
  use map_m
  use math_m
  implicit none

  private
  public normalization, norm, denorm, init_normconst, destruct_normconst

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
    module procedure :: norm_scalar
    module procedure :: norm_array
    module procedure :: norm_matrix
  end interface

  interface denorm
    module procedure :: denorm_scalar
    module procedure :: denorm_array
    module procedure :: denorm_matrix
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

  function norm_scalar(value, unit, n) result(nvalue)
    !! normalize scalar

    real,                          intent(in) :: value
      !! value to normalize
    character(*),                  intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalue
      !! return normalized value

    nvalue = value / get_norm_value(unit, n)
  end function

  function denorm_scalar(value, unit, n) result(nvalue)
    !! denormalize scalar

    real,                          intent(in) :: value
      !! value to denormalize
    character(*),                  intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalue
      !! return denormalized value

    nvalue = value * get_norm_value(unit, n)
  end function

  function norm_array(values, unit, n) result(nvalues)
    !! normalize array

    real,                          intent(in) :: values(:)
      !! values to normalize
    character(*),                  intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalues(size(values))
      !! return normalized values

    nvalues = values / get_norm_value(unit, n)
  end function

  function denorm_array(values, unit, n) result(nvalues)
    !! denormalize array

    real,                          intent(in) :: values(:)
      !! values to denormalize
    character(*),                  intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalues(size(values))
      !! return denormalized values

    nvalues = values * get_norm_value(unit, n)
  end function

  function norm_matrix(values, unit, n) result(nvalues)
    !! normalize matrix

    real,                          intent(in) :: values(:,:)
      !! values to normalize
    character(*),                  intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalues(size(values,1),size(values,2))
      !! return normalized values

    nvalues = values / get_norm_value(unit, n)
  end function

  function denorm_matrix(values, unit, n) result(nvalues)
    !! denormalize matrix

    real,                          intent(in) :: values(:,:)
      !! values to denormalize
    character(*),                  intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalues(size(values,1),size(values,2))
      !! return denormalized values

    nvalues = values * get_norm_value(unit, n)
  end function

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
      node => n%unit_const%find(string(unit))
    else
      node => normconst%unit_const%find(string(unit))
    end if
    if (.not. associated(node)) then
      print *, "Unit: " // trim(unit)
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
    real, parameter :: EC     = 1.602176462e-19                        ! Electron charge           [ As      ]
    real, parameter :: EM     = 9.10938188e-31                         ! Electron rest mass        [ kg      ]
    real, parameter :: PLANCK = 6.58211889e-16                         ! reduced Planck's constant [ eVs     ]
    real, parameter :: BOLTZ  = 8.617333262145e-5                      ! Boltzmann's constant      [ eV/K    ]
    real, parameter :: EPS0   = 8.854187817e-12                        ! Vacuum permittivity       [ As/(Vm) ]

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

    call this%unit_const%insert(string("1"), 1.0)

    call this%unit_const%insert(string("m"   ), meter          )
    call this%unit_const%insert(string("m^2" ), meter**2       )
    call this%unit_const%insert(string("m^3" ), meter**3       )
    call this%unit_const%insert(string("cm"  ), meter    * 1e2 )
    call this%unit_const%insert(string("cm^2"), meter**2 * 1e4 )
    call this%unit_const%insert(string("cm^3"), meter**3 * 1e6 )
    call this%unit_const%insert(string("mm"  ), meter    * 1e3 )
    call this%unit_const%insert(string("mm^2"), meter**2 * 1e6 )
    call this%unit_const%insert(string("mm^3"), meter**3 * 1e9 )
    call this%unit_const%insert(string("um"  ), meter    * 1e6 )
    call this%unit_const%insert(string("um^2"), meter**2 * 1e12)
    call this%unit_const%insert(string("um^3"), meter**3 * 1e18)
    call this%unit_const%insert(string("nm"  ), meter    * 1e9 )
    call this%unit_const%insert(string("nm^2"), meter**2 * 1e18)
    call this%unit_const%insert(string("nm^3"), meter**3 * 1e27)
    call this%unit_const%insert(string("pm"  ), meter    * 1e12)
    call this%unit_const%insert(string("pm^2"), meter**2 * 1e24)
    call this%unit_const%insert(string("pm^3"), meter**3 * 1e36)
    call this%unit_const%insert(string("fm"  ), meter    * 1e15)
    call this%unit_const%insert(string("fm^2"), meter**2 * 1e30)
    call this%unit_const%insert(string("fm^3"), meter**3 * 1e45)

    call this%unit_const%insert(string("1/m"   ), 1.0 / (meter          ))
    call this%unit_const%insert(string("1/m^2" ), 1.0 / (meter**2       ))
    call this%unit_const%insert(string("1/m^3" ), 1.0 / (meter**3       ))
    call this%unit_const%insert(string("1/cm"  ), 1.0 / (meter    * 1e2 ))
    call this%unit_const%insert(string("1/cm^2"), 1.0 / (meter**2 * 1e4 ))
    call this%unit_const%insert(string("1/cm^3"), 1.0 / (meter**3 * 1e6 ))
    call this%unit_const%insert(string("1/mm"  ), 1.0 / (meter    * 1e3 ))
    call this%unit_const%insert(string("1/mm^2"), 1.0 / (meter**2 * 1e6 ))
    call this%unit_const%insert(string("1/mm^3"), 1.0 / (meter**3 * 1e9 ))
    call this%unit_const%insert(string("1/um"  ), 1.0 / (meter    * 1e6 ))
    call this%unit_const%insert(string("1/um^2"), 1.0 / (meter**2 * 1e12))
    call this%unit_const%insert(string("1/um^3"), 1.0 / (meter**3 * 1e18))
    call this%unit_const%insert(string("1/nm"  ), 1.0 / (meter    * 1e9 ))
    call this%unit_const%insert(string("1/nm^2"), 1.0 / (meter**2 * 1e18))
    call this%unit_const%insert(string("1/nm^3"), 1.0 / (meter**3 * 1e27))
    call this%unit_const%insert(string("1/pm"  ), 1.0 / (meter    * 1e12))
    call this%unit_const%insert(string("1/pm^2"), 1.0 / (meter**2 * 1e24))
    call this%unit_const%insert(string("1/pm^3"), 1.0 / (meter**3 * 1e36))
    call this%unit_const%insert(string("1/fm"  ), 1.0 / (meter    * 1e15))
    call this%unit_const%insert(string("1/fm^2"), 1.0 / (meter**2 * 1e30))
    call this%unit_const%insert(string("1/fm^3"), 1.0 / (meter**3 * 1e45))

    call this%unit_const%insert(string("s" ), second       )
    call this%unit_const%insert(string("ms"), second * 1e3 )
    call this%unit_const%insert(string("us"), second * 1e6 )
    call this%unit_const%insert(string("ns"), second * 1e9 )
    call this%unit_const%insert(string("ps"), second * 1e12)
    call this%unit_const%insert(string("fs"), second * 1e15)

    call this%unit_const%insert(string("1/s" ), 1.0 / (second       ))
    call this%unit_const%insert(string("1/ms"), 1.0 / (second * 1e3 ))
    call this%unit_const%insert(string("1/us"), 1.0 / (second * 1e6 ))
    call this%unit_const%insert(string("1/ns"), 1.0 / (second * 1e9 ))
    call this%unit_const%insert(string("1/ps"), 1.0 / (second * 1e12))
    call this%unit_const%insert(string("1/fs"), 1.0 / (second * 1e15))

    call this%unit_const%insert(string("Hz" ), 1.0 / (second       ))
    call this%unit_const%insert(string("kHz"), 1.0 / (second * 1e3 ))
    call this%unit_const%insert(string("MHz"), 1.0 / (second * 1e6 ))
    call this%unit_const%insert(string("GHz"), 1.0 / (second * 1e9 ))
    call this%unit_const%insert(string("THz"), 1.0 / (second * 1e12))

    call this%unit_const%insert(string("m/s" ), meter / second      )
    call this%unit_const%insert(string("cm/s"), meter / second * 1e2)

    call this%unit_const%insert(string("1/cm^2/s"), 1.0 / (meter**2 * 1e2  * second))
    call this%unit_const%insert(string("1/um^2/s"), 1.0 / (meter**2 * 1e12 * second))

    call this%unit_const%insert(string("kg"    ), kilogram                         )
    call this%unit_const%insert(string("kg/m^3"), kilogram / meter**3              )
    call this%unit_const%insert(string("g/cm^3"), kilogram * 1e3 / (meter**3 * 1e6))

    call this%unit_const%insert(string("V"    ), volt                       )
    call this%unit_const%insert(string("mV"   ), volt * 1e3                 )
    call this%unit_const%insert(string("V/m"  ), volt / meter               )
    call this%unit_const%insert(string("V/cm" ), volt / (meter * 1e2)       )
    call this%unit_const%insert(string("kV/cm"), volt * 1e-3 / (meter * 1e2))
    call this%unit_const%insert(string("1/V"  ), 1.0 / volt                 )

    call this%unit_const%insert(string("eV"       ), volt                         )
    call this%unit_const%insert(string("meV"      ), volt * 1e3                   )
    call this%unit_const%insert(string("eV/cm"    ), volt / (meter * 1e2)         )
    call this%unit_const%insert(string("eV/cm^3"  ), volt / (meter * 1e4)         )
    call this%unit_const%insert(string("eV/cm^2/s"), volt / (meter * 1e2 * second))
    call this%unit_const%insert(string("1/eV"     ), 1.0 / volt                   )
    call this%unit_const%insert(string("1/cm^3/eV"), 1.0 / (meter**3 * 1e6 * volt))

    call this%unit_const%insert(string("cm^2/V/s"            ), meter**2*1e4/volt/second                              )
    call this%unit_const%insert(string("V/s"                 ), volt/second                                           )
    call this%unit_const%insert(string("s/V"                 ), second/volt                                           )
    call this%unit_const%insert(string("(V/cm)^(-2/3)*K*cm/s"), (volt/(meter*1e2))**(-2.0/3.0)*kelvin*meter*1e2/second)

    call this%unit_const%insert(string("A"      ), ampere)
    call this%unit_const%insert(string("As"     ), ampere * second)
    call this%unit_const%insert(string("A/m"    ), ampere / meter)
    call this%unit_const%insert(string("A/cm"   ), ampere / (meter * 1e2))
    call this%unit_const%insert(string("A/cm/V" ), ampere / (meter * 1e2 * volt))
    call this%unit_const%insert(string("A/cm^2" ), ampere / (meter**2 * 1e4))
    call this%unit_const%insert(string("mA/um^2"), ampere * 1e3 / (meter**2 * 1e12))

    call this%unit_const%insert(string("A/V" ), ampere / volt)
    call this%unit_const%insert(string("V/A" ), volt / ampere)
    call this%unit_const%insert(string("Ohm" ), volt / ampere)
    call this%unit_const%insert(string("mOhm"), volt / ampere * 1e3 )
    call this%unit_const%insert(string("uOhm"), volt / ampere * 1e6 )
    call this%unit_const%insert(string("nOhm"), volt / ampere * 1e9 )
    call this%unit_const%insert(string("pOhm"), volt / ampere * 1e12)
    call this%unit_const%insert(string("fOhm"), volt / ampere * 1e15)
    call this%unit_const%insert(string("kOhm"), volt / ampere * 1e-3)
    call this%unit_const%insert(string("MOhm"), volt / ampere * 1e-6)

    call this%unit_const%insert(string("F" ),  ampere * second / volt       )
    call this%unit_const%insert(string("mF"),  ampere * second / volt * 1e3 )
    call this%unit_const%insert(string("uF"),  ampere * second / volt * 1e6 )
    call this%unit_const%insert(string("nF"),  ampere * second / volt * 1e9 )
    call this%unit_const%insert(string("pF"),  ampere * second / volt * 1e12)
    call this%unit_const%insert(string("fF"),  ampere * second / volt * 1e15)

    call this%unit_const%insert(string("H" ), volt * second / ampere       )
    call this%unit_const%insert(string("mH"), volt * second / ampere * 1e3 )
    call this%unit_const%insert(string("uH"), volt * second / ampere * 1e6 )
    call this%unit_const%insert(string("nH"), volt * second / ampere * 1e9 )
    call this%unit_const%insert(string("pH"), volt * second / ampere * 1e12)
    call this%unit_const%insert(string("fH"), volt * second / ampere * 1e15)

    call this%unit_const%insert(string("eps0"), diel)
  end subroutine

  subroutine normalization_destruct(this)
    !! destruct normalization constants (release memory)
    class(normalization), intent(inout) :: this

    ! destruct unit constants
    call this%unit_const%destruct()
  end subroutine

end module
