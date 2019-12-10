module normalization_m
  use error_m
  use vector_m
  implicit none

  type normalization
    type(vector_string) :: unit
      !! phyiscal unit tokens (e.g. "eV" or "kV/cm")
    type(vector_real)   :: const
      !! corresponding normalization constants
  contains
    procedure :: init        => normalization_init
    procedure :: search_unit => normalization_search_unit
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
      !! temperature

    call normconst%init(T)
  end subroutine

  function norm_scalar(value, unit, n) result(nvalue)
    !! normalize scalar

    real,                          intent(in) :: value
      !! value to normalize
    character(len=*),              intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalue
      !! return normalized value

    ! local variables
    integer :: i

    nvalue = 0

    ! search for unit and normalize value
    if (present(n)) then
      i = n%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalue = value / n%const%d(i)
    else
      i = normconst%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalue = value / normconst%const%d(i)
    end if
  end function

  function denorm_scalar(value, unit, n) result(nvalue)
    !! denormalize scalar

    real,                          intent(in) :: value
      !! value to denormalize
    character(len=*),              intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalue
      !! return denormalized value

    ! local variables
    integer :: i

    nvalue = 0

    ! search for unit and denormalize value
    if (present(n)) then
      i = n%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalue = value * n%const%d(i)
    else
      i = normconst%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalue = value * normconst%const%d(i)
    end if
  end function

  function norm_array(values, unit, n) result(nvalues)
    !! normalize array

    real,                          intent(in) :: values(:)
      !! values to normalize
    character(len=*),              intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalues(size(values))
      !! return normalized values

    ! local variables
    integer :: i

    nvalues = 0

    ! search for unit and normalize values
    if (present(n)) then
      i = n%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalues = values / n%const%d(i)
    else
      i = normconst%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalues = values / normconst%const%d(i)
    end if
  end function

  function denorm_array(values, unit, n) result(nvalues)
    !! denormalize array

    real,                          intent(in) :: values(:)
      !! values to denormalize
    character(len=*),              intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalues(size(values))
      !! return denormalized values

    ! local variables
    integer :: i

    nvalues = 0

    ! search for unit and normalize values
    if (present(n)) then
      i = n%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalues = values * n%const%d(i)
    else
      i = normconst%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalues = values * normconst%const%d(i)
    end if
  end function

  function norm_matrix(values, unit, n) result(nvalues)
    !! normalize matrix

    real,                          intent(in) :: values(:,:)
      !! values to normalize
    character(len=*),              intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalues(size(values,1),size(values,2))
      !! return normalized values

    ! local variables
    integer :: i

    nvalues = 0

    ! search for unit and normalize values
    if (present(n)) then
      i = n%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalues = values / n%const%d(i)
    else
      i = normconst%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalues = values / normconst%const%d(i)
    end if
  end function

  function denorm_matrix(values, unit, n) result(nvalues)
    !! denormalize matrix

    real,                          intent(in) :: values(:,:)
      !! values to denormalize
    character(len=*),              intent(in) :: unit
      !! physical unit token
    type(normalization), optional, intent(in) :: n
      !! optional normalization object (default: use global normconst)
    real                                      :: nvalues(size(values,1),size(values,2))
      !! return denormalized values

    ! local variables
    integer :: i

    nvalues = 0

    ! search for unit and normalize values
    if (present(n)) then
      i = n%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalues = values * n%const%d(i)
    else
      i = normconst%search_unit(unit)
      if (i < 1) then
        print *, trim(unit)
        call program_error("Unit not found in normalization object!")
      end if
      nvalues = values * normconst%const%d(i)
    end if
  end function

  subroutine normalization_init(this, T)
    !! initialize normalization constants

    class(normalization), intent(out) :: this
    real,                 intent(in)  :: T
      !! temperature

    ! constants
    real, parameter :: PI     = 3.141592653589793238462643
    real, parameter :: EC     = 1.602176462e-19                        ! Electron charge      [ As ]
    real, parameter :: EM     = 9.10938188e-31                         ! Electron rest mass   [ kg ]
    real, parameter :: PLANCK = 6.58211889e-16                         ! Planck's constant    [ eVs ]
    real, parameter :: CLIGHT = 2.99792458d8                           ! Speed of light       [ m/s ]
    real, parameter :: AVOGA  = 6.02214199e23                          ! Avogadro's constant  [ 1/mol ]
    real, parameter :: BOLTZ  = 8.617342e-5                            ! Boltzmann's constant [ eV/K ]
    real, parameter :: EPS0   = 8.854187817e-12                        ! Dielectric constant  [ As/(Vm) ]
    real, parameter :: FSC    = EC / (PLANCK * CLIGHT * 4 * PI * EPS0) ! Fine structure constant ~ 1/137

    ! local variables
    real :: eV, hq, rmom, spr, spk, time, velo, pot, field, conc, dpc, scrt, curr, resist, cap, ind, cvr, dk, mass, &
      &     norm_permittivity, mass_dens

    ! init norm constants
    mass   = EM              ! mass normed to em0 [kg]
    eV     = BOLTZ*T         ! energy [eV]
    hq     = PLANCK          ! Planck's constant [eVs]
    rmom   = dsqrt(EM/EC*eV) ! momentum [eVs/m]
    spr    = hq/rmom         ! r-space [m]
    spk    = 1d0/spr         ! k-space [1/m]
    time   = hq/eV           ! time    [s]
    velo   = spr/time        ! velocity [m/s]
    pot    = eV              ! el. potential [V]
    field  = pot/spr         ! el. field [V/m]
    conc   = spk*spk*spk     ! concentration [1/m**3]
    dpc    = field           ! deformation potential constant [eV/m]
    scrt   = 1d0/time        ! scattering rate [1/s]
    curr   = EC/(time*spr)   ! current [A/m]
    resist = eV*time/EC      ! resistance [Ohm = V/A]
    cap    = EC/eV           ! capacitance [As/V]
    ind    = resist * time   ! inductance [Ohm*s]
    cvr    = CLIGHT/velo
    dk     = 4d0*PI*FSC*cvr  ! dielectric constant (eps_vacuum /= 1)
    norm_permittivity = EC / eV / spr
    mass_dens         =  mass*conc      ! mass density  [kg/m**3]

    ! initialize vectors
    call this%unit%init( n = 0, c = 64)
    call this%const%init(n = 0, c = 64)

    ! insert units
    call insert("1",                    1.0)
    call insert("m",                    1.0 / spr)
    call insert("m^2",                  1.0 / spr / spr)
    call insert("m^3",                  1.0 / spr / spr / spr)
    call insert("cm",                   1e-2 / spr)
    call insert("um",                   1e-6 / spr)
    call insert("um^2",                 1e-12 / spr / spr)
    call insert("um^3",                 1e-18 / spr / spr / spr)
    call insert("nm",                   1e-9 / spr)
    call insert("1/cm^3",               1e6 / conc)
    call insert("1/cm^2",               1e4 / spk / spk)
    call insert("1/cm",                 1e2 / spk)
    call insert("eV",                   1.0 / eV)
    call insert("meV",                  1e-3 / eV)
    call insert("eV/cm",                1e2 / eV / spk)
    call insert("eV/cm^3",              1e6 / ev / conc)
    call insert("eV/cm^2/s",            1e4 / ev / spk / spk / scrt)
    call insert("1/cm^3/eV",            1e6 / conc * eV)
    call insert("1/eV",                 1e0 * eV)
    call insert("V",                    1.0 / pot)
    call insert("mV",                   1e-3 / pot)
    call insert("V/m",                  1.0 / field)
    call insert("V/cm",                 1e2 / field)
    call insert("kV/cm",                1e5 / field)
    call insert("1/V",                  pot)
    call insert("kg",                   1e0 / mass)
    call insert("s",                    1.0 / time)
    call insert("fs",                   1e-15 / time)
    call insert("1/s",                  1.0 / scrt)
    call insert("Hz",                   1.0 / scrt)
    call insert("cm^2/V/s",             1e-4 / (velo / field))
    call insert("V/s",                  1.0 / (pot / time))
    call insert("m/s",                  1.0 / velo)
    call insert("cm/s",                 1e-2 / velo)
    call insert("(V/cm)^(-2/3)*K*cm/s", (10.0 ** (-10.0 / 3.0)) / (field ** (- 2.0 / 3.0) * T * velo))
    call insert("A",                    1.0 / curr / spr)
    call insert("As",                   1.0 / EC)
    call insert("A/m",                  1.0 / curr)
    call insert("A/cm",                 1e2 / curr)
    call insert("A/cm/V",               1e2 / curr * pot)
    call insert("A/cm^2",               1e4 / curr * spr)
    call insert("mA/um^2",              1e9 / curr * spr)
    call insert("1/cm^2/s",             1e4 / spk / spk / scrt)
    call insert("1/um^2/s",             1e12 / spk / spk / scrt)
    call insert("A/V",                  resist)
    call insert("Ohm",                  1.0 / resist)
    call insert("V/A",                  1.0 / resist)
    call insert("F",                    1.0 / cap)
    call insert("s/V",                  1.0 / time * pot)
    call insert("H",                    1.0 / ind)
    call insert("eps0",                 1.0 / dk)
    call insert("C/V/m",                1e0 / norm_permittivity)
    call insert("g/cm^3",               1e3 / mass_dens)

  contains
    subroutine insert(unit, const)
      character(len=*), intent(in) :: unit
      real,             intent(in) :: const
      type(string)                 :: s

      s%s = unit

      call this%unit%push(s)
      call this%const%push(const)
    end subroutine
  end subroutine

  function normalization_search_unit(this, unit) result(i)
    !! Search for physical unit token and return its index

    class(normalization), intent(in) :: this
    character(len=*),     intent(in) :: unit
      !! physical unit token to search for
    integer                          :: i
      !! if found: return index; if not found: return -1

    do i = 1, this%unit%n
      if (this%unit%d(i)%s == unit) return
    end do

    ! not found
    i = -1
  end function

end module
