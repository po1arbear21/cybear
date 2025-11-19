m4_include(util/macro.f90.inc)

module semiconductor_m

  use distribution_table_m, only: distribution_table
  use error_m,              only: program_error
  use fukushima_m,          only: fd1h, fdm1h, fdm3h, fdm5h, fdm7h, ifd1h
  use ieee_arithmetic,      only: ieee_is_finite, ieee_value, ieee_positive_inf, ieee_negative_inf
  use math_m,               only: PI, expm1
  use newton_m,             only: newton1D_opt, newton1D
  use util_m,               only: int2str

  implicit none

  integer,      parameter :: CR_ELEC      = 1
  integer,      parameter :: CR_HOLE      = 2
  character(*), parameter :: CR_NAME(2)   = [ "n", "p" ]
  real,         parameter :: CR_CHARGE(2) = [ -1.0, 1.0 ]

  integer,      parameter :: DOP_DCON      = 1
  integer,      parameter :: DOP_ACON      = 2
  character(*), parameter :: DOP_NAME(2)   = [ "D", "A" ]
  real,         parameter :: DOP_CHARGE(2) = [ 1.0, -1.0]

  integer,      parameter :: DOS_PARABOLIC      = 1
    !! Parabolic band structure
  integer,      parameter :: DOS_PARABOLIC_TAIL = 2
    !! Parabolic band structure with exponential tail
  integer,      parameter :: DOS_GAUSS          = 3

  integer,      parameter :: DIST_MAXWELL   = 1
    !! Maxwell-Boltzmann distribution density
  integer,      parameter :: DIST_FERMI     = 2
    !! Fermi-Dirac distribution density
  integer,      parameter :: DIST_FERMI_REG = 3
    !! Use Fermi-Dirac integral with regularization directly


  type semiconductor
    real    :: edos(2)
      !! effective density of states Nc, Nv
    real    :: band_gap
      !! band gap
    real    :: band_edge(2)
      !! conduction/valence band edge

    integer                  :: dos
      !! density of states (DOS_PARABOLIC or DOS_PARABOLIC_TAIL)
    real                     :: dos_params(1)
      !! parameters for density of states
    integer                  :: dist
      !! distribution density (DIST_MAXWELL, DIST_FERMI or DIST_FERMI_REG)
    real                     :: dist_params(2)
      !! parameters for distribution
    type(distribution_table) :: dist_tab
      !! get cumulative distribution from table

    logical           :: mob
      !! enable/disable mobility saturation
    real, allocatable :: alpha(:)
      !! Caughey-Thomas alpha parameter
    real, allocatable :: beta(:)
      !! Caughey-Thomas beta parameter
    real, allocatable :: mob_min(:)
      !! Caughey-Thomas minimal mobility
    real, allocatable :: mob_max(:)
      !! Caughey-Thomas maximal mobility
    real, allocatable :: N_ref(:)
      !! Caughey-Thomas reference density
    real, allocatable :: v_sat(:)
      !! Caughey-Thomas saturation velocity

    logical           :: incomp_ion
      !! enable/disable incomplete ionization (Altermatt-Schenk model)
    real, allocatable :: ii_tau(:)
      !! generation/recombination time constant (carrier index)
    real, allocatable :: ii_E_dop0(:)
      !! dopant energy relative to the carrier band for a single dopant
    real, allocatable :: ii_g(:)
      !! degeneracy factor of dopant states
    real, allocatable :: ii_N_crit(:)
      !! cricital concentration (alternative to altermatt-schenk)
    real, allocatable :: ii_dop_th(:)
      !! full ionization if doping > threshold

  contains
    procedure :: init_dist => semiconductor_init_dist
    procedure :: get_dist  => semiconductor_get_dist
    procedure :: get_idist => semiconductor_get_idist
  end type

contains

  subroutine semiconductor_init_dist(this)
    class(semiconductor), intent(inout) :: this

    if (this%dos == DOS_PARABOLIC) return

    call this%dist_tab%init("dist_" // int2str(this%dos) // int2str(this%dist), dos, ieee_value(1.0, ieee_negative_inf), ieee_value(1.0, ieee_positive_inf), dist, .false., -100.0, 1000.0, 3)

  contains

    function dos(t) result(Z)
      !! get density of states
      real(kind=16), intent(in) :: t
        !! energy relative to band edge (in units of k_B T)
      real(kind=16)             :: Z
        !! return density of states

      real :: t0, sigma

      select case (this%dos)
      case (DOS_PARABOLIC)
        if (t > 0) then
          Z = 2 * sqrt(t / PI)
        else
          Z = 0
        end if

      case (DOS_PARABOLIC_TAIL)
        t0 = this%dos_params(1)
        if (t > t0) then
          Z = 2 * sqrt(t / PI)
        else
          Z = 2 * sqrt(t0 / PI) * exp(0.5 * (t - t0) / t0)
        end if

      case (DOS_GAUSS)
        sigma = this%dos_params(1)
        Z = 1 / (sqrt(2 * PI) * sigma) * exp(- (t / sigma)**2 / 2)

      end select
    end function

    function dist(u, k) result(f)
      !! k-th derivative of distribution density (e.g. fermi-dirac or maxwell-boltzmann)
      real(kind=16),    intent(in) :: u
        !! energy relative to chemical potential/fermi level (in units of k_B T)
      integer,          intent(in) :: k
        !! k-th derivative (possible values from 0 to kmax)
      real(kind=16)                :: f
        !! return k-th distribution density

      real(kind=16) :: e, f0

      select case (this%dist)
      case (DIST_MAXWELL)
        f = exp(-u)
        if (mod(k, 2) /= 0) f = - f

      case (DIST_FERMI)
        e  = exp(u)
        f0 = 1 / (1 + e)
        f  = 0

        select case (k)
        case (0)
          f = f0
        case (1)
          if (ieee_is_finite(e)) f = (- e * f0) * f0
        case (2)
          if (ieee_is_finite(e)) f = ((e * f0) * (expm1(u) * f0)) * f0
        case (3)
          if (ieee_is_finite(e)) f = (e * f0) * (- 1 + 6 * (e * f0) * (1 - e * f0)) * f0
        case (4)
          if (ieee_is_finite(e)) f = (- (e * f0) + 14 * (e * f0)**2 - 36 * (e * f0)**3 + 24 * (e * f0)**4) * f0
        end select
      end select
    end function

  end subroutine

  subroutine semiconductor_get_dist(this, eta, k, F, dFdeta)
    !! get k-th derivative cumulative distribution function
    class(semiconductor), intent(in)  :: this
    real,                 intent(in)  :: eta
      !! chemical potential
    integer,              intent(in)  :: k
      !! derivative (can be 0)
    real,                 intent(out) :: F
      !! output cumulative distribution function
    real,                 intent(out) :: dFdeta
      !! output derivative of F wrt eta

    real, parameter :: G0 = 1.0 / gamma(1.5)
    real, parameter :: G1 = 1.0 / gamma(0.5)
    real, parameter :: G2 = 1.0 / gamma(-0.5)
    real, parameter :: G3 = 1.0 / gamma(-1.5)
    real, parameter :: G4 = 1.0 / gamma(-2.5)
    real, parameter :: G5 = 1.0 / gamma(-3.5)

    real :: A, B

    if (this%dos == DOS_PARABOLIC) then
      select case (this%dist)
      case (DIST_MAXWELL)
        ! Maxwell-Boltzmann integral
        F      = exp(eta)
        dFdeta = F

      case (DIST_FERMI)
        ! Fermi-Dirac integral
        select case (k)
        case (0)
          F      = fd1h(eta) * G0
          dFdeta = fdm1h(eta) * G1
        case (1)
          F      = fdm1h(eta) * G1
          dFdeta = fdm3h(eta) * G2
        case (2)
          F      = fdm3h(eta) * G2
          dFdeta = fdm5h(eta) * G3
        case (3)
          F      = fdm5h(eta) * G3
          dFdeta = fdm7h(eta) * G4
        end select

      case (DIST_FERMI_REG)
        ! Fermi-Dirac integral with regularization
        A = this%dist_params(1)
        B = this%dist_params(2)
        select case (k)
        case (0)
          F      = (fd1h( eta) + A*     fd1h(B*eta)) * G0
          dFdeta = (fdm1h(eta) + A*B*   fdm1h(B*eta)) * G1
        case (1)
          F      = (fdm1h(eta) + A*B*   fdm1h(B*eta)) * G1
          dFdeta = (fdm3h(eta) + A*B**2*fdm3h(B*eta)) * G2
        case (2)
          F      = (fdm3h(eta) + A*B**2*fdm3h(B*eta)) * G2
          dFdeta = (fdm5h(eta) + A*B**3*fdm5h(B*eta)) * G3
        case (3)
          F      = (fdm5h(eta) + A*B**3*fdm5h(B*eta)) * G3
          dFdeta = (fdm7h(eta) + A*B**4*fdm7h(B*eta)) * G4
        end select

      end select
    else
      ! use table
      call this%dist_tab%get(eta, k, F, dFdeta)
    end if
  end subroutine

  subroutine semiconductor_get_idist(this, F, eta, detadF)
    !! get inverse of cumulative distribution function
    class(semiconductor), intent(in)  :: this
    real,                 intent(in)  :: F
      !! value of cumulative distribution function
    real,                 intent(out) :: eta
      !! output corresponding chemical potential
    real,                 intent(out) :: detadF
      !! output derivative of eta wrt F

    real, parameter :: G = gamma(1.5)
    real, parameter :: C = 6.6e-11

    real               :: eta0, p(1), detadp(1)
    type(newton1D_opt) :: nopt

    ! safety check
    if (F <= 0) then
      print "(A,ES25.16E3)", "F = ", F
      call program_error("F must be positive")
    end if

    if (this%dos == DOS_PARABOLIC) then
      select case (this%dist)
      case (DIST_MAXWELL)
        eta    = log(F)
        detadF = 1 / F

      case (DIST_FERMI)
        call ifd1h(F * G, eta, detadF)
        detadF = detadF * G

      case (DIST_FERMI_REG)
        ! initial guess
        if (F < C) then
          eta0 = log(F / this%dist_params(1)) / this%dist_params(2)
        else
          call ifd1h(F * G, eta0, detadF)
        end if

        ! solve with Newton
        call nopt%init()
        p(1) = F
        call newton1D(fun, p, nopt, eta0, eta, dxdp = detadp)
        detadF = detadp(1)

      end select
    else
      ! use table
      call this%dist_tab%inv(F, eta, detadF)
    end if

    ! second safety check
    if(.not. (ieee_is_finite(eta) .and. ieee_is_finite(detadF))) then
      print "(A,ES25.16E3)", "eta    = ", eta
      print "(A,ES25.16E3)", "detadF = ", detadF
      call program_error("eta or detadF not finite")
    end if

  contains

    subroutine fun(eta, p, res, dresdeta, dresdp)
      real,              intent(in)  :: eta
        !! argument
      real,              intent(in)  :: p(:)
        !! parameters: F
      real,              intent(out) :: res
        !! output function value
      real,    optional, intent(out) :: dresdeta
        !! optional output derivative of f wrt x
      real,    optional, intent(out) :: dresdp(:)
        !! optional output derivatives of f wrt p

      real :: dresdeta_

      call this%get_dist(eta, 0, res, dresdeta_)
      res       = res - p(1)
      if (present(dresdeta)) dresdeta = dresdeta_
      if (present(dresdp)) dresdp(1) = - 1
    end subroutine

  end subroutine

end module
