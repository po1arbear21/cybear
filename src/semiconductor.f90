m4_include(util/macro.f90.inc)

module semiconductor_m

  use, intrinsic :: iso_fortran_env, only: real128
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_positive_inf, ieee_negative_inf

  use lookup_table_m, only: par_int_table, def_int_table
  use error_m,        only: program_error
  use fukushima_m,    only: fd1h, fdm1h, fdm3h, fdm5h, fdm7h, ifd1h
  use math_m,         only: expm1
  use newton_m,       only: newton_opt, newton
  use util_m,         only: int2str

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

  integer,      parameter :: STAB_SG    = 1
    !! Scharfetter-Gummel (not thermodynamically consistent)
  integer,      parameter :: STAB_ED    = 2
    !! enhanced diffusion (thermodynamically consistent)
  integer,      parameter :: STAB_MED   = 3
    !! modified enhanced diffusion, more accurate (thermodynamically consistent)
  integer,      parameter :: STAB_EXACT = 4
    !! exact (numerical) solution of integral equation (thermodynamically consistent)

  type semiconductor
    real :: edos(2)
      !! effective density of states Nc, Nv
    real :: band_gap
      !! band gap
    real :: band_edge(2)
      !! conduction/valence band edge

    integer                  :: dos
      !! density of states (DOS_PARABOLIC or DOS_PARABOLIC_TAIL)
    real                     :: dos_params(1)
      !! parameters for density of states
    integer                  :: dist
      !! distribution density (DIST_MAXWELL, DIST_FERMI or DIST_FERMI_REG)
    real                     :: dist_params(2)
      !! parameters for distribution
    integer                  :: stab
      !! stabilization method for degenerate case (STAB_SG, STAB_ED, STAB_ED2, STAB_EXACT)
    type(par_int_table)      :: Ftab(0:4)
      !! distribution function + derivatives
    type(def_int_table)      :: Itab(-5:5)
      !! integral_0^{eta} F(eta')^k deta'; k = -5 to 5

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

    logical :: srh = .false.
      !! enable/disable SRH recombination
    real    :: srh_tau_n = 0.0
      !! SRH electron lifetime (normalized)
    real    :: srh_tau_p = 0.0
      !! SRH hole lifetime (normalized)

  contains
    procedure :: init_dist    => semiconductor_init_dist
    generic   :: get_dist     => semiconductor_get_dist, &
      &                          semiconductor_get_dist_16
    generic   :: get_inv_dist => semiconductor_get_inv_dist, &
      &                          semiconductor_get_inv_dist_16
    generic   :: get_int_dist => semiconductor_get_int_dist, &
      &                          semiconductor_get_int_dist_16

    procedure, private :: semiconductor_get_dist
    procedure, private :: semiconductor_get_dist_16
    procedure, private :: semiconductor_get_inv_dist
    procedure, private :: semiconductor_get_inv_dist_16
    procedure, private :: semiconductor_get_int_dist
    procedure, private :: semiconductor_get_int_dist_16
  end type

contains

  subroutine semiconductor_init_dist(this, tabledir)
    class(semiconductor), intent(inout) :: this
    character(*),         intent(in)    :: tabledir
      !! file system table directory

    real(real128), parameter :: ETA0 = -5e3_real128, ETA1 = 5e3_real128
    real(real128), parameter :: RTOL = 5e-16_real128, ATOL = 1e-24_real128, XTOL = 1e-8_real128

    character(255)            :: project_root
    character(:), allocatable :: fname, fdir
    integer                   :: k, i, j
    logical                   :: exist
    real(real128)             :: inf, p(1)

    inf = ieee_value(inf, IEEE_POSITIVE_INF)

    ! get absolute path to table directory
    ! FIXME: use environment variable provided by fargo (not yet implemented)
    call get_environment_variable("PWD", project_root)
    i = len_trim(project_root)
    do
      inquire (file = project_root(1:i) // "/fargo.toml", exist = exist)
      if (exist) exit
      do j = i-1, 1, -1
        if (project_root(j:j) == '/') exit
      end do
      i = j
      if (j < 1) call program_error("fargo.toml not found, could not determine project root dir")
    end do
    fdir = project_root(1:i) // '/' // tabledir
    call system("mkdir -p " // fdir)

    do k = lbound(this%Ftab, 1), ubound(this%Ftab, 1)
      fname = "dist_" // int2str(this%dos) // int2str(this%dist) // "_" // int2str(k)
      if (this%dos == DOS_PARABOLIC) then
        call this%Ftab(k)%init(fdir, fname, ETA0, ETA1, Z_dist, 0.0_real128, inf, RTOL, ATOL, XTOL)
      else
        call this%Ftab(k)%init(fdir, fname, ETA0, ETA1, Z_dist, -inf, inf, RTOL, ATOL, XTOL)
      end if
    end do

    do k = -5, 5
      if (k == 0) cycle
      fname = "I_" // int2str(this%dos) // int2str(this%dist) // "_" // int2str(k)
      p(1) = k
      call this%Itab(k)%init(fdir, fname, ETA0, ETA1, 3, Fk, p, RTOL, ATOL, XTOL)
    end do

  contains

    subroutine Z_dist(t, p, Zf, dZfdt, dZfdp)
      real(real128), intent(in)  :: t
      real(real128), intent(in)  :: p(:)
      real(real128), intent(out) :: Zf
      real(real128), intent(out) :: dZfdt
      real(real128), intent(out) :: dZfdp(:)

      real(real128), parameter :: PI16 = 4 * atan(1.0_real128)

      real          :: t0, sigma
      real(real128) :: Z, f(0:1)

      ! density of states
      select case (this%dos)
      case (DOS_PARABOLIC)
        if (t > 0) then
          Z = 2 * sqrt(t / PI16)
        else
          Z = 0
        end if

      case (DOS_PARABOLIC_TAIL)
        t0 = this%dos_params(1)
        if (t > t0) then
          Z = 2 * sqrt(t / PI16)
        else
          Z = 2 * sqrt(t0 / PI16) * exp(0.5 * (t - t0) / t0)
        end if

      case (DOS_GAUSS)
        sigma = this%dos_params(1)
        Z = 1 / (sqrt(2 * PI16) * sigma) * exp(- (t / sigma)**2 / 2)

      end select

      ! distribution density
      select case (this%dist)
      case (DIST_MAXWELL)
        f(0) = exp(t - p(1))
        if (mod(k, 2) /= 0) f(0) = - f(0)
        f(1) = - f(0)

      case (DIST_FERMI, DIST_FERMI_REG)
        call fermi_dirac(t - p(1), k, f)
        if (mod(k, 2) == 0) then
          f(1) = - f(1)
        else
          f(0) = - f(0)
        end if

      end select

      ! combine
      Zf       = Z * f(0)
      dZfdt    = 0
      dZfdp(1) = Z * f(1)
    end subroutine

    subroutine Fk(eta, p, F, dFdeta, dFdp)
      real(real128), intent(in)  :: eta
      real(real128), intent(in)  :: p(:)
      real(real128), intent(out) :: F
      real(real128), intent(out) :: dFdeta
      real(real128), intent(out) :: dFdp(:)

      integer :: k

      m4_ignore(dFdp)

      k = nint(p(1))

      call this%get_dist(eta, 0, F, dFdeta)

      if (k /= 1) then
        dFdeta = k * F**(k-1) * dFdeta
        F      = F**k
      end if
    end subroutine

    subroutine fermi_dirac(u, k, f)
      real(real128), intent(in)  :: u
        !! argument u = t - eta
      integer,       intent(in)  :: k
        !! derivative index
      real(real128), intent(out) :: f(0:1)
        !! output f^{(k)} and f^{(k+1)}

      real(real128), parameter :: uc1 = acosh(2.0_real128), uc2 = acosh(5.0_real128)

      integer       :: j
      real(real128) :: c, s1, s2

      do j = 0, 1
        select case (k+j)
        case (0)
          f(j) = 1 / (1 + exp(u))
        case (1)
          f(j) = - 1 / (2 * (cosh(u) + 1))
        case (2)
          if (abs(u) < 1e-16) then
            f(j) = u*(0.125 - u**2 / 24)
          else
            s1 = sinh(0.5 * u)
            s2 = sinh(u)
            if (ieee_is_finite(s1) .and. ieee_is_finite(s2)) then
              f(j) = 2 * s1 * (s1 / s2)**3
            else
              f(j) = 0
            end if
          end if
        case (3)
          s1 = sinh(0.5 * (u + uc1))
          s2 = sinh(0.5 * (u - uc1))
          c  = cosh(0.5 * u)
          if (ieee_is_finite(s1) .and. ieee_is_finite(s2) .and. ieee_is_finite(c)) then
            f(j) = -0.25 * (s1 / c) * (s2 / c) / c**2
          else
            f(j) = 0
          end if
        case (4)
          s1 = sinh(0.5 * (u + uc2))
          s2 = sinh(0.5 * (u - uc2))
          c  = cosh(0.5 * u)
          if (ieee_is_finite(s1) .and. ieee_is_finite(s2) .and. ieee_is_finite(c)) then
            f(j) = 0.25 * (s1 / c) * (s2 / c) * tanh(0.5 * u) / c**2
          else
            f(j) = 0
          end if
        case (5)
          s1 = cosh(u)
          s2 = cosh(2*u)
          c  = cosh(0.5*u)
          if (ieee_is_finite(s1) .and. ieee_is_finite(s2) .and. ieee_is_finite(c)) then
            f(j) = (26 * s1/c - s2/c - 33/c) / 32 / c**5
          else
            f(j) = 0
          end if
        case default
          call program_error("fermi-dirac not implemented for k+j = " // int2str(k + j))
        end select
      end do
    end subroutine

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

    real(real128) :: eta16, F16, dFdeta16

    eta16  = eta
    call this%get_dist(eta16, k, F16, dFdeta16)
    F      = real(F16)
    dFdeta = real(dFdeta16)
  end subroutine

  subroutine semiconductor_get_dist_16(this, eta, k, F, dFdeta)
    !! get k-th derivative cumulative distribution function in quad precision
    class(semiconductor), intent(in)  :: this
    real(real128),        intent(in)  :: eta
      !! chemical potential
    integer,              intent(in)  :: k
      !! derivative (can be 0)
    real(real128),        intent(out) :: F
      !! output cumulative distribution function
    real(real128),        intent(out) :: dFdeta
      !! output derivative of F wrt eta

    real(real128) :: A, B, F2, dF2deta

    if (this%dist == DIST_FERMI_REG) then
      A = this%dist_params(1)
      B = this%dist_params(2)
      call this%Ftab(k)%get(eta, F, dFdeta)
      call this%Ftab(k)%get(B*eta, F2, dF2deta)
      F      = F      + A * F2      * B**k
      dFdeta = dFdeta + A * dF2deta * B**(k+1)
    else
      call this%Ftab(k)%get(eta, F, dFdeta)
    end if
  end subroutine

  subroutine semiconductor_get_inv_dist(this, F, eta, detadF)
    !! get inverse of cumulative distribution function
    class(semiconductor), intent(in)  :: this
    real,                 intent(in)  :: F
      !! value of cumulative distribution function
    real,                 intent(out) :: eta
      !! output corresponding chemical potential
    real,                 intent(out) :: detadF
      !! output derivative of eta wrt F

    real(real128) :: F16, eta16, detadF16

    F16    = F
    call this%get_inv_dist(F16, eta16, detadF16)
    eta    = real(eta16)
    detadF = real(detadF16)
  end subroutine

  subroutine semiconductor_get_inv_dist_16(this, F, eta, detadF)
    !! get inverse of cumulative distribution function in quad precision
    class(semiconductor), intent(in)  :: this
    real(real128),        intent(in)  :: F
      !! value of cumulative distribution function
    real(real128),        intent(out) :: eta
      !! output corresponding chemical potential
    real(real128),        intent(out) :: detadF
      !! output derivative of eta wrt F

    real(real128), parameter :: C = 6.6e-11_real128
    real(real128), parameter :: RTOL = 1e-28_real128, ATOL = 1e-32_real128

    real(real128) :: err, F2, dF2deta

    ! safety check
    if (F <= 0) then
      print "(A,ES41.32E3)", "F = ", F
      call program_error("F must be positive")
    end if

    if (this%dist == DIST_FERMI_REG) then
      ! initial guess
      if (F < C) then
        eta = log(F / this%dist_params(1)) / this%dist_params(2)
      else
        call this%Ftab(0)%inv(F, eta, detadF)
      end if

      ! newton iteration
      err = huge(err)
      do while (err > RTOL)
        call this%get_dist(eta, 0, F2, dF2deta)
        F2  = F2 - F
        err = - F2 / dF2deta
        eta = eta + err
        err = abs(err) / (abs(eta) + ATOL)
      end do
      detadF = 1 / dF2deta
    else
      call this%Ftab(0)%inv(F, eta, detadF)
    end if

    ! second safety check
    if(.not. (ieee_is_finite(eta) .and. ieee_is_finite(detadF))) then
      print "(A,ES41.32E3)", "eta    = ", eta
      print "(A,ES41.32E3)", "detadF = ", detadF
      call program_error("eta or detadF not finite")
    end if
  end subroutine

  subroutine semiconductor_get_int_dist(this, eta, k, I, dIdeta)
    !! integral_eta(1)^eta(2) dist(eta)^k deta
    class(semiconductor), intent(in)  :: this
    real,                 intent(in)  :: eta(2)
      !! integration bounds
    integer,              intent(in)  :: k
      !! exponent
    real,                 intent(out) :: I
      !! output integral over dist^k
    real,                 intent(out) :: dIdeta(2)
      !! output derivatives of I wrt eta

    real(real128) :: eta16(2), I16, dIdeta16(2)

    eta16  = eta
    call this%get_int_dist(eta16, k, I16, dIdeta16)
    I      = real(I16)
    dIdeta = real(dIdeta16)
  end subroutine

  subroutine semiconductor_get_int_dist_16(this, eta, k, I, dIdeta)
    !! integral_eta(1)^eta(2) dist(eta)^k deta in quad precision
    class(semiconductor), intent(in)  :: this
    real(real128),        intent(in)  :: eta(2)
      !! integration bounds
    integer,              intent(in)  :: k
      !! exponent
    real(real128),        intent(out) :: I
      !! output integral over dist^k
    real(real128),        intent(out) :: dIdeta(2)
      !! output derivatives of I wrt eta

    if (k == 0) then
      I      = eta(2) - eta(1)
      dIdeta = [-1, 1]
    else
      call this%Itab(k)%get(Fk, eta(1), eta(2), I, dIdeta(1), dIdeta(2))
    end if

  contains

    subroutine Fk(eta, p, F, dFdeta, dFdp)
      real(real128), intent(in)  :: eta
      real(real128), intent(in)  :: p(:)
      real(real128), intent(out) :: F
      real(real128), intent(out) :: dFdeta
      real(real128), intent(out) :: dFdp(:)

      integer :: k

      m4_ignore(dFdp)

      k = nint(p(1))

      call this%get_dist(eta, 0, F, dFdeta)

      if (k /= 1) then
        dFdeta = k * F**(k-1) * dFdeta
        F      = F**k
      end if
    end subroutine

  end subroutine

end module
