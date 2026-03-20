m4_include(util/macro.f90.inc)

module semiconductor_m

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_positive_inf, ieee_negative_inf

  use lookup_table_m, only: par_int_table, def_int_table
  use error_m,        only: assert_failed, program_error
  use fukushima_m,    only: fermi_dirac_integral, ifd1h
  use math_m,         only: expm1
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
  character(*), parameter :: RATE_NAME(2)  = [ "DC", "AV" ]

  integer,      parameter :: DOS_PARABOLIC      = 1
  integer,      parameter :: DOS_PARABOLIC_TAIL = 2

  integer,      parameter :: DIST_MAXWELL = 1
    !! Maxwell-Boltzmann distribution density
  integer,      parameter :: DIST_FERMI   = 2
    !! Fermi-Dirac distribution density

  integer,      parameter :: STAB_SG    = 1
    !! Scharfetter-Gummel (not thermodynamically consistent)
  integer,      parameter :: STAB_ED    = 2
    !! enhanced diffusion (thermodynamically consistent)
  integer,      parameter :: STAB_EXACT = 3
    !! exact (numerical) solution of integral equation (thermodynamically consistent)

  type semiconductor
    integer :: ci0, ci1
      !! enabled carrier index range (maximal: CR_ELEC..CR_HOLE)

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
    logical                  :: reg
      !! use regularization in DOS
    real                     :: reg_params(2)
      !! regularization parameters
    integer                  :: dist
      !! distribution density (DIST_MAXWELL or DIST_FERMI)
    integer                  :: stab
      !! stabilization method for degenerate case (STAB_SG, STAB_ED, STAB_ED2, STAB_EXACT)
    type(par_int_table)      :: Ftab(0:4)
      !! distribution function + derivatives
    type(def_int_table)      :: Itab(-5:5)
      !! integral_0^{eta} F(eta')^k deta'; k = -5 to 5

    !! velocity saturation parameters (Caughey-Thomas model)
    logical           :: mob_sat
      !! enable/disable mobility saturation
    real, allocatable :: beta(:)
      !! Caughey-Thomas beta parameter
    real, allocatable :: v_sat(:)
      !! Caughey-Thomas saturation velocity

    !! incomplete ionization parameters
    logical           :: incomp_ion
      !! enable/disable incomplete ionization (Pearson-Bardeen model)
    real, allocatable :: ii_tau(:)
      !! generation/recombination time constant (carrier index)
    real, allocatable :: ii_E_dop0(:)
      !! dopant energy relative to the carrier band for a single dopant
    real, allocatable :: ii_g(:)
      !! degeneracy factor of dopant states
    real, allocatable :: ii_N_crit(:)
      !! critical concentration (alternative to altermatt-schenk)
    real, allocatable :: ii_dop_th(:)
      !! full ionization if doping > threshold
    logical           :: ii_tun
      !! enable/disable tunneling ionization
    real, allocatable :: ii_tau_tun(:)
      !! tunneling time constant
    real, allocatable :: ii_m_tun(:)
      !! tunneling effective mass
    logical           :: ii_pf
      !! enable/disable poole-frenkel effect
    real              :: ii_pf_a
      !! smoothing parameter for softmin function in PF barrier lowering
    real              :: ii_ef_min
      !! minimum electric field to ensure differentiability: F = sqrt(Fx^2 + Fy^2 + Fz^2 + ii_ef_min^2)

  contains
    procedure :: init_dist    => semiconductor_init_dist
    procedure :: get_dist     => semiconductor_get_dist
    procedure :: get_inv_dist => semiconductor_get_inv_dist
    procedure :: get_int_dist => semiconductor_get_int_dist
  end type

contains

  subroutine semiconductor_init_dist(this, tabledir)
    class(semiconductor), intent(inout) :: this
    character(*),         intent(in)    :: tabledir
      !! file system table directory

    real, parameter :: ETA0 = -1e4, ETA1 = 5e3

    character(255)            :: project_root
    character(:), allocatable :: fname, fdir
    integer                   :: k, i, j
    logical                   :: exist
    real                      :: inf, p(1)

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

    inf = ieee_value(inf, IEEE_POSITIVE_INF)
    if ((this%dist == DIST_FERMI) .and. (this%dos /= DOS_PARABOLIC)) then
      print "(A)", "generate Ftab for distribution function"
      do k = lbound(this%Ftab, 1), ubound(this%Ftab, 1)
        fname = "dist_" // int2str(this%dos) // int2str(this%dist) // "_" // int2str(k)
        call this%Ftab(k)%init(fdir, fname, ETA0, ETA1, Z_dist, [-inf, inf], 1e-12, 1e-15, 1e-6)
      end do
    end if

    if ((this%dist /= DIST_MAXWELL) .or. this%reg) then
      print "(A)", "generate Itab for generalized Scharfetter-Gummel"
      do k = -5, 5
        if (k == 0) cycle
        p(1) = k
        call this%Itab(k)%init(ETA0, ETA1, 3, Fk, p, 5e-15, 0.0, 1e-10)
      end do
    end if

  contains

    subroutine Z_dist(t, p, Zf, dZfdt, dZfdp)
      use math_m, only: PI

      real, intent(in)  :: t
      real, intent(in)  :: p(:)
      real, intent(out) :: Zf
      real, intent(out) :: dZfdt
      real, intent(out) :: dZfdp(:)

      real :: t0, Z, eta, f(0:1)

      ! density of states
      m4_assert(this%dos == DOS_PARABOLIC_TAIL)
      t0 = this%dos_params(1)
      if (t > t0) then
        Z = 2 * sqrt(t / PI)
      else
        Z = 2 * sqrt(t0 / PI) * exp(0.5 * t / t0) * exp(-0.5)
      end if

      ! distribution density
      eta = p(1)
      call fermi_dirac(t - eta, k, f)
      if (mod(k, 2) == 0) then
        f(1) = - f(1)
      else
        f(0) = - f(0)
      end if

      ! combine
      Zf       = Z * f(0)
      dZfdt    = 0
      dZfdp(1) = Z * f(1)
    end subroutine

    subroutine Fk(eta, p, F, dFdeta, dFdp)
      real, intent(in)  :: eta
      real, intent(in)  :: p(:)
      real, intent(out) :: F
      real, intent(out) :: dFdeta
      real, intent(out) :: dFdp(:)

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
      real,    intent(in)  :: u
        !! argument u = t - eta
      integer, intent(in)  :: k
        !! derivative index
      real,    intent(out) :: f(0:1)
        !! output f^{(k)} and f^{(k+1)}

      real, parameter :: uc1 = acosh(2.0), uc2 = acosh(5.0)

      integer :: j
      real    :: c, s1, s2

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
          s2 = cosh(2 * u)
          c  = cosh(0.5 * u)
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
      !! derivative index (k = 0...4)
    real,                 intent(out) :: F
      !! output cumulative distribution function
    real,                 intent(out) :: dFdeta
      !! output derivative of F wrt eta

    integer :: j
    real    :: A, B, t, c, ff(0:1)

    if (this%dist == DIST_MAXWELL) then
      ! Maxwell-Boltzmann
      F      = exp(eta)
      dFdeta = F
    elseif (this%dist == DIST_FERMI) then
      if (this%dos == DOS_PARABOLIC) then
        ! Fermi-Dirac integral
        F      = fermi_dirac_integral(eta, k)
        dFdeta = fermi_dirac_integral(eta, k + 1)
      else
        ! use table for more complicated cases
        call this%Ftab(k)%get(eta, F, dFdeta)
      end if
    end if

    ! regularization
    if (this%reg) then
      A = this%reg_params(1)
      B = this%reg_params(2)

      t = tanh(0.5 * B * eta)
      c = 1 / cosh(0.5 * B * eta)**2

      do j = 0, 1
        select case (k + j)
        case (0)
          ff(j) =   2 * exp(B * eta) / (exp(B * eta) + 1)
        case (1)
          ff(j) =   0.5 * c
        case (2)
          ff(j) = - 0.5 * c * t
        case (3)
          ff(j) =   0.5 * c * (t**2 - 0.5 * c)
        case (4)
          ff(j) = - 0.5 * c * t * (t**2 - 2 * c)
        case (5)
          ff(j) =   0.5 * c * (t**4 + c**2 - 5.5 * t**2 * c)
        end select
        ff(j) = ff(j) * 0.5 * A * B**(k + j)
      end do
      F      = F + ff(0)
      dFdeta = dFdeta + ff(1)
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

    integer, parameter :: MAX_IT = 15
    real,    parameter :: RTOL = 5e-14, ATOL = 5e-16
    real,    parameter :: G = gamma(1.5)

    integer :: it
    real    :: A, B, d, FF, dFFdeta, err

    ! safety check
    if (F <= 0) then
      print "(A,ES41.32E3)", "F = ", F
      call program_error("F must be positive")
    end if

    A = this%reg_params(1)
    B = this%reg_params(2)

    if (this%reg .and. (F <= 0.4 * A)) then
      eta = log(2 * F / A / (2 - 2 * F / A)) / B
    else
      if (this%dist == DIST_MAXWELL) then
        eta = log(F)
        detadF = 1 / F
      elseif (this%dos == DOS_PARABOLIC) then
        call ifd1h(G * F, eta, detadF)
        detadF = detadF * G
      else
        call this%Ftab(0)%inv(F, eta, detadF)
      end if
    end if

    ! Newton iteration if regularization is active
    if (this%reg) then
      it = 0
      err = huge(1.0)
      do while (err > max(RTOL * abs(eta), ATOL))
        it = it + 1
        if (it > MAX_IT) then
          print "(A,ES25.16E3)", "F = ", F
          call program_error("No convergence after "//int2str(MAX_IT)//" iterations")
        end if

        call this%get_dist(eta, 0, FF, dFFdeta)

        d = - (FF - F) / dFFdeta

        err = abs(d)
        eta = eta + d
      end do

      detadF = 1 / dFFdeta
    end if

    ! second safety check
    if (.not. (ieee_is_finite(eta) .and. ieee_is_finite(detadF)) .or. (detadF <= 0)) then
      print "(A,ES25.16E3)", "F      = ", F
      print "(A,ES25.16E3)", "eta    = ", eta
      print "(A,ES25.16E3)", "detadF = ", detadF
      call program_error("not finite or detadF non-positive")
    end if
  end subroutine

  subroutine semiconductor_get_int_dist(this, eta, k, rtol, I, dIdeta)
    !! integral_eta(1)^eta(2) dist(eta)^k deta
    class(semiconductor), intent(in)  :: this
    real,                 intent(in)  :: eta(2)
      !! integration bounds
    integer,              intent(in)  :: k
      !! exponent
    real,                 intent(in)  :: rtol
      !! relative error tolerance
    real,                 intent(out) :: I
      !! output integral over dist^k
    real,                 intent(out) :: dIdeta(2)
      !! output derivatives of I wrt eta

    if (k == 0) then
      I      = eta(2) - eta(1)
      dIdeta = [-1, 1]
    else
      if ((this%dist == DIST_MAXWELL) .and. (.not. this%reg)) then
        I = exp(k*eta(1)) * expm1(k*(eta(2) - eta(1))) / k
        dIdeta(1) = - exp(k*eta(1))
        dIdeta(2) =   exp(k*eta(2))
      else
        call this%Itab(k)%get(Fk, eta(1), eta(2), rtol, I, dIdeta(1), dIdeta(2))
      end if
    end if

  contains

    subroutine Fk(eta, p, F, dFdeta, dFdp)
      real, intent(in)  :: eta
      real, intent(in)  :: p(:)
      real, intent(out) :: F
      real, intent(out) :: dFdeta
      real, intent(out) :: dFdp(:)

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
