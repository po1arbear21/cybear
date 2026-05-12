m4_include(util/macro.f90.inc)

module contact_m

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  use error_m,         only: assert_failed, program_error
  use grid_m,          only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME
  use grid_data_m,     only: grid_data_real
  use high_precision_m
  use newton_m,        only: newton, newton_opt
  use normalization_m, only: norm, denorm
  use semiconductor_m, only: CR_ELEC, CR_HOLE, CR_CHARGE, DOP_DCON, DOP_ACON, semiconductor

  implicit none

  private
  public :: CT_OHMIC, CT_GATE, CT_SCHOTTKY, CT_REALOHMIC, contact

  ! contact types
  integer, parameter :: CT_OHMIC     = 1
  integer, parameter :: CT_GATE      = 2
  integer, parameter :: CT_SCHOTTKY  = 3
  integer, parameter :: CT_REALOHMIC = 4

  type contact
    !! device contact

    character(:), allocatable :: name
      !! contact name
    integer                   :: type
      !! type of contact (CT_OHMIC, CT_GATE, CT_SCHOTTKY, CT_REALOHMIC)
    real                      :: phims
      !! metal-semiconductor workfunction difference
    real                      :: phi_b
      !! Schottky barrier height for ELECTRONS (phi_Bn = E_C - E_F_metal, normalized to kT).
      !! Hole barrier is phi_Bp = band_gap - phi_b. Must satisfy 0 <= phi_b <= band_gap.
    real                      :: A_richardson_n
      !! Richardson constant for electrons [A/cm^2/K^2]
    real                      :: A_richardson_p
      !! Richardson constant for holes [A/cm^2/K^2]
    logical                   :: ifbl = .false.
      !! image force barrier lowering flag
    logical                   :: tunneling = .false.
      !! enable Tsu-Esaki tunneling model
    real                      :: m_tunnel_n = 1.0
      !! tunneling effective mass ratio for electrons (m*/m0)
    real                      :: m_tunnel_p = 1.0
      !! tunneling effective mass ratio for holes (m*/m0)
    real                      :: vrec
      !! metal-semiconductor recombination velocity (real ohmic contact)
  contains
    procedure :: set_phims_ohmic    => contact_set_phims_ohmic
    procedure :: set_phims_schottky => contact_set_phims_schottky
  end type

contains

  subroutine contact_set_phims_ohmic(this, smc, dop, ii_E_dop)
    !! set phims for ohmic contacts by assuming charge neutrality
    class(contact),      intent(inout) :: this
    type(semiconductor), intent(in)    :: smc
      !! semiconductor parameters
    real,                intent(in)    :: dop(2)
      !! donor/acceptor concentration
    real, optional,      intent(in)    :: ii_E_dop(2)
      !! dopant energy level (for incomplete ionization)

    real             :: dum(0), fmin, fmax, xmin, xmax, x0
    type(newton_opt) :: opt

    m4_assert(present(ii_E_dop) .eqv. smc%incomp_ion)

    ! coarse search for bounds
    if (dop(1) > dop(2)) then
      xmin = 0
      xmax = 10
      call phims_newton(xmin, dum, fmin)
      call phims_newton(xmax, dum, fmax)

      do while (sign(1.0, fmin) == sign(1.0, fmax))
        xmax = xmax + 10.0
        call phims_newton(xmax, dum, fmax)
      end do
      x0 = xmax
    else
      xmin = -10
      xmax = 0
      call phims_newton(xmin, dum, fmin)
      call phims_newton(xmax, dum, fmax)
      do while (sign(1.0, fmin) == sign(1.0, fmax))
        xmin = xmin - 10.0
        call phims_newton(xmin, dum, fmin)
      end do
      x0 = xmin
    end if

    call opt%init(xmin = xmin, xmax = xmax, dx_lim = 10.0)
    call newton(phims_newton, dum, opt, x0, this%phims)

    ! safety check
    m4_assert(ieee_is_finite(this%phims))

  contains

    subroutine phims_newton(phims, p, rho, drhodphims, drhodp)
      !! residual for newton iteration: rho = rho = p - n + ND^{+} - NA^{-}
      real,              intent(in)  :: phims
        !! argument (phims)
      real,              intent(in)  :: p(:)
        !! parameters (empty)
      real,              intent(out) :: rho
        !! output function value
      real,    optional, intent(out) :: drhodphims
        !! optional output derivative of rho wrt phims
      real,    optional, intent(out) :: drhodp(:)
        !! optional output derivatives of rho wrt p

      integer       :: ci
      real          :: ch, dion, F, dF, drho, t
      type(hp_real) :: hrho, e, fac, ion

      m4_ignore(p)

      ! work with high precision residual, double precision derivative
      hrho  = real_to_hp(0.0)
      drho = 0.0

      ! densities
      do ci = smc%ci0, smc%ci1
        ch = CR_CHARGE(ci)
        call smc%get_dist(ch * (smc%band_edge(ci) - phims), 0, F, dF)
        hrho = hrho + ch * smc%edos(ci) *  F
        drho = drho -      smc%edos(ci) * dF
      end do

      ! doping
      do ci = DOP_DCON, DOP_ACON
        ch = CR_CHARGE(ci)
        if (smc%incomp_ion) then
          if (dop(ci) < smc%ii_dop_th(ci)) then
            e = TwoSum(ii_E_dop(ci) + ch*smc%band_edge(ci), -ch*phims)
            t = hp_to_real(e)
            if (t > 300) then
              fac = real_to_hp(0.0)
            else
              e   = exp(e)
              fac = 1.0 / (1.0 + smc%ii_g(ci) * e)
            end if
            ion  = fac
            dion = hp_to_real(fac * (1.0 - fac))

            hrho = hrho - ch * dop(ci) * ion
            drho = drho -      dop(ci) * dion
          else
            ! fully ionized
            hrho = hrho - ch * dop(ci)
          end if
        else
          ! fully ionized
          hrho = hrho - ch * dop(ci)
        end if
      end do

      rho = hp_to_real(hrho)
      if (present(drhodphims)) drhodphims = drho
      if (present(drhodp)) then
        m4_ignore(drhodp)
      end if
    end subroutine
  end subroutine

  subroutine contact_set_phims_schottky(this, smc)
    !! Set phims for Schottky contacts using Sentaurus Device approach:
    !! phims = -phi_b + ln(N_c / n_i)
    !! where n_i = sqrt(N_c * N_v) * exp(-E_g / 2kT)
    class(contact),      intent(inout) :: this
    type(semiconductor), intent(in)    :: smc

    real :: ni_term

    ! validate input convention: phi_b is the electron barrier and must lie within the gap
    if (this%phi_b < 0.0 .or. this%phi_b > smc%band_gap) then
      call program_error("contact "//this%name//": phi_b (electron Schottky barrier) "// &
        & "must lie in [0, band_gap]; check the .ini convention (phi_b = E_C - E_F_metal)")
    end if

    ! ln(N_c/n_i) = 0.5*ln(N_c/N_v) + E_g/2  (in normalized units)
    ni_term = 0.5 * log(smc%edos(CR_ELEC) / smc%edos(CR_HOLE)) + 0.5 * smc%band_gap

    this%phims = -this%phi_b + ni_term

    m4_assert(ieee_is_finite(this%phims))
  end subroutine

end module
