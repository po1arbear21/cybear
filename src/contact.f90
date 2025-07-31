m4_include(util/macro.f90.inc)

module contact_m

  use ieee_arithmetic, only: ieee_is_finite
  use error_m,         only: assert_failed, program_error
  use grid_m,          only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME
  use grid_data_m,     only: grid_data_real
  use high_precision_m
  use newton_m,        only: newton1D, newton1D_opt
  use normalization_m, only: norm, denorm
  use semiconductor_m, only: CR_ELEC, CR_HOLE, CR_CHARGE, DOP_DCON, DOP_ACON, DOP_CHARGE, DOS_PARABOLIC, DIST_MAXWELL, semiconductor

  implicit none

  private
  public :: CT_OHMIC, CT_GATE, CT_SCHOTTKY, contact

  ! contact types
  integer, parameter :: CT_OHMIC   = 1
  integer, parameter :: CT_GATE    = 2
  integer, parameter :: CT_SCHOTTKY = 3

  type contact
    !! device contact

    character(:), allocatable :: name
      !! contact name
    integer                   :: type
      !! type of contact (CT_OHMIC, CT_GATE, CT_SCHOTTKY)
    real                      :: phims
      !! metal-semiconductor workfunction difference
    real                      :: barrier_height
      !! Schottky barrier height Φ_B (eV)
    real                      :: richardson_const
      !! Richardson constant A* (A/cm²/K²)
    real                      :: surf_recomb_vel(2)
      !! surface recombination velocities S_n, S_p (cm/s)
    logical                   :: tunneling_enabled
      !! enable field emission tunneling
  contains
    procedure :: set_phims_ohmic => contact_set_phims_ohmic
  end type

contains

  subroutine contact_set_phims_ohmic(this, ci0, ci1, dop, ii_E_dop, smc)
    !! set phims for ohmic contacts by assuming charge neutrality
    class(contact),        intent(inout) :: this
    integer,               intent(in)    :: ci0, ci1
      !! lower/upper carrier index (CR_ELEC, CR_HOLE)
    real,                  intent(in)    :: dop(2)
      !! donor/acceptor concentration
    real,                  intent(in)    :: ii_E_dop(2)
      !! dopant energy level (for incomplete ionization)
    type(semiconductor),   intent(in)    :: smc
      !! semiconductor parameters

    real               :: dum(0), fmin, fmax, xmin, xmax, x0
    type(newton1D_opt) :: opt

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
    call newton1D(phims_newton, dum, opt, x0, this%phims)

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
      do ci = ci0, ci1
        ch = CR_CHARGE(ci)

        if ((smc%dos == DOS_PARABOLIC) .and. (smc%dist == DIST_MAXWELL)) then
          hrho = hrho + ch * sqrt(smc%edos(1) * smc%edos(2)) * exp(- ch * phims - 0.5 * smc%band_gap)
          drho = drho -      sqrt(smc%edos(1) * smc%edos(2)) * exp(- ch * phims - 0.5 * smc%band_gap)
        else
          call smc%get_dist(ch * (smc%band_edge(ci) - phims), 0, F, dF)
          hrho = hrho + ch * smc%edos(ci) *  F
          drho = drho -      smc%edos(ci) * dF
        end if
      end do

      ! doping
      do ci = DOP_DCON, DOP_ACON
        ch = CR_CHARGE(ci)
        if (smc%incomp_ion .and. (dop(ci) < smc%ii_dop_th(ci))) then
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
      end do

      rho = hp_to_real(hrho)
      if (present(drhodphims)) drhodphims = drho
      if (present(drhodp)) then
        m4_ignore(drhodp)
      end if
    end subroutine
  end subroutine

end module
