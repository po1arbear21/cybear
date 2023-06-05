m4_include(util/macro.f90.inc)

module contact_m

  use ieee_arithmetic, only: ieee_is_finite
  use distributions_m, only: fermi_dirac_integral_1h, inv_fermi_dirac_integral_1h
  use error_m,         only: assert_failed, program_error
  use grid_m,          only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME
  use grid_data_m,     only: grid_data_real
  use newton_m,        only: newton1D, newton1D_opt
  use semiconductor_m, only: CR_ELEC, CR_HOLE, DOP_DCON, DOP_ACON, DOP_CHARGE, semiconductor

  implicit none

  private
  public :: CT_OHMIC, CT_GATE, contact

  ! contact types
  integer, parameter :: CT_OHMIC = 1
  integer, parameter :: CT_GATE  = 2

  type contact
    !! device contact

    character(:), allocatable :: name
      !! contact name
    integer                   :: type
      !! type of contact (CT_OHMIC, CT_GATE)
    real                      :: phims
      !! metal-semiconductor workfunction difference
  contains
    procedure :: set_phims_ohmic => contact_set_phims_ohmic
  end type

contains

  subroutine contact_set_phims_ohmic(this, ci0, ci1, dop, smc)
    !! set phims for ohmic contacts by assuming charge neutrality
    class(contact),        intent(inout) :: this
    integer,               intent(in)    :: ci0, ci1
      !! lower/upper carrier index (CR_ELEC, CR_HOLE)
    real,                  intent(in)    :: dop(2)
      !! donor/acceptor concentration
    type(semiconductor),   intent(in)    :: smc
      !! semiconductor parameters

    real               :: IF12, dIF12, phims0, dop_eff, L
    type(newton1D_opt) :: opt

    ! effective doping (dcon - acon)
    dop_eff = dot_product(DOP_CHARGE, dop)

    if ((ci0 == CR_ELEC) .and. (ci1 == CR_HOLE)) then
      if (smc%degen) then
        ! initial guess
        if (dop(DOP_DCON) > dop(DOP_ACON)) then
          call inv_fermi_dirac_integral_1h(dop(DOP_DCON) / smc%edos(CR_ELEC), IF12, dIF12)
          phims0 = smc%band_edge(CR_ELEC) + IF12
        else
          call inv_fermi_dirac_integral_1h(dop(DOP_ACON) / smc%edos(CR_HOLE), IF12, dIF12)
          phims0 = smc%band_edge(CR_HOLE) - IF12
        end if

        ! newton iteration
        call opt%init()
        call newton1D(phims_newton, [dop_eff], opt, phims0, this%phims)
      else
        ! analytic solution
        if (dop_eff == 0) then
          this%phims = 0
        else
          L = 0.5 * smc%band_gap + log(abs(dop_eff) / sqrt(smc%edos(1) * smc%edos(2)))
          if (L < 9) then
            this%phims = asinh(0.5 * exp(L))
          else
            ! avoid overflow, approximation is exact up to machine precision for L >= 9
            this%phims = L + exp(-2 * L)
          end if
          this%phims = sign(this%phims, dop_eff)
        end if
      end if
    elseif (ci0 == CR_ELEC) then
      m4_assert(dop_eff > 0)

      if (smc%degen) then
        call inv_fermi_dirac_integral_1h(dop_eff / smc%edos(CR_ELEC), IF12, dIF12)
        this%phims = smc%band_edge(CR_ELEC) + IF12
      else
        this%phims = 0.5 * smc%band_gap + log(dop_eff / sqrt(smc%edos(1) * smc%edos(2)))
      end if
    elseif (ci1 == CR_HOLE) then
      m4_assert(dop_eff < 0)

      if (smc%degen) then
        call inv_fermi_dirac_integral_1h(dop_eff / smc%edos(CR_HOLE), IF12, dIF12)
        this%phims = smc%band_edge(CR_HOLE) - IF12
      else
        this%phims = -0.5 * smc%band_gap - log(- dop_eff / sqrt(smc%edos(1) * smc%edos(2)))
      end if
    end if

    ! safety check
    m4_assert(ieee_is_finite(this%phims))

  contains

    subroutine phims_newton(x, p, f, dfdx, dfdp)
      !! residual for newton iteration: f = n - p - ND + NA
      real,              intent(in)  :: x
        !! argument (phims)
      real,              intent(in)  :: p(:)
        !! parameters (dcon - acon)
      real,              intent(out) :: f
        !! output function value
      real,    optional, intent(out) :: dfdx
        !! optional output derivative of f wrt x
      real,    optional, intent(out) :: dfdp(:)
        !! optional output derivatives of f wrt p

      real :: Fn, dFn, Fp, dFp

      call fermi_dirac_integral_1h(x - smc%band_edge(CR_ELEC), Fn, dFn)
      call fermi_dirac_integral_1h(smc%band_edge(CR_HOLE) - x, Fp, dFp)

      ! f = NC * F12(eta_n) - NV * F12(eta_p) - ND + NA
      f = smc%edos(CR_ELEC) * Fn - smc%edos(CR_HOLE) * Fp - p(1)
      if (present(dfdx)) dfdx = smc%edos(CR_ELEC) * dFn + smc%edos(CR_HOLE) * dFp
      if (present(dfdp)) dfdp(1) = -1
    end subroutine
  end subroutine

end module
