m4_include(util/macro.f90.inc)

module contact_m

  use ieee_arithmetic, only: ieee_is_finite
  use distributions_m, only: fermi_dirac_integral_1h, inv_fermi_dirac_integral_1h
  use error_m,         only: assert_failed, program_error
  use grid_m,          only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME
  use grid_data_m,     only: grid_data_real
  use high_precision_m
  use newton_m,        only: newton1D, newton1D_opt
  use normalization_m, only: norm, denorm
  use semiconductor_m, only: CR_ELEC, CR_HOLE, CR_CHARGE, DOP_DCON, DOP_ACON, DOP_CHARGE, semiconductor

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

  subroutine contact_set_phims_ohmic(this, ci0, ci1, dop, asb, edop, smc)
    !! set phims for ohmic contacts by assuming charge neutrality
    class(contact),        intent(inout) :: this
    integer,               intent(in)    :: ci0, ci1
      !! lower/upper carrier index (CR_ELEC, CR_HOLE)
    real,                  intent(in)    :: dop(2)
      !! donor/acceptor concentration
    real,                  intent(in)    :: asb(2)
      !! Altermatt-Schenk b parameter
    real,                  intent(in)    :: edop(2)
      !! Altermatt-Schenk doping energy level
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

    subroutine phims_newton(x, p, f, dfdx, dfdp)
      !! residual for newton iteration: f = rho = p - n + ND^{+} - NA^{-}
      real,              intent(in)  :: x
        !! argument (phims)
      real,              intent(in)  :: p(:)
        !! parameters (empty)
      real,              intent(out) :: f
        !! output function value
      real,    optional, intent(out) :: dfdx
        !! optional output derivative of f wrt x
      real,    optional, intent(out) :: dfdp(:)
        !! optional output derivatives of f wrt p

      integer       :: ci
      real          :: ch, dion, F1h, dF1h, dff, t
      type(hp_real) :: ff, e, fac, ion

      m4_ignore(p)

      ! work with high precision residual, double precision derivative
      ff  = real_to_hp(0.0)
      dff = 0.0

      ! densities
      do ci = ci0, ci1
        ch = CR_CHARGE(ci)

        if (smc%degen) then
          call fermi_dirac_integral_1h(ch * (smc%band_edge(ci) - x), F1h, dF1h)
          ff  =  ff + ch * smc%edos(ci) * F1h
          dff = dff -      smc%edos(ci) * dF1h
        else
          ff  =  ff + ch * sqrt(smc%edos(1) * smc%edos(2)) * exp(- ch * x - 0.5 * smc%band_gap)
          dff = dff -      sqrt(smc%edos(1) * smc%edos(2)) * exp(- ch * x - 0.5 * smc%band_gap)
        end if
      end do

      ! doping
      do ci = DOP_DCON, DOP_ACON
        ch = CR_CHARGE(ci)
        if (smc%incomp_ion) then
          e = exp(ch * TwoSum(x, - (smc%band_edge(ci) + ch * edop(ci))))
          t = hp_to_real(e)
          if (.not. ieee_is_finite(t) .or. (t > 1e300)) then
            fac = real_to_hp(0.0)
          else
            fac  = 1.0 / (1.0 + smc%g_dop(ci) * e)
          end if
          ion  = 1.0 - asb(ci) * fac
          dion = hp_to_real(asb(ci) * fac * (1.0 - fac))

          ff  = ff   - ch * dop(ci) * ion
          dff = dff -      dop(ci) * dion
        else
          ff = ff - ch * dop(ci)
        end if
      end do

      f = hp_to_real(ff)
      if (present(dfdx)) dfdx = dff
      if (present(dfdp)) then
        m4_ignore(dfdp)
      end if
    end subroutine
  end subroutine

end module
