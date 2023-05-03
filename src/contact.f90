module contact_m

  use semiconductor_m, only: CR_ELEC, CR_HOLE, semiconductor
  use newton_m,        only: newton1D, newton1D_opt
  use grid_data_m,     only: grid_data_real
  use grid_m,          only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, IDX_NAME
  use distributions_m, only: fermi_dirac_integral_1h, inv_fermi_dirac_integral_1h

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

  subroutine contact_set_phims_ohmic(this, ci0, ci1, dcon, acon, smc)
    class(contact),        intent(inout) :: this
    integer,               intent(in)    :: ci0, ci1
    real,                  intent(in)    :: dcon, acon
    type(semiconductor),   intent(in)    :: smc

    real               :: IF12, dIF12, phims0
    type(newton1D_opt) :: opt

    if ((ci0 == CR_ELEC) .and. (ci1 == CR_HOLE)) then
      if (smc%degen) then
        if (dcon > acon) then
          call inv_fermi_dirac_integral_1h(dcon / smc%edos(CR_ELEC), IF12, dIF12)
          phims0 = smc%band_edge(CR_ELEC) + IF12
        else
          call inv_fermi_dirac_integral_1h(acon / smc%edos(CR_HOLE), IF12, dIF12)
          phims0 = smc%band_edge(CR_HOLE) - IF12
        end if
        call opt%init()
        call newton1D(phims_newton, [dcon, acon], opt, phims0, this%phims)
      else
        this%phims = asinh(0.5 * (dcon - acon) / smc%n_intrin)
      end if
    elseif (ci0 == CR_ELEC) then
      if (smc%degen) then
        call inv_fermi_dirac_integral_1h(dcon / smc%edos(CR_ELEC), IF12, dIF12)
        this%phims = smc%band_edge(CR_ELEC) + IF12
      else
        this%phims = log(dcon / smc%n_intrin)
      end if
    elseif (ci1 == CR_HOLE) then
      if (smc%degen) then
        call inv_fermi_dirac_integral_1h(acon / smc%edos(CR_HOLE), IF12, dIF12)
        this%phims = smc%band_edge(CR_HOLE) - IF12
      else
        this%phims = -log(acon / smc%n_intrin)
      end if
    end if

  contains

    subroutine phims_newton(x, p, f, dfdx, dfdp)
      real,              intent(in)  :: x
        !! argument (phims)
      real,              intent(in)  :: p(:)
        !! parameters (dcon, acon)
      real,              intent(out) :: f
        !! output function value
      real,    optional, intent(out) :: dfdx
        !! optional output derivative of f wrt x
      real,    optional, intent(out) :: dfdp(:)
        !! optional output derivatives of f wrt p

      real :: Fn, dFn, Fp, dFp

      call fermi_dirac_integral_1h(x - smc%band_edge(CR_ELEC), Fn, dFn)
      call fermi_dirac_integral_1h(smc%band_edge(CR_HOLE) - x, Fp, dFp)

      f = smc%edos(CR_ELEC) * Fn - smc%edos(CR_HOLE) * Fp - p(1) + p(2)
      if (present(dfdx)) then
        dfdx = smc%edos(CR_ELEC) * dFn + smc%edos(CR_HOLE) * dFp
      end if
      if (present(dfdp)) then
        dfdp(1) = -1
        dfdp(2) = +1
      end if
    end subroutine
  end subroutine

end module
