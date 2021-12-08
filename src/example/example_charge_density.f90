#include "../util/macro.f90.inc"

module example_charge_density_m

  use equation_m,        only: equation
  use example_density_m, only: dens
  use example_device_m,  only: dop_v, grd
  use grid_m,            only: grid_data1_real, IDX_VERTEX
  use jacobian_m,        only: jacobian, jacobian_ptr
  use stencil_m,         only: dirichlet_stencil
  use variable_m,        only: variable_real

  implicit none

  private
  public calc_charge_dens, charge_dens

  type, extends(variable_real) :: charge_density
    !! charge density
    real, pointer :: x(:) => null()
  contains
    procedure :: init => charge_density_init
  end type

  type, extends(equation) :: calc_charge_density
  contains
    procedure :: init => calc_charge_density_init
    procedure :: eval => calc_charge_density_eval
  end type

  type(charge_density)      :: charge_dens
  type(calc_charge_density) :: calc_charge_dens

contains

  subroutine charge_density_init(this)
    class(charge_density), intent(out) :: this

    type(grid_data1_real), pointer :: p => null()

    call this%variable_init("rho", "C/cm^3", g = grd, idx_type = IDX_VERTEX, idx_dir = 0)

    ! get pointer to data
    p      => this%data%get_ptr1()
    this%x => p%data
  end subroutine

  subroutine calc_charge_density_init(this)
    class(calc_charge_density), intent(out) :: this

    integer                 :: i_dep, i_prov, i
    type(jacobian), pointer :: jaco

    ! init equation
    call this%equation_init("charge_density_calc")

    ! provides charge_density and depends on density
    i_prov = this%provide(charge_dens)
    i_dep  = this%depend(dens)

    ! init jaco
    jaco => this%init_jaco(i_prov, i_dep, const = .true.)

    ! set jacobian entries
    do i = 1, size(grd%x)
      call jaco%set([i], [i], -1.0)
    end do
    call this%init_final()
  end subroutine

  subroutine calc_charge_density_eval(this)
    class(calc_charge_density), intent(inout) :: this

    IGNORE(this)

    ! calculating the charge density
    call charge_dens%set([dop_v%get() - dens%get()])
  end subroutine

end module
