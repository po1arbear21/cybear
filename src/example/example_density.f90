module example_density_m

  use example_device_m,    only: grd
  use example_potential_m, only: pot
  use grid_m,              only: grid_data1_real, IDX_VERTEX
  use variable_m,          only: variable_real

  implicit none

  private
  public dens

  type, extends(variable_real) :: density
    !! density
    real, pointer :: x(:) => null()
  contains
    procedure :: init => density_init
  end type

  type(density) :: dens

contains

  subroutine density_init(this)
    class(density), intent(out) :: this

    type(grid_data1_real), pointer :: p => null()

    call this%variable_init("n", "1/cm^3", g = grd, idx_type = IDX_VERTEX, idx_dir = 0)

    ! get pointer to data
    p      => this%data%get_ptr1()
    this%x => p%data
  end subroutine

end module
