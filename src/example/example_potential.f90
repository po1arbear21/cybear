module example_potential_m

  use example_device_m,  only: grd
  use grid_m,            only: grid_data1_real, IDX_VERTEX
  use variable_m,        only: variable_real

  implicit none

  private
  public pot, potential

  type, extends(variable_real) :: potential
    !! electric potential
    real, pointer :: x(:) => null()
  contains
    procedure :: init => potential_init
  end type

  type(potential) :: pot

contains

  subroutine potential_init(this)
    class(potential), intent(out) :: this

    type(grid_data1_real), pointer :: p => null()

    call this%variable_init("pot", "V", g = grd, idx_type = IDX_VERTEX, idx_dir = 0)

    ! get pointer to data
    p      => this%data%get_ptr1()
    this%x => p%data
  end subroutine

end module
