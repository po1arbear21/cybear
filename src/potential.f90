module potential_m

  use device_params_m, only: device_params
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data2_real
  use variable_m,      only: variable_real

  implicit none

  private
  public potential

  type, extends(variable_real) :: potential
    !! electrostatic potential
    real, pointer :: x(:,:) => null()
      !! direct pointer to data for easy access
  contains
    procedure :: init => potential_init
  end type

contains

  subroutine potential_init(this, par)
    !! initialize electrostatic potential
    class(potential),    intent(out) :: this
    type(device_params), intent(in)  :: par
      !! device parameters

    type(grid_data2_real), pointer :: p => null()

    ! init base
    call this%variable_init("pot", "V", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)

    ! get pointer to data
    p      => this%data%get_ptr2()
    this%x => p%data
  end subroutine

end module
