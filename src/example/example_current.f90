module example_current_m

  use grid_m,                   only: grid_data0_real
  use variable_m,               only: variable

  implicit none

  private
  public current

  type, extends(variable) :: current
    !! electric current
    real, pointer :: x => null()
  contains
    procedure :: init => current_init
  end type

contains

  subroutine current_init(this, cont_name)
    class(current),   intent(out) :: this
    character(*),     intent(in)  :: cont_name

    type(grid_data0_real), pointer :: p

    call this%variable_init("current_"//cont_name, "A/um^2")

    ! get pointer to data
    p      => this%data%get_ptr0()
    this%x => p%data
  end subroutine

end module
