module current_m

  use grid_m,     only: grid_data0_real
  use variable_m, only: variable_real

  implicit none

  private
  public current

  type, extends(variable_real) :: current
    !! current through component
    real, pointer :: x => null()
      !! pointer to value
  contains
    procedure :: init => current_init
  end type

contains

  subroutine current_init(this, name)
    !! initialize current variable
    class(current), intent(out) :: this
    character(*),   intent(in)  :: name
      !! variable name

    type(grid_data0_real), pointer :: p

    ! init base
    call this%variable_init(name, "A")

    ! set pointer to value
    p      => this%data%get_ptr0()
    this%x => p%data
  end subroutine

end module
