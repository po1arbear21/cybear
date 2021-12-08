module voltage_m

  use grid_m,     only: grid_data0_real
  use variable_m, only: variable_real

  implicit none

  private
  public voltage

  type, extends(variable_real) :: voltage
    !! node voltage (potential)
    real, pointer :: x => null()
      !! pointer to value
  contains
    procedure :: init => voltage_init
  end type

contains

  subroutine voltage_init(this, name)
    !! initialize voltage variable
    class(voltage), intent(out) :: this
    character(*),   intent(in)  :: name
      !! variable name

    type(grid_data0_real), pointer :: p

    ! init base
    call this%variable_init(name, "V")

    ! set pointer to value
    p      => this%data%get_ptr0()
    this%x => p%data
  end subroutine

end module
