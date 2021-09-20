module example_voltage_m

  use grid_m,     only: grid_data0_real
  use variable_m, only: variable

  implicit none

  private
  public voltage
  !public volt, voltage

  type, extends(variable) :: voltage
    !! electric voltage
    real, pointer :: x => null()
  contains
    procedure :: init => voltage_init
  end type

  type(voltage), allocatable :: volt(:)

contains

  subroutine voltage_init(this, cont_name)
    class(voltage),   intent(out) :: this
    character(*),     intent(in)  :: cont_name

    type(grid_data0_real), pointer :: p

    call this%variable_init("voltage_"//cont_name, "V")

    ! get pointer to data
    p      => this%data%get_ptr0()
    this%x => p%data
  end subroutine

end module
