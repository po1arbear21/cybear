module example_voltage_m

  use input_m,    only: input_file
  use variable_m, only: variable

  implicit none

  private
  public init_cont_v
  public cont_v, voltage

  type, extends(variable) :: voltage
    !! electric voltage
    real, pointer :: x => null()
  contains
    procedure init => voltage_init
  end type

  type(voltage), allocatable :: cont_v(:)

contains

  subroutine voltage_init(this, f, sid)
    class(voltage),   intent(out) :: this
    type(input_file), intent(in)  :: f
    integer,          intent(in)  :: sid

    real                           :: v
    type(grid_data0_real), pointer :: p

    call this%variable_init("voltage", "V")

    ! get pointer to data
    p      => this%data%get_ptr0()
    this%x => p%data

    ! set value for voltage
    call f%get(sid, "V", v)
    call this%set(v)
  end subroutine

end module
