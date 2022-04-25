module vsource_m

  use circuit_m,      only: circuit, component, terminal
  use current_m,      only: current
  use jacobian_m,     only: jacobian_ptr
  use res_equation_m, only: res_equation
  use string_m,       only: new_string
  use voltage_m,      only: voltage

  implicit none

  private
  public vsource

  type, extends(res_equation) :: vsource_equ
    !! voltage source equation (V0 + V1 - V2 == 0)

    type(voltage), pointer :: volt0 => null()
      !! input voltage
    type(voltage), pointer :: volt1 => null()
      !! first terminal voltage
    type(voltage), pointer :: volt2 => null()
      !! second terminal voltage
    type(current), pointer :: curr1 => null()
      !! first terminal current (main variable)
    type(current), pointer :: curr2 => null()
      !! second terminal current (= - curr1; provided)

    type(jacobian_ptr) :: jaco_volt(0:2)
      !! jacobians of f wrt voltages
    type(jacobian_ptr) :: jaco_curr1
      !! jacobian of curr2 wrt curr1
  contains
    procedure :: init => vsource_equ_init
    procedure :: eval => vsource_equ_eval
  end type

  type vsource
    !! voltage source

    type(voltage)            :: V
      !! input voltage
    type(component), pointer :: comp => null()
      !! pointer to circuit component
    type(vsource_equ)        :: equ
      !! voltage source equation
  contains
    procedure :: init => vsource_init
  end type

contains

  subroutine vsource_equ_init(this, name, term1, term2, volt)
    !! initialize voltage source equation
    class(vsource_equ),     intent(out) :: this
    character(*),           intent(in)  :: name
      !! voltage source name
    type(terminal), target, intent(in)  :: term1
      !! first terminal
    type(terminal), target, intent(in)  :: term2
      !! second terminal
    type(voltage),  target, intent(in)  :: volt
      !! input voltage (determines difference between volt1 and volt2)

    integer :: ivolt0, ivolt1, ivolt2, icurr1, icurr2, dum(0)

    ! init base
    call this%equation_init(name)

    ! save pointers
    this%volt0 => volt
    this%volt1 => term1%nd%volt
    this%volt2 => term2%nd%volt
    this%curr1 => term1%curr
    this%curr2 => term2%curr

    ! set main variable
    call this%init_f(this%curr1)

    ! depend on voltages and curr1; provide curr2
    ivolt0 = this%depend(this%volt0)
    ivolt1 = this%depend(this%volt1)
    ivolt2 = this%depend(this%volt2)
    icurr1 = this%depend(this%curr1)
    icurr2 = this%provide(this%curr2)

    ! init voltage jacobians
    this%jaco_volt(0)%p => this%init_jaco_f(ivolt0, const = .true.)
    this%jaco_volt(1)%p => this%init_jaco_f(ivolt1, const = .true.)
    this%jaco_volt(2)%p => this%init_jaco_f(ivolt2, const = .true.)
    call this%jaco_volt(0)%p%set(dum, dum,  1.0)
    call this%jaco_volt(1)%p%set(dum, dum,  1.0)
    call this%jaco_volt(2)%p%set(dum, dum, -1.0)

    ! init current jacobian
    this%jaco_curr1%p => this%init_jaco(icurr2, icurr1, const = .true.)
    call this%jaco_curr1%p%set(dum, dum, -1.0)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine vsource_equ_eval(this)
    !! evaluate voltage source equation
    class(vsource_equ),  intent(inout) :: this

    call this%f%set(this%volt0%get() + this%volt1%get() - this%volt2%get())
    call this%curr2%set(- this%curr1%get())
  end subroutine

  subroutine vsource_init(this, crt, name, node1, node2)
    !! initialize voltage source
    class(vsource), intent(out)   :: this
    type(circuit),  intent(inout) :: crt
      !! circuit to add this to
    character(*),   intent(in)    :: name
      !! component name
    character(*),   intent(in)    :: node1
      !! name of first node this is connected to
    character(*),   intent(in)    :: node2
      !! name of first node this is connected to

    ! init voltage variable and provide as input variable
    call this%V%init(name // ".V")
    call crt%sys%provide(this%V, input = .true.)

    ! create new component
    this%comp => crt%add_component(name, [new_string(node1), new_string(node2)])

    ! initialize equation and add to system
    call this%equ%init(name, this%comp%terms(1), this%comp%terms(2), this%V)
    call crt%sys%add_equation(this%equ)
  end subroutine

end module
