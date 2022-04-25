module resistor_m

  use circuit_m,  only: circuit, component, terminal
  use current_m,  only: current
  use equation_m, only: equation
  use jacobian_m, only: jacobian_ptr
  use string_m,   only: new_string
  use voltage_m,  only: voltage

  implicit none

  private
  public resistor

  type, extends(equation) :: resistor_equ
    !! resistor equation

    type(voltage), pointer :: volt1 => null()
      !! first terminal voltage
    type(voltage), pointer :: volt2 => null()
      !! second terminal voltage
    type(current), pointer :: curr1 => null()
      !! first terminal current
    type(current), pointer :: curr2 => null()
      !! second terminal current

    real, pointer :: R => null()
      !! resistance

    type(jacobian_ptr) :: jaco_volt(2,2)
      !! jacobians wrt voltages
  contains
    procedure :: init => resistor_equ_init
    procedure :: eval => resistor_equ_eval
  end type

  type resistor
    type(component), pointer :: comp => null()
      !! pointer to circuit component
    type(resistor_equ)       :: equ
      !! resistor equation
  contains
    procedure :: init => resistor_init
  end type

contains

  subroutine resistor_equ_init(this, name, term1, term2, R)
    !! initialize resistor equation
    class(resistor_equ),    intent(out) :: this
    character(*),           intent(in)  :: name
      !! resistor name
    type(terminal), target, intent(in)  :: term1
      !! first terminal
    type(terminal), target, intent(in)  :: term2
      !! second terminal
    real,           target, intent(in)  :: R
      !! resistance

    integer :: ivolt1, ivolt2, icurr1, icurr2, dum(0)

    ! init base
    call this%equation_init(name)

    ! save pointers to variables
    this%volt1 => term1%nd%volt
    this%volt2 => term2%nd%volt
    this%curr1 => term1%curr
    this%curr2 => term2%curr
    this%R     => R

    ! depend on voltages
    ivolt1 = this%depend(this%volt1)
    ivolt2 = this%depend(this%volt2)

    ! provide currents
    icurr1 = this%provide(this%curr1)
    icurr2 = this%provide(this%curr2)

    ! init jacobians
    this%jaco_volt(1,1)%p => this%init_jaco(icurr1, ivolt1, const = .true.)
    this%jaco_volt(1,2)%p => this%init_jaco(icurr1, ivolt2, const = .true.)
    this%jaco_volt(2,1)%p => this%init_jaco(icurr2, ivolt1, const = .true.)
    this%jaco_volt(2,2)%p => this%init_jaco(icurr2, ivolt2, const = .true.)

    ! I1 = (V1 - V2) / R
    call this%jaco_volt(1,1)%p%set(dum, dum,   1.0 / R)
    call this%jaco_volt(1,2)%p%set(dum, dum, - 1.0 / R)

    ! I2 = (V2 - V1) / R
    call this%jaco_volt(2,1)%p%set(dum, dum, - 1.0 / R)
    call this%jaco_volt(2,2)%p%set(dum, dum,   1.0 / R)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine resistor_equ_eval(this)
    !! evaluate resistor equation
    class(resistor_equ), intent(inout) :: this

    call this%curr1%set((this%volt1%get() - this%volt2%get()) / this%R)
    call this%curr2%set((this%volt2%get() - this%volt1%get()) / this%R)
  end subroutine

  subroutine resistor_init(this, crt, name, node1, node2, R)
    !! initialize resistor
    class(resistor), intent(out)   :: this
    type(circuit),   intent(inout) :: crt
      !! circuit to add this to
    character(*),    intent(in)    :: name
      !! component name
    character(*),    intent(in)    :: node1
      !! name of first node this is connected to
    character(*),    intent(in)    :: node2
      !! name of first node this is connected to
    real,            intent(in)    :: R
      !! resistance

    ! create new component
    this%comp => crt%add_component(name, [new_string(node1), new_string(node2)])

    ! initialize equation and add to system
    call this%equ%init(name, this%comp%terms(1), this%comp%terms(2), R)
    call crt%sys%add_equation(this%equ)
  end subroutine

end module
