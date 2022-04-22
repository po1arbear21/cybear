module capacitor_m

  use circuit_m,      only: circuit, component, terminal
  use current_m,      only: current
  use jacobian_m,     only: jacobian_ptr
  use res_equation_m, only: res_equation
  use string_m,       only: new_string
  use voltage_m,      only: voltage

  implicit none

  private
  public capacitor

  type, extends(res_equation) :: capacitor_equ
    !! capacitor equation

    type(voltage), pointer :: volt1 => null()
      !! first terminal voltage
    type(voltage), pointer :: volt2 => null()
      !! second terminal voltage
    type(current), pointer :: curr1 => null()
      !! first terminal current
    type(current), pointer :: curr2 => null()
      !! second terminal current

    type(jacobian_ptr) :: jaco_volt(2)
    type(jacobian_ptr) :: jaco_curr(2)
  contains
    procedure :: init => capacitor_equ_init
    procedure :: eval => capacitor_equ_eval
  end type

  type capacitor
    type(component), pointer :: comp => null()
      !! pointer to circuit component
    type(capacitor_equ)      :: equ
      !! capacitor equation
  contains
    procedure :: init     => capacitor_init
    procedure :: destruct => capacitor_destruct
  end type

contains

  subroutine capacitor_equ_init(this, name, term1, term2, C)
    !! initialize capacitor equation
    class(capacitor_equ),   intent(out) :: this
    character(*),           intent(in)  :: name
      !! capacitor name
    type(terminal), target, intent(in)  :: term1
      !! first terminal
    type(terminal), target, intent(in)  :: term2
      !! second terminal
    real,                   intent(in)  :: C
      !! capacitance

    integer :: ivolt1, ivolt2, icurr1, icurr2, dum(0)

    ! init base
    call this%equation_init(name)

    ! save pointers
    this%volt1 => term1%nd%volt
    this%volt2 => term2%nd%volt
    this%curr1 => term1%curr
    this%curr2 => term2%curr

    ! set main variable
    call this%init_f(this%curr1)

    ! depend on voltages and curr1; provide curr2
    ivolt1 = this%depend(this%volt1)
    ivolt2 = this%depend(this%volt2)
    icurr1 = this%depend(this%curr1)
    icurr2 = this%provide(this%curr2)

    ! init voltage jacobians
    this%jaco_volt(1)%p => this%init_jaco_f(ivolt1, const = .true., dtime = .true.)
    this%jaco_volt(2)%p => this%init_jaco_f(ivolt2, const = .true., dtime = .true.)
    call this%jaco_volt(1)%p%set(dum, dum,  C)
    call this%jaco_volt(2)%p%set(dum, dum, -C)

    ! init current jacobians
    this%jaco_curr(1)%p => this%init_jaco_f(      icurr1, const = .true.)
    this%jaco_curr(2)%p => this%init_jaco(icurr2, icurr1, const = .true.)
    call this%jaco_curr(1)%p%set(dum, dum, -1.0)
    call this%jaco_curr(2)%p%set(dum, dum, -1.0)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine capacitor_equ_eval(this)
    !! evaluate capacitor equation
    class(capacitor_equ), intent(inout) :: this

    call this%f%set(- this%curr1%get())
    call this%curr2%set(- this%curr1%get())
  end subroutine

  subroutine capacitor_init(this, crt, name, node1, node2, C)
    !! initialize capacitor
    class(capacitor), intent(out)   :: this
    type(circuit),    intent(inout) :: crt
      !! circuit to add this to
    character(*),     intent(in)    :: name
      !! component name
    character(*),     intent(in)    :: node1
      !! name of first node this is connected to
    character(*),     intent(in)    :: node2
      !! name of first node this is connected to
    real,             intent(in)    :: C
      !! capacitance

    ! create new component
    this%comp => crt%add_component(name, [new_string(node1), new_string(node2)])

    ! initialize equation and add to system
    call this%equ%init(name, this%comp%terms(1), this%comp%terms(2), C)
    call crt%sys%add_equation(this%equ)
  end subroutine

  subroutine capacitor_destruct(this)
    !! destruct capacitor
    class(capacitor), intent(inout) :: this

    call this%equ%destruct()
  end subroutine

end module
