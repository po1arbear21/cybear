module csource_m

  use circuit_m,      only: circuit, component, terminal
  use current_m,      only: current
  use equation_m,     only: equation
  use jacobian_m,     only: jacobian_ptr
  use res_equation_m, only: res_equation
  use string_m,       only: new_string
  use voltage_m,      only: voltage

  implicit none

  private
  public csource

  type, extends(equation) :: csource_equ
    !! current source equation (I1 = I0; I2 = -I0)

    type(current), pointer :: curr0 => null()
      !! input current
    type(current), pointer :: curr1 => null()
      !! first terminal current (= curr0)
    type(current), pointer :: curr2 => null()
      !! second terminal current (= -curr0)

    type(jacobian_ptr) :: jaco_curr(2)
      !! jacobian of curr1/curr2 wrt curr0
  contains
    procedure :: init => csource_equ_init
    procedure :: eval => csource_equ_eval
  end type

  type csource
    !! current source

    type(current)            :: I
      !! input current
    type(component), pointer :: comp => null()
      !! pointer to circuit component
    type(csource_equ)        :: equ
      !! current source equation
  contains
    procedure :: init => csource_init
  end type

contains

  subroutine csource_equ_init(this, name, term1, term2, curr)
    !! initialize current source equation
    class(csource_equ),     intent(out) :: this
    character(*),           intent(in)  :: name
      !! resistor name
    type(terminal), target, intent(in)  :: term1
      !! first terminal
    type(terminal), target, intent(in)  :: term2
      !! second terminal
    type(current),  target, intent(in)  :: curr
      !! input current

    integer :: icurr0, icurr1, icurr2, dum(0)

    ! init base
    call this%equation_init(name)

    ! save pointers
    this%curr0 => curr
    this%curr1 => term1%curr
    this%curr2 => term2%curr

    ! depend on curr0, provide curr1 and curr2
    icurr0 = this%depend(this%curr0)
    icurr1 = this%provide(this%curr1)
    icurr2 = this%provide(this%curr2)

    ! init jacobians
    this%jaco_curr(1)%p => this%init_jaco(icurr1, icurr0, const = .true.)
    this%jaco_curr(2)%p => this%init_jaco(icurr2, icurr0, const = .true.)

    ! I1 = I0; I2 = - I0
    call this%jaco_curr(1)%p%set(dum, dum,  1.0)
    call this%jaco_curr(2)%p%set(dum, dum, -1.0)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine csource_equ_eval(this)
    !! evaluate resistor equation
    class(csource_equ), intent(inout) :: this

    call this%curr1%set(  this%curr0%get())
    call this%curr2%set(- this%curr0%get())
  end subroutine

  subroutine csource_init(this, crt, name, node1, node2)
    !! initialize current source
    class(csource), intent(out)   :: this
    type(circuit),  intent(inout) :: crt
      !! circuit to add this to
    character(*),   intent(in)    :: name
      !! component name
    character(*),   intent(in)    :: node1
      !! name of first node this is connected to
    character(*),   intent(in)    :: node2
      !! name of first node this is connected to

    ! init current variable and provide as input variable
    call this%I%init(name // ".I")
    call crt%sys%provide(this%I, input = .true.)

    ! create new component
    this%comp => crt%add_component(name, [new_string(node1), new_string(node2)])

    ! initialize equation and add to system
    call this%equ%init(name, this%comp%terms(1), this%comp%terms(2), this%I)
    call crt%sys%add_equation(this%equ)
  end subroutine

end module
