module csource_m

  use circuit_m,  only: circuit, component, terminal
  use current_m,  only: current
  use esystem_m,  only: esystem
  use jacobian_m, only: jacobian, jacobian_ptr
  use math_m,     only: PI
  use equation_m, only: equation
  use string_m,   only: new_string

  implicit none

  private
  public csource, new_csource

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

  type, extends(component) :: csource
    !! circuit component: current source

    type(current) :: I
      !! input current

    type(csource_equ) :: equ
      !! current source equation
  contains
    procedure :: init_final => csource_init_final
    procedure :: destruct   => csource_destruct
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

  function new_csource(crt, name) result(csrc)
    !! create new current source
    type(circuit), intent(inout) :: crt
      !! circuit
    character(*),  intent(in)    :: name
      !! component name
    type(csource), pointer       :: csrc
      !! return pointer to new current source

    ! allocate and initialize component
    allocate (csrc)
    call csrc%init(crt, name, 2)

    ! init input current
    call csrc%I%init(name)
  end function

  subroutine csource_init_final(this, sys)
    !! initialize current source equation and add to equation system
    class(csource), target, intent(inout) :: this
    type(esystem),          intent(inout) :: sys

    call this%equ%init(this%name, this%terms(1), this%terms(2), this%I)
    call sys%add_equation(this%equ)
    call sys%provide(this%I, input = .true.)
  end subroutine

  subroutine csource_destruct(this)
    !! destruct voltage source
    class(csource), intent(inout) :: this

    call this%equ%destruct()
  end subroutine

end module
