module ground_m

  use circuit_m,      only: circuit, component, terminal
  use current_m,      only: current
  use jacobian_m,     only: jacobian
  use res_equation_m, only: res_equation
  use string_m,       only: new_string
  use voltage_m,      only: voltage

  implicit none

  private
  public ground

  type, extends(res_equation) :: ground_equ
    !! circuit ground equation

    type(voltage), pointer :: volt => null()
      !! node voltage
    type(current), pointer :: curr => null()
      !! ground current

    type(jacobian), pointer :: jaco_volt
      !! jacobians wrt voltages
  contains
    procedure :: init => ground_equ_init
    procedure :: eval => ground_equ_eval
  end type

  type ground
    type(component), pointer :: comp => null()
      !! pointer to circuit component
    type(ground_equ) :: equ
      !! ground equation
  contains
    procedure :: init => ground_init
  end type

contains

  subroutine ground_equ_init(this, name, term)
    !! initialize ground equation
    class(ground_equ),      intent(out) :: this
    character(*),           intent(in)  :: name
      !! ground name
    type(terminal), target, intent(in)  :: term
      !! terminal

    integer :: dum(0)

    ! init base
    call this%equation_init(name)

    ! save pointers to variables
    this%volt => term%nd%volt
    this%curr => term%curr

    ! set main variable
    call this%init_f(this%curr)

    ! f = V == 0
    this%jaco_volt => this%init_jaco_f(this%depend(this%volt), const = .true.)
    call this%jaco_volt%set(dum, dum, 1.0)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine ground_equ_eval(this)
    !! evaluate ground equation
    class(ground_equ), intent(inout) :: this

    call this%f%set(this%volt%get())
  end subroutine

  subroutine ground_init(this, crt, name, node1)
    !! initialize ground
    class(ground), intent(out)   :: this
    type(circuit), intent(inout) :: crt
      !! circuit to add this to
    character(*),  intent(in)    :: name
      !! component name
    character(*),  intent(in)    :: node1
      !! name of node this is connected to

    ! create new component
    this%comp => crt%add_component(name, [new_string(node1)])

    ! initialize equation and add to system
    call this%equ%init(name, this%comp%terms(1))
    call crt%sys%add_equation(this%equ)
  end subroutine

end module
