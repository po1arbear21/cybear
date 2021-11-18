module ground_m

  use circuit_m,      only: circuit, component, terminal
  use current_m,      only: current
  use res_equation_m, only: res_equation
  use esystem_m,      only: esystem
  use jacobian_m,     only: jacobian
  use string_m,       only: new_string
  use voltage_m,      only: voltage

  implicit none

  private
  public ground, new_ground

  type, extends(res_equation) :: ground_equ
    !! ground equation

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

  type, extends(component) :: ground
    !! circuit component: ground

    type(ground_equ) :: equ
      !! ground equation
  contains
    procedure :: init_final => ground_init_final
    procedure :: destruct   => ground_destruct
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

    ! f = V = 0
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

  function new_ground(crt, name) result(gnd)
    !! create new ground and add to circuit
    type(circuit), intent(inout) :: crt
      !! circuit
    character(*),  intent(in)    :: name
      !! ground name
    type(ground), pointer        :: gnd
      !! return pointer to new ground

    ! allocate and initialize component
    allocate (gnd)
    call gnd%init(crt, name, 1)
  end function

  subroutine ground_init_final(this, sys)
    !! initialize ground equation and add to equation system
    class(ground), target, intent(inout) :: this
    type(esystem),         intent(inout) :: sys

    ! initialize equation
    call this%equ%init(this%name, this%terms(1))
    call sys%add_equation(this%equ)
  end subroutine

  subroutine ground_destruct(this)
    !! destruct ground component
    class(ground), intent(inout) :: this

    call this%equ%destruct()
  end subroutine

end module
