module inductor_m

  use circuit_m,      only: circuit, component, terminal
  use current_m,      only: current
  use esystem_m,      only: esystem
  use jacobian_m,     only: jacobian_ptr
  use res_equation_m, only: res_equation
  use string_m,       only: new_string
  use voltage_m,      only: voltage

  implicit none

  private
  public inductor, new_inductor

  type, extends(res_equation) :: inductor_equ
    !! inductor equation

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
    procedure :: init => inductor_equ_init
    procedure :: eval => inductor_equ_eval
  end type

  type, extends(component) :: inductor
    !! circuit component: inductor

    real :: L
      !! inductance

    type(inductor_equ) :: equ
      !! inductor equation
  contains
    procedure :: init_final => inductor_init_final
    procedure :: destruct   => inductor_destruct
  end type

contains

  subroutine inductor_equ_init(this, name, term1, term2, L)
    !! initialize inductor equation
    class(inductor_equ),    intent(out) :: this
    character(*),           intent(in)  :: name
      !! inductor name
    type(terminal), target, intent(in)  :: term1
      !! first terminal
    type(terminal), target, intent(in)  :: term2
      !! second terminal
    real,           target, intent(in)  :: L
      !! inductance

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
    this%jaco_volt(1)%p => this%init_jaco_f(ivolt1, const = .true.)
    this%jaco_volt(2)%p => this%init_jaco_f(ivolt2, const = .true.)
    call this%jaco_volt(1)%p%set(dum, dum, -1.0)
    call this%jaco_volt(2)%p%set(dum, dum,  1.0)

    ! init current jacobians
    this%jaco_curr(1)%p => this%init_jaco_f(      icurr1, const = .true., dtime = .true.)
    this%jaco_curr(2)%p => this%init_jaco(icurr2, icurr1, const = .true.)
    call this%jaco_curr(1)%p%set(dum, dum, L)
    call this%jaco_curr(2)%p%set(dum, dum, -1.0)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine inductor_equ_eval(this)
    !! evaluate inductor equation
    class(inductor_equ), intent(inout) :: this

    call this%f%set(this%volt2%get() - this%volt1%get())
    call this%curr2%set(- this%curr1%get())
  end subroutine

  function new_inductor(crt, name, L) result(ind)
    !! create new inductor and add to circuit
    type(circuit), intent(inout) :: crt
      !! circuit
    character(*),  intent(in)    :: name
      !! inductor name
    real,          intent(in)    :: L
      !! inductance
    type(inductor), pointer     :: ind
      !! return pointer to new inductor

    ! allocate and initialize component
    allocate (ind)
    call ind%init(crt, name, 2)

    ! set inductance
    ind%L = L
  end function

  subroutine inductor_init_final(this, sys)
    !! initialize inductor equation and add to equation system
    class(inductor), target, intent(inout) :: this
    type(esystem),           intent(inout) :: sys

    ! initialize equation
    call this%equ%init(this%name, this%terms(1), this%terms(2), this%L)
    call sys%add_equation(this%equ)
  end subroutine

  subroutine inductor_destruct(this)
    !! destruct inductor
    class(inductor), intent(inout) :: this

    call this%equ%destruct()
  end subroutine

end module
