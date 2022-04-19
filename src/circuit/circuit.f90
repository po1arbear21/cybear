m4_include(../util/macro.f90.inc)

module circuit_m

  use error_m,        only: program_error
  use esystem_m,      only: esystem
  use current_m,      only: current
  use jacobian_m,     only: jacobian_ptr
  use map_m,          only: map_string_int, mapnode_string_int
  use res_equation_m, only: res_equation
  use string_m,       only: string, new_string
  use vector_m,       only: vector_int
  use voltage_m,      only: voltage

  implicit none

  private
  public node, terminal, component, circuit

  type node_ptr
    !! pointer to circuit node
    type(node), pointer :: p => null()
  end type

  type terminal_ptr
    !! pointer to circuit component terminal
    type(terminal), pointer :: p => null()
  end type

  type component_ptr
    !! polymorphic pointer to component
    class(component), pointer :: p => null()
  end type

  m4_define({T},{node_ptr})
  m4_include(../util/vector_def.f90.inc)

  m4_define({T},{terminal_ptr})
  m4_include(../util/vector_def.f90.inc)

  m4_define({T},{component_ptr})
  m4_include(../util/vector_def.f90.inc)

  type, extends(res_equation) :: kcl_equ
    !! Kirchhoff Current Law equation

    type(node), pointer :: nd => null()
      !! pointer to circuit node

    type(jacobian_ptr), allocatable :: jaco_curr(:)
      !! jacobians wrt currents
  contains
    procedure :: init => kcl_equ_init
    procedure :: eval => kcl_equ_eval
  end type

  type node
    !! circuit node

    integer                   :: inode
      !! node index
    type(voltage)             :: volt
      !! node voltage variable
    type(vector_terminal_ptr) :: terms
      !! connected terminals
    type(kcl_equ)             :: kcl
      !! Kirchhoff Current Law equation
  contains
    procedure :: init       => node_init
    procedure :: connect    => node_connect
    procedure :: init_final => node_init_final
    procedure :: destruct   => node_destruct
  end type

  type terminal
    !! circuit component terminal

    integer             :: iterm
      !! terminal index
    type(current)       :: curr
      !! current variable
    type(node), pointer :: nd => null()
      !! connected node
  contains
    procedure :: init => terminal_init
  end type

  type, abstract :: component
    !! circuit component

    character(:),   allocatable :: name
      !! name of component (unique)
    integer                     :: icomp
      !! component index
    type(terminal), allocatable :: terms(:)
      !! current terminals
  contains
    procedure :: init => component_init
    procedure(component_init_final), deferred :: init_final
    procedure(component_destruct),   deferred :: destruct
  end type

  abstract interface
    subroutine component_init_final(this, sys)
      import component, esystem
      class(component), target, intent(inout) :: this
      type(esystem),            intent(inout) :: sys
    end subroutine
    subroutine component_destruct(this)
      import component
      class(component), intent(inout) :: this
    end subroutine
  end interface

  type circuit
    !! electrical circuit

    character(:), allocatable  :: name
    type(vector_component_ptr) :: comps
      !! circuit components
    type(map_string_int)       :: comp_names
      !! component names
    type(vector_node_ptr)      :: nodes
      !! circuit nodes
    type(vector_int)           :: free_nodes
      !! list of free circuit nodes (previously deleted)
    type(esystem)              :: sys
      !! equation system
  contains
    procedure :: init       => circuit_init
    procedure :: add        => circuit_add
    procedure :: get        => circuit_get
    procedure :: connect    => circuit_connect
    procedure :: init_final => circuit_init_final
    procedure :: destruct   => circuit_destruct
  end type

contains

  m4_define({T},{node_ptr})
  m4_include(../util/vector_imp.f90.inc)

  m4_define({T},{terminal_ptr})
  m4_include(../util/vector_imp.f90.inc)

  m4_define({T},{component_ptr})
  m4_include(../util/vector_imp.f90.inc)

  subroutine kcl_equ_init(this, nd)
    !! initialize KCL equation
    class(kcl_equ),     intent(out) :: this
    type(node), target, intent(in)  :: nd
      !! circuit node

    integer       :: iterm, dum(0)
    character(16) :: name

    ! init base
    write (name, "(A,I0)") "KCL", nd%inode
    call this%equation_init(trim(name))

    ! save pointer to circuit node
    this%nd => nd

    ! use node voltage as main variable
    call this%init_f(nd%volt)

    ! depend on currents from all terminals (f = sum I = 0)
    allocate (this%jaco_curr(nd%terms%n))
    do iterm = 1, nd%terms%n
      this%jaco_curr(iterm)%p => this%init_jaco_f(this%depend(nd%terms%d(iterm)%p%curr), const = .true.)
      call this%jaco_curr(iterm)%p%set(dum, dum, 1.0)
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine kcl_equ_eval(this)
    !! evaluate KCL equation
    class(kcl_equ), intent(inout) :: this

    integer :: iterm
    real    :: tmp(1)

    tmp = 0
    do iterm = 1, this%nd%terms%n
      tmp = tmp + this%nd%terms%d(iterm)%p%curr%x
    end do
    call this%f%set(tmp)
  end subroutine

  subroutine node_init(this, inode)
    !! initialize circuit node
    class(node),  intent(out) :: this
    integer,      intent(in)  :: inode
      !! node index

    character(16) :: name

    ! set node index
    this%inode = inode

    ! init voltage variable
    write (name, "(A,I0)") "V", inode
    call this%volt%init(trim(name))

    ! init terminals vector
    call this%terms%init(0, c = 4)
  end subroutine

  subroutine node_connect(this, term)
    !! connect node to component terminal
    class(node),    target, intent(inout) :: this
    type(terminal), target, intent(inout) :: term

    type(terminal_ptr) :: tp

    tp%p => term
    call this%terms%push(tp)
    term%nd => this
  end subroutine

  subroutine node_init_final(this, sys)
    !! initialize KCL equation and add to system
    class(node), target, intent(inout) :: this
    type(esystem),       intent(inout) :: sys
      !! equation system

    call this%kcl%init(this)
    call sys%add_equation(this%kcl)
  end subroutine

  subroutine node_destruct(this)
    !! destruct circuit node
    class(node), intent(inout) :: this

    call this%kcl%destruct()
  end subroutine

  subroutine terminal_init(this, comp_name, iterm)
    !! initialize circuit component terminal
    class(terminal), intent(out) :: this
    character(*),    intent(in)  :: comp_name
      !! component name
    integer,         intent(in)  :: iterm
      !! terminal index

    character(len(comp_name)+16) :: curr_name

    ! set terminal index
    this%iterm = iterm

    ! init current
    write (curr_name, "(A,I0)") trim(comp_name)//".I", iterm
    call this%curr%init(trim(curr_name))
  end subroutine

  subroutine component_init(this, crt, name, nterm)
    !! initialize circuit component
    class(component), target, intent(out)   :: this
    type(circuit),            intent(inout) :: crt
      !! circuit
    character(*),             intent(in)    :: name
      !! component name
    integer,                  intent(in)    :: nterm
      !! number of terminals

    integer :: iterm

    ! save name
    this%name = name

    ! initialize terminals
    allocate (this%terms(nterm))
    do iterm = 1, nterm
      call this%terms(iterm)%init(name, iterm)
    end do

    ! add to circuit
    this%icomp = crt%add(this)
  end subroutine

  subroutine circuit_init(this, name)
    !! initialize circuit
    class(circuit), intent(out) :: this
    character(*),   intent(in)  :: name
      !! circuit name

    this%name = name
    call this%comps%init(0, c = 8)
    call this%comp_names%init()
    call this%nodes%init(     0, c = 8)
    call this%free_nodes%init(0, c = 8)
  end subroutine

  function circuit_add(this, comp) result(icomp)
    !! add component to circuit
    class(circuit),            intent(inout) :: this
    class(component), pointer, intent(in)    :: comp
      !! component to add
    integer                                  :: icomp
      !! return component index

    logical :: status

    call this%comps%push(component_ptr(comp))
    icomp = this%comps%n

    call this%comp_names%insert(new_string(comp%name), icomp, status = status)
    if (.not. status) call program_error("component with name "//comp%name//" already exists")
  end function

  function circuit_get(this, name) result(comp)
    !! get component by name
    class(circuit), intent(in) :: this
    character(*),   intent(in) :: name
      !! component name
    class(component), pointer  :: comp
      !! return pointer to component or null if not found

    type(mapnode_string_int), pointer :: p

    ! lookup name in map
    p => this%comp_names%find(new_string(name))

    ! return component or null
    if (associated(p)) then
      comp => this%comps%d(p%value)%p
    else
      nullify(comp)
    end if
  end function

  subroutine circuit_connect(this, term1, term2)
    !! connect two terminals, automatically create or merge nodes if necessary
    class(circuit), target, intent(inout) :: this
    type(terminal), target, intent(inout) :: term1
      !! first terminal
    type(terminal), target, intent(inout) :: term2
      !! second terminal

    integer             :: iterm, inode
    type(node), pointer :: node1, node2
    type(node_ptr)      :: ptr

    ! get nodes
    node1 => term1%nd
    node2 => term2%nd

    if (associated(node1) .and. associated(node2)) then
      if (associated(node1, node2)) return

      ! merge existing nodes (connect all terminals of node2 to node1)
      do iterm = 1, node2%terms%n
        call node1%connect(node2%terms%d(iterm)%p)
      end do

      ! delete node2
      call this%free_nodes%push(node2%inode)
      call node2%destruct()
      deallocate (this%nodes%d(node2%inode)%p)
    elseif (associated(node1)) then
      call node1%connect(term2)
    elseif (associated(node2)) then
      call node2%connect(term1)
    else
      ! get free node index
      if (this%free_nodes%n > 0) then
        inode = this%free_nodes%back()
        call this%free_nodes%pop()
      else
        call this%nodes%push(ptr)
        inode = this%nodes%n
      end if

      ! create new node
      allocate (node1)
      call node1%init(inode)
      this%nodes%d(inode)%p => node1

      ! connect new node to both terminals
      call node1%connect(term1)
      call node1%connect(term2)
    end if
  end subroutine

  subroutine circuit_init_final(this)
    !! finish circuit initialization (create equation system)
    class(circuit), intent(inout) :: this

    integer                   :: icomp, iterm, inode
    class(component), pointer :: comp
    type(node),       pointer :: nd

    ! initialize equation system
    call this%sys%init(this%name)

    ! add component equations to system
    do icomp = 1, this%comps%n
      comp => this%comps%d(icomp)%p
      do iterm = 1, size(comp%terms)
        if (.not. associated(comp%terms(iterm)%nd)) then
          print "(A,A)", "component name: ", comp%name
          print "(A,I0)", "terminal index: ", iterm
          call program_error("component not connected")
        end if
      end do
      call comp%init_final(this%sys)
    end do

    ! add node equations to system
    do inode = 1, this%nodes%n
      nd => this%nodes%d(inode)%p
      if (.not. associated(nd)) cycle
      call nd%init_final(this%sys)
    end do

    ! finish system
    call this%sys%init_final()

    call this%sys%g%output(this%name)
  end subroutine

  subroutine circuit_destruct(this)
    !! destruct electrical circuit (deallocate all nodes and components)
    class(circuit), intent(inout) :: this

    integer :: i

    if (allocated(this%comps%d)) then
      do i = 1, this%comps%n
        if (associated(this%comps%d(i)%p)) then
          call this%comps%d(i)%p%destruct()
          deallocate (this%comps%d(i)%p)
        end if
      end do
    end if
    if (allocated(this%nodes%d)) then
      do i = 1, this%nodes%n
        if (associated(this%nodes%d(i)%p)) then
          call this%nodes%d(i)%p%destruct()
          deallocate (this%nodes%d(i)%p)
        end if
      end do
    end if
    call this%sys%destruct()
  end subroutine

end module
