m4_include(../util/macro.f90.inc)

module circuit_m

  use error_m,        only: assert_failed, program_error
  use esystem_m,      only: esystem
  use current_m,      only: current
  use hashmap_m,      only: hashmap_int
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
    class(node), pointer :: p => null()
  end type

  type terminal_ptr
    !! pointer to circuit component terminal
    class(terminal), pointer :: p => null()
  end type

  type component_ptr
    !! pointer to component
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

    character(:), allocatable :: name
      !! name of node (unique)
    integer                   :: inode
      !! node index
    type(voltage)             :: volt
      !! node voltage variable
    type(vector_terminal_ptr) :: terms
      !! list of connected terminals
    type(kcl_equ)             :: kcl
      !! Kirchhoff Current Law equation (sum_k I_k = 0)
  contains
    procedure :: init       => node_init
    procedure :: connect    => node_connect
    procedure :: init_final => node_init_final
  end type

  type terminal
    !! circuit component terminal

    integer             :: iterm
      !! terminal index
    type(node), pointer :: nd => null()
      !! connected node
    type(current)       :: curr
      !! current variable
  contains
    procedure :: init    => terminal_init
    procedure :: get_ptr => terminal_get_ptr
  end type

  type component
    !! circuit component

    character(:),   allocatable :: name
      !! name of component (unique)
    integer                     :: icomp
      !! component index
    type(terminal), allocatable :: terms(:)
      !! current terminals
  contains
    procedure :: init => component_init
  end type

  type circuit
    !! electrical circuit

    character(:), allocatable  :: name
      !! circuit name

    type(vector_component_ptr) :: comps
      !! circuit components
    type(map_string_int)       :: comp_map
      !! component name => component index

    type(vector_node_ptr) :: nodes
      !! circuit nodes
    type(map_string_int)  :: node_map
      !! node name => node index

    type(esystem) :: sys
      !! equation system including all node and component equations
  contains
    procedure :: init          => circuit_init
    procedure :: add_component => circuit_add_component
    procedure :: init_final    => circuit_init_final

    final :: circuit_destruct
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

    integer :: iterm, idep, dum(0)

    ! init base
    call this%equation_init(nd%name // ".KCL")

    ! save pointer to circuit node
    this%nd => nd

    ! use node voltage as main variable
    call this%init_f(nd%volt)

    ! depend on currents from all terminals (f = sum I == 0)
    allocate (this%jaco_curr(nd%terms%n))
    do iterm = 1, nd%terms%n
      idep = this%depend(nd%terms%d(iterm)%p%curr)
      this%jaco_curr(iterm)%p => this%init_jaco_f(idep, const = .true.)
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

    ! set residual to sum of all currents
    tmp = 0
    do iterm = 1, this%nd%terms%n
      tmp = tmp + this%nd%terms%d(iterm)%p%curr%x
    end do
    call this%f%set(tmp)
  end subroutine

  subroutine node_init(this, name, inode)
    !! initialize circuit node
    class(node),  intent(out) :: this
    character(*), intent(in)  :: name
      !! node name
    integer,      intent(in)  :: inode
      !! node index

    ! set members (KCL equation is initialized in init_final, after terminals are connected)
    this%name  = name
    this%inode = inode
    call this%volt%init(name // ".V")
    call this%terms%init(0, c = 4)
  end subroutine

  subroutine node_connect(this, term)
    !! connect node to terminal
    class(node),            intent(inout) :: this
    type(terminal), target, intent(in)    :: term
      !! terminal to connect to

    ! save pointer to terminal in list
    call this%terms%push(term%get_ptr())
  end subroutine

  subroutine node_init_final(this)
    !! finish node initialization
    class(node), target, intent(inout) :: this

    ! initialize KCL equation
    call this%kcl%init(this)
  end subroutine

  subroutine terminal_init(this, comp_name, iterm, nd)
    !! initialize component terminal
    class(terminal),    intent(out)   :: this
    character(*),       intent(in)    :: comp_name
      !! name of circuit component
    integer,            intent(in)    :: iterm
      !! terminal index
    type(node), target, intent(inout) :: nd
      !! connected circuit node

    character(len(comp_name)+16) :: curr_name

    this%iterm = iterm
    this%nd => nd

    ! initialize current variable
    write (curr_name, "(A,I0)") comp_name // ".I", iterm
    call this%curr%init(trim(curr_name))

    ! connect node to this terminal
    call nd%connect(this)
  end subroutine

  function terminal_get_ptr(this) result(ptr)
    !! returns pointer type to this terminal
    class(terminal), target, intent(in) :: this
    type(terminal_ptr)                  :: ptr

    ptr%p => this
  end function

  subroutine component_init(this, name, icomp, nodes)
    !! initialize circuit component
    class(component), intent(out) :: this
    character(*),     intent(in)  :: name
      !! component name
    integer,          intent(in)  :: icomp
      !! component index
    type(node_ptr),   intent(in)  :: nodes(:)
      !! pointers to connnected nodes

    integer :: iterm

    this%name  = name
    this%icomp = icomp

    ! initialize terminals
    allocate (this%terms(size(nodes)))
    do iterm = 1, size(nodes)
      call this%terms(iterm)%init(name, iterm, nodes(iterm)%p)
    end do
  end subroutine

  subroutine circuit_init(this, name)
    !! initialize circuit
    class(circuit), intent(out) :: this
    character(*),   intent(in)  :: name
      !! name of circuit

    this%name = name

    call this%comps%init(0, c = 8)
    call this%comp_map%init()

    call this%nodes%init(0, c = 8)
    call this%node_map%init()

    call this%sys%init(name)
  end subroutine

  function circuit_add_component(this, comp_name, node_names) result(comp)
    !! add new component to circuit
    class(circuit), intent(inout) :: this
    character(*),   intent(in)    :: comp_name
      !! name of new component
    type(string),   intent(in)    :: node_names(:)
      !! names of connected nodes
    type(component), pointer      :: comp
      !! return pointer to newly created component

    integer                           :: i
    logical                           :: status
    type(node_ptr), allocatable       :: n(:)
    type(mapnode_string_int), pointer :: p
    type(component_ptr)               :: cptr

    ! get node pointers
    allocate (n(size(node_names)))
    do i = 1, size(node_names)
      ! lookup name of node in map
      p => this%node_map%find(node_names(i))

      if (associated(p)) then
        ! set pointer to existing node
        n(i)%p => this%nodes%d(p%value)%p
      else
        ! create new node
        allocate (n(i)%p)
        call this%nodes%push(n(i))
        call n(i)%p%init(node_names(i)%s, this%nodes%n)
        call this%node_map%insert(node_names(i), this%nodes%n)
      end if
    end do

    ! create new component
    allocate (cptr%p)
    comp => cptr%p
    call this%comps%push(cptr)
    call comp%init(comp_name, this%comps%n, n)
    call this%comp_map%insert(new_string(comp_name), this%comps%n, status = status)
    m4_assert(status)
  end function

  subroutine circuit_init_final(this)
    !! finish circuit initialization
    class(circuit), intent(inout) :: this

    integer :: i

    ! finish initialization of nodes and add equations to system
    do i = 1, this%nodes%n
      call this%nodes%d(i)%p%init_final()
      call this%sys%add_equation(this%nodes%d(i)%p%kcl)
    end do

    ! finish esystem initialization
    call this%sys%init_final()
  end subroutine

  subroutine circuit_destruct(this)
    !! circuit destructor
    type(circuit), intent(inout) :: this

    integer :: i

    if (allocated(this%comps%d)) then
      do i = 1, this%comps%n
        if (associated(this%comps%d(i)%p)) deallocate (this%comps%d(i)%p)
      end do
    end if

    if (allocated(this%nodes%d)) then
      do i = 1, this%nodes%n
        if (associated(this%nodes%d(i)%p)) deallocate (this%nodes%d(i)%p)
      end do
    end if
  end subroutine

end module
