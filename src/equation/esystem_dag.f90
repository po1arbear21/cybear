module esystem_dag_m
  use res_equation_m
  use jacobian_chain_m
  use vector_m
  implicit none

  ! node status
  integer, parameter :: NDSTATUS_DEP   = 0 ! created as dependency (temporary)
  integer, parameter :: NDSTATUS_PROV  = 1 ! provided
  integer, parameter :: NDSTATUS_MAIN  = 2 ! main var
  integer, parameter :: NDSTATUS_RES   = 3 ! residual

  type dag_node_ptr
    type(dag_node), pointer :: p  => null()
  end type

#define T dag_node_ptr
#define TT type(dag_node_ptr)
#include "../util/vector_def.f90.inc"

  type dag_node
    !! node in directed acyclic graph

    integer :: id
      !! index for dag%nodes
    integer :: mvar_id
      !! index for dag%mvari (NDSTATUS_MAIN only)
    integer :: status
      !! node status

    logical :: visited
      !! visited flag for traversal
    logical :: analyzed
      !! analyzed flag for traversal
    logical :: const
      !! constant flag

    class(vselector), pointer :: v => null()
      !! pointer to var selector
    type(dag_equ),    pointer :: e => null()
      !! pointer to DAG equation

    type(vector_dag_node_ptr) :: parents
      !! parent nodes
    type(vector_dag_node_ptr) :: parents_t
      !! parent nodes for time derivatives of vars (NDSTATUS_RES only)

    type(vector_jacobian_matrix_ptr) :: partial_jaco
      !! partial derivatives wrt parents
    type(vector_jacobian_matrix_ptr) :: partial_jaco_t
      !! partial derivatives wrt parents_t (NDSTATUS_RES only)

    type(jacobian_matrix_ptr), allocatable :: total_jaco(:)
      !! total derivatives wrt main vars
    type(jacobian_matrix_ptr), allocatable :: total_jaco_t(:)
      !! total derivatives wrt time derivatives of main vars (NDSTATUS_RES only)

    type(vector_jacobian_chain_ptr), allocatable :: jchain(:)
      !! jacobian chains for computation of total derivatives wrt main vars
    type(vector_jacobian_chain_ptr), allocatable :: jchain_t(:)
      !! jacobian chains for computation of total derivatives wrt time der. of main vars
  contains
    procedure :: init     => dag_node_init
    procedure :: destruct => dag_node_destruct
    procedure :: analyze  => dag_node_analyze
    procedure :: eval     => dag_node_eval
  end type

  type dag_equ_ptr
    type(dag_equ), pointer :: p => null()
  end type

#define T dag_equ_ptr
#define TT type(dag_equ_ptr)
#include "../util/vector_def.f90.inc"

  type dag_equ
    !! equation with additional DAG information
    integer                         :: id
      !! index for dag%equs
    class(equation),    pointer     :: e    => null()
      !! pointer to equation
    type(dag_node),     pointer     :: main => null()
      !! pointer to main var (only used for residual equation)
    type(dag_node),     pointer     :: res  => null()
      !! pointer to residual (only used for residual equation)
    type(dag_node_ptr), allocatable :: prov(:)
      !! pointer to provided nodes
    logical                         :: evaluated
      !! evaluated flag for traversal
  contains
    procedure :: init => dag_equ_init
  end type

  type dag
    ! directed acyclic graph
    type(vector_dag_equ_ptr)  :: equs
      !! equations
    type(vector_dag_node_ptr) :: nodes
      !! all nodes
    type(vector_int)          :: resi
      !! residual equation indices
    type(vector_int)          :: mvari
      !! main var indices
    type(vector_int)          :: ev
      !! evaluation list
  contains
    procedure :: init     => dag_init
    procedure :: destruct => dag_destruct
    procedure :: add_equ  => dag_add_equ
    procedure :: add_node => dag_add_node
    procedure :: analyze  => dag_analyze
  end type

contains

#define T dag_node_ptr
#define TT type(dag_node_ptr)
#include "../util/vector_imp.f90.inc"

#define T dag_equ_ptr
#define TT type(dag_equ_ptr)
#include "../util/vector_imp.f90.inc"

  recursive subroutine dag_node_init(this, id, mvar_id, v, e, prov_i, status, d)
    !! initialize DAG node
    class(dag_node),  target, intent(out)   :: this
    integer,                  intent(in)    :: id
      !! index for d%nodes
    integer,                  intent(in)    :: mvar_id
      !! index for d%mvari (only used for NDSTATUS_MAIN)
    class(vselector), target, intent(in)    :: v
      !! var selector
    class(dag_equ),   target, intent(in)    :: e
      !! DAG equation
    integer,                  intent(in)    :: prov_i
      !! index for e%e%vprov%d (only used for NDSTATUS_PROV)
    integer,                  intent(in)    :: status
      !! node status
    type(dag),        target, intent(inout) :: d
      !! DAG

    integer            :: i
    type(dag_node_ptr) :: np

    ! set members
    this%id       = id
    this%status   = status
    this%visited  = (status == NDSTATUS_MAIN) ! main vars are implicitly visited
    this%analyzed = (status == NDSTATUS_MAIN) ! main vars are implicitly analyzed
    this%const    =  .false.
    this%v        => v

    ! main var index
    if (status == NDSTATUS_MAIN) then
      this%mvar_id = mvar_id
    else
      this%mvar_id = 0
    end if

    ! set rest of members only for non-dependency node
    if (status == NDSTATUS_DEP) return

    ! save DAG equation pointer
    this%e => e

    ! init vectors
    call this%parents%init(    0, 8)
    call this%parents_t%init(  0, 8)
    call this%partial_jaco%init(  0, 8)
    call this%partial_jaco_t%init(0, 8)

    ! get parents, parents_t
    if (status == NDSTATUS_PROV) then
      ! collect dependencies
      do i = 1, e%e%vdep%n
        if (associated(e%e%jaco(prov_i,i)%p)) then
          ! add node as dependency to DAG
          np%p => d%add_node(e%e%vdep%d(i)%p, e, 0, NDSTATUS_DEP)

          ! save node pointer
          call this%parents%push(np)

          ! save jacobian matrix pointer for partial derivatives
          call this%partial_jaco%push(jacobian_matrix_ptr(e%e%jaco(prov_i,i)%p%matr))
        end if
      end do
    elseif (status == NDSTATUS_RES) then
      ! cast equation to residual equation
      select type (r => e%e)
        class is (res_equation)
          do i = 1, r%vdep%n
            if (associated(r%jaco_f(i)%p)) then
              ! add node as dependency to DAG
              np%p => d%add_node(r%vdep%d(i)%p, e, 0, NDSTATUS_DEP)

              ! save node pointer
              call this%parents%push(np)

              ! save jacobian matrix pointer for partial derivatives
              call this%partial_jaco%push(jacobian_matrix_ptr(r%jaco_f(i)%p%matr))
            end if
            if (associated(r%jaco_ft(i)%p)) then
              ! add node as dependency to DAG
              np%p => d%add_node(r%vdep%d(i)%p, e, 0, NDSTATUS_DEP)

              ! save node pointer
              call this%parents_t%push(np)

              ! save jacobian matrix pointer for partial derivatives
              call this%partial_jaco_t%push(jacobian_matrix_ptr(r%jaco_ft(i)%p%matr))
            end if
          end do

        class default
          call program_error("status = NDSTATUS_RES but equation is no residual equation")
      end select
    end if
  end subroutine

  subroutine dag_node_destruct(this)
    !! free memory
    class(dag_node), intent(inout) :: this

    integer :: i, j

    if (allocated(this%jchain)) then
      do i = 1, size(this%jchain)
        do j = 1, this%jchain(i)%n
          if (associated(this%jchain(i)%d(j)%p)) deallocate (this%jchain(i)%d(j)%p)
        end do
      end do
    end if

    if (allocated(this%jchain_t)) then
      do i = 1, size(this%jchain_t)
        do j = 1, this%jchain_t(i)%n
          if (associated(this%jchain_t(i)%d(j)%p)) deallocate (this%jchain_t(i)%d(j)%p)
        end do
      end do
    end if
  end subroutine

  recursive subroutine dag_node_analyze(this, d)
    !! initialize jacobian chains and total derivatives
    class(dag_node), target, intent(inout) :: this
    type(dag),       target, intent(inout) :: d

    integer :: i

    ! check if this node was already visited
    if (this%visited) then
      call program_error("circular dependency detected")
    else
      this%visited = .true.
    end if

    ! set constant flag, reset later if any parent is not constant (main vars are not constant and already analyzed)
    this%const = .true.

    ! analyze parents which are not already analyzed
    do i = 1, this%parents%n
      if (.not. this%parents%d(i)%p%analyzed) then
        call this%parents%d(i)%p%analyze(d)
      end if

      ! reset const flag if parent is not constant
      if (.not. this%parents%d(i)%p%const) then
        this%const = .false.
      end if
    end do

    ! analyze parents_t which are not already analyzed
    do i = 1, this%parents_t%n
      if (.not. this%parents_t%d(i)%p%analyzed) then
        call this%parents_t%d(i)%p%analyze(d)
      end if

      ! reset const flag if parent is not constant
      if (.not. this%parents_t%d(i)%p%const) then
        this%const = .false.
      end if
    end do

    ! init jacobian chains
    if (this%const) then
      ! evaluate provider equation if constant node (no derivatives needed, just provided values)
      call this%e%e%eval()
    else
      ! allocate total derivative jacobian matrices and chain vectors
      allocate (this%total_jaco(  d%mvari%n))
      allocate (this%total_jaco_t(d%mvari%n))
      allocate (this%jchain(      d%mvari%n))
      allocate (this%jchain_t(    d%mvari%n))

      ! init jacobian chain vectors
      do i = 1, d%mvari%n
        call this%jchain(i)%init(0, 8)
        call this%jchain_t(i)%init(0, 8)
      end do

      ! add jacobian chains for parents
      call init_jchains(this%parents, this%partial_jaco, this%total_jaco, this%jchain, .false.)

      ! add jacobian chains for parents_t
      call init_jchains(this%parents_t, this%partial_jaco_t, this%total_jaco_t, this%jchain_t, .true.)
    end if

    ! set analyzed flag
    this%analyzed = .true.

    ! analyze other nodes of equation
    if (associated(this%e%res)) then
      if (.not. this%e%res%analyzed) then
        call this%e%res%analyze(d)
      end if
    end if
    do i = 1, size(this%e%prov)
      if (.not. this%e%prov(i)%p%analyzed) then
        call this%e%prov(i)%p%analyze(d)
      end if
    end do

    ! add equation to evaluation list
    if (.not. this%e%evaluated) then
      call d%ev%push(this%e%id)
      this%e%evaluated = .true.
    end if

  contains

    subroutine init_jchains(this_parents, this_partial_jaco, this_total_jaco, this_jchain, time)
      type(vector_dag_node_ptr),        intent(in)    :: this_parents
        !! parents or parents_t
      type(vector_jacobian_matrix_ptr), intent(in)    :: this_partial_jaco
        !! partial_jaco or partial_jaco_t
      type(jacobian_matrix_ptr),        intent(inout) :: this_total_jaco(:)
        !! total_jaco or total_jaco_t
      type(vector_jacobian_chain_ptr),  intent(inout) :: this_jchain(:)
        !! jchain or jchain_t
      logical,                          intent(in)    :: time
        !! time derivative flag

      integer                                :: i, j, k
      type(jacobian_matrix),     pointer     :: n_total_jaco
      type(jacobian_matrix_ptr)              :: jmul(2)
      type(jacobian_matrix_ptr), allocatable :: jadd(:)
      type(jacobian_chain_ptr)               :: chain

      ! add multiplication chains for parents
      do i = 1, this_parents%n
        associate (n => this_parents%d(i)%p)
          if (n%const) then
            ! do nothing for const node
            cycle
          elseif (n%status == NDSTATUS_MAIN) then
            ! set total derivatives
            this_total_jaco(n%mvar_id)%p => this_partial_jaco%d(i)%p
          else
            ! chain rule for all main vars
            do j = 1, d%mvari%n
              if (time) then
                n_total_jaco => n%total_jaco_t(j)%p
              else
                n_total_jaco => n%total_jaco(j)%p
              end if

              ! do nothing if no dependency on this main var
              if (.not. associated(n_total_jaco)) cycle

              ! new multiplication chain
              nullify (chain%p)
              allocate (jacobian_mul_chain :: chain%p)

              ! init multiplication chain
              jmul(1)%p => this_partial_jaco%d(i)%p
              jmul(2)%p => n_total_jaco
              call chain%p%init(jmul)

              ! save chain
              call this_jchain(j)%push(chain)
            end do
          endif
        end associate
      end do

      ! add addition chains for parents (only after all multiplication chains have been added)
      do i = 1, d%mvari%n
        if (this_jchain(i)%n == 0) then ! no chain for this main var
          cycle
        elseif ((.not. associated(this_total_jaco(i)%p)) .and. this_jchain(i)%n == 1) then ! only one chain, set total derivatives
          this_total_jaco(i)%p => this_jchain(i)%d(1)%p%result
        else ! more than one chain, add results
          ! allocate temporary jacobian matrix pointers for addition
          if (associated(this_total_jaco(i)%p)) then
            ! partial derivative wrt main var exists
            allocate (jadd(1 + this_jchain(i)%n))
            jadd(1)%p => this_total_jaco(i)%p
            k = 1
          else
            ! no partial derivative wrt main var exists
            allocate (jadd(this_jchain(i)%n))
            k = 0
          end if

          ! save rest of jacobian matrix pointers
          do j = 1, this_jchain(i)%n
            jadd(k + j)%p => this_jchain(i)%d(j)%p%result
          end do

          ! init addition chain
          nullify (chain%p)
          allocate (jacobian_add_chain :: chain%p)
          call chain%p%init(jadd)

          ! deallocate temporary jacobian matrix pointers
          deallocate (jadd)

          ! save chain
          call this_jchain(i)%push(chain)

          ! set/update total derivative jacobian matrix
          this_total_jaco(i)%p => this_jchain(i)%d(this_jchain(i)%n)%p%result
        end if
      end do
    end subroutine

  end subroutine

  subroutine dag_node_eval(this)
    !! evaluate non-constant parts of jacobian chains
    class(dag_node), intent(inout) :: this ! DAG node

    integer :: i, j

    ! loop over main vars
    if (allocated(this%jchain)) then
      do i = 1, size(this%jchain)
        ! evaluate non-constant parts of jacobian chains
        do j = 1, this%jchain(i)%n
          call this%jchain(i)%d(j)%p%eval(const = .false.)
        end do
      end do
    end if
  end subroutine

  subroutine dag_equ_init(this, id, e, d)
    !! initialize DAG equation
    class(dag_equ),  target, intent(out)   :: this
    integer,                 intent(in)    :: id
      !! index for d%equs
    class(equation), target, intent(in)    :: e
      !! underlying equation
    type(dag),       target, intent(inout) :: d
      !! DAG

    integer :: i

    ! set members
    this%id        =  id
    this%e         => e
    this%evaluated =  .false.

    ! main var, residuals
    select type (e)
      class is (res_equation)
        this%main => d%add_node(e%mvar, this, 0, NDSTATUS_MAIN)
        this%res  => d%add_node(e%f,    this, 0, NDSTATUS_RES)
    end select

    ! provided (dependent var selectors are automatically added by dag_node_init)
    allocate (this%prov(e%vprov%n))
    do i = 1, e%vprov%n
      this%prov(i)%p => d%add_node(e%vprov%d(i)%p, this, i, NDSTATUS_PROV)
    end do
  end subroutine

  subroutine dag_init(this)
    !! initialize DAG (directed acyclic graph)
    class(dag), intent(out) :: this

    ! init vectors
    call this%equs%init( 0, 32)
    call this%nodes%init(0, 32)
    call this%resi%init( 0, 32)
    call this%mvari%init(0, 32)
    call this%ev%init(   0, 32)
  end subroutine

  subroutine dag_destruct(this)
    !! release memory
    class(dag), intent(inout) :: this

    integer :: i

    do i = 1, this%equs%n
      if (associated(this%equs%d(i)%p)) deallocate (this%equs%d(i)%p)
    end do

    do i = 1, this%nodes%n
      if (associated(this%nodes%d(i)%p)) then
        call this%nodes%d(i)%p%destruct()
        deallocate (this%nodes%d(i)%p)
      end if
    end do
  end subroutine

  subroutine dag_add_equ(this, e)
    !! add new equation to DAG
    class(dag),      target, intent(inout) :: this
    class(equation), target, intent(in)    :: e
      !! equation to add

    type(dag_equ_ptr) :: ep

    ! allocate DAG equation pointer
    allocate (ep%p)

    ! save DAG equation
    call this%equs%push(ep)

    ! init DAG equation
    call ep%p%init(this%equs%n, e, this)

    ! check for residual equation
    select type (e)
      class is (res_equation)
        call this%resi%push(this%equs%n)
    end select
  end subroutine

  recursive function dag_add_node(this, v, e, provi, status) result(n)
    !! add new node to DAG
    class(dag),       target, intent(inout) :: this
    class(vselector), target, intent(in)    :: v
      !! var selector
    type(dag_equ),    target, intent(in)    :: e
      !! DAG equation
    integer,                  intent(in)    :: provi
      !! index in e%vprov (only used for NDSTAT_PROV)
    integer,                  intent(in)    :: status
      !! node status
    type(dag_node), pointer                 :: n
      !! return pointer to node

    integer            :: i
    type(dag_node_ptr) :: np

    ! search for var selector
    do i = 1, this%nodes%n
      associate (m => this%nodes%d(i)%p)
        ! compare v with existing node
        if (.not. v%compare(m%v)) cycle

        ! node was found
        n => this%nodes%d(i)%p

        ! update node if not added as dependency
        if (status /= NDSTATUS_DEP) then
          ! can not add var selector multiple times as non-dependency
          if (m%status /= NDSTATUS_DEP) then
            print *, "vselector:"
            call v%print()
            print *, "equations:"
            print *, e%e%name
            print *, m%e%e%name
            call program_error("variable selector provided multiple times")
          end if

          ! update node
          if (status == NDSTATUS_MAIN) then
            call this%mvari%push(i)
            call n%init(i, this%mvari%n, v, e, 0, status, this)
          else
            call n%init(i, 0, v, e, provi, status, this)
          end if
        end if

        ! done (node was found and potentially updated)
        return
      end associate
    end do

    ! node not found, create new one and save pointer
    nullify (n)
    allocate (n)
    np%p => n
    call this%nodes%push(np)

    ! initialize node
    if (status == NDSTATUS_MAIN) then
      call this%mvari%push(this%nodes%n)
      call n%init(this%nodes%n, this%mvari%n, v, e, 0, status, this)
    else
      call n%init(this%nodes%n, 0, v, e, provi, status, this)
    end if
  end function

  subroutine dag_analyze(this)
    !! analyze dag, get eval list etc.
    class(dag), target, intent(inout) :: this ! directed acyclic graph (DAG)

    integer :: i

    ! analyze nodes, start with residuals
    do i = 1, this%resi%n
      associate (e => this%equs%d(this%resi%d(i))%p)
        if (.not. e%res%analyzed) then
          call e%res%analyze(this)
        end if
      end associate
    end do
  end subroutine

end module