module esystem_depgraph_m

  use equation_m,        only: equation
  use error_m,           only: program_error
  use hashmap_m,         only: hashmap_int
  use jacobian_chain_m,  only: jacobian_chain, jacobian_add_chain, jacobian_mul_chain, jacobian_chain_ptr, vector_jacobian_chain_ptr
  use jacobian_matrix_m, only: jacobian_matrix, jacobian_matrix_ptr, vector_jacobian_matrix_ptr
  use res_equation_m,    only: res_equation
  use vector_m,          only: vector_int
  use vselector_m,       only: vselector

  implicit none

  private
  public depgraph
  public STATUS_DEP, STATUS_PROV, STATUS_MAIN, STATUS_RES

  ! node status
  integer, parameter :: STATUS_DEP  = 0 ! created as dependency (temporary)
  integer, parameter :: STATUS_PROV = 1 ! provided
  integer, parameter :: STATUS_MAIN = 2 ! main var
  integer, parameter :: STATUS_RES  = 3 ! residual

  type node
    !! node in dependency graph
    integer :: id
      !! index for depgraph%nodes
    integer :: iimvar
      !! index for depgraph%imvar (depgraph%imvar%d(iimvar) == id)
    integer :: status
      !! node status

    logical :: visited
      !! visited flag for traversal
    logical :: analyzed
      !! analyzed flag for traversal
    logical :: const
      !! constant flag

    type(vselector), pointer :: v => null()
      !! pointer to corresponding vselector
    integer                  :: iequ
      !! index of depgraph equation this was added from

    type(vector_int) :: parents
      !! indices of parent nodes
    type(vector_int) :: parents_t
      !! indices of parent nodes for time derivatives (only used when status == STATUS_RES)

    type(vector_jacobian_matrix_ptr) :: partial_jaco
      !! partial derivatives wrt parents
    type(vector_jacobian_matrix_ptr) :: partial_jaco_p
      !! partial preconditioner derivatives wrt parents
    type(vector_jacobian_matrix_ptr) :: partial_jaco_t
      !! partial derivatives wrt parents_t (only used when status == STATUS_RES)

    type(jacobian_matrix_ptr), allocatable :: total_jaco(:)
      !! total derivatives wrt main vars
    type(jacobian_matrix_ptr), allocatable :: total_jaco_p(:)
      !! total preconditioner derivatives wrt main vars
    type(jacobian_matrix_ptr), allocatable :: total_jaco_t(:)
      !! total derivatives wrt time derivatives of main vars (STATUS_RES only)

    type(vector_jacobian_chain_ptr), allocatable :: jchain(:)
      !! jacobian chains for computation of total derivatives wrt main vars
    type(vector_jacobian_chain_ptr), allocatable :: jchain_p(:)
      !! jacobian chains for computation of total precondtioner derivatives wrt main vars
    type(vector_jacobian_chain_ptr), allocatable :: jchain_t(:)
      !! jacobian chains for computation of total derivatives wrt time der. of main vars
  contains
    procedure :: init     => node_init
    procedure :: destruct => node_destruct
    procedure :: analyze  => node_analyze
    procedure :: eval     => node_eval
  end type

  type node_ptr
    type(node), pointer :: p => null()
  end type

#define T node_ptr
#define TT type(node_ptr)
#include "../util/vector_def.f90.inc"

  type depgraph_equ
    !! equation with additional information
    integer                  :: id
      !! index for depgraph%equs
    class(equation), pointer :: e => null()
      !! pointer to corresponding equation
    integer                  :: imain
      !! index of main variable node (only used for residual equations)
    integer                  :: ires
      !! index of residual variable node (only used for residual equations)
    integer, allocatable     :: iprov(:)
      !! indices of provided variable nodes
    logical                  :: evaluated
      !! evaluation flag for traversal
  contains
    procedure :: init => depgraph_equ_init
  end type

#define T depgraph_equ
#define TT type(depgraph_equ)
#include "../util/vector_def.f90.inc"

  type depgraph
    !! equation system dependency graph
    type(vector_depgraph_equ) :: equs
      !! equations
    type(vector_node_ptr)     :: nodes
      !! all nodes
    type(hashmap_int)         :: hnodes
      !! hashmap for vselector => node index
    type(vector_int)          :: ires
      !! residual equation indices
    type(vector_int)          :: imvar
      !! main variable indices
    type(vector_int)          :: ieval
      !! equation evaluation list
    logical                   :: prec
      !! preconditioner flag
  contains
    procedure :: init     => depgraph_init
    procedure :: destruct => depgraph_destruct
    procedure :: add_equ  => depgraph_add_equ
    procedure :: add_node => depgraph_add_node
    procedure :: analyze  => depgraph_analyze
    procedure :: output   => depgraph_output
  end type

contains

#define T node_ptr
#define TT type(node_ptr)
#include "../util/vector_imp.f90.inc"

#define T depgraph_equ
#define TT type(depgraph_equ)
#include "../util/vector_imp.f90.inc"

  recursive subroutine node_init(this, g, id, e, v, iprov, iimvar, status)
    !! initialize dependency graph node
    class(node),             intent(out)   :: this
    type(depgraph),          intent(inout) :: g
      !! dependency graph
    integer,                 intent(in)    :: id
      !! index this node is saved at in g%nodes
    type(depgraph_equ),      intent(in)    :: e
      !! corresponding dependency graph equation
    type(vselector), target, intent(in)    :: v
      !! corresponding vselector
    integer,                 intent(in)    :: iprov
      !! index for e%e%vprov (only used when status == STATUS_PROV)
    integer,                 intent(in)    :: iimvar
      !! index for g%imvar (only used when status == STATUS_MAIN)
    integer,                 intent(in)    :: status
      !! node status (cf. STATUS_DEP/PROV/MAIN/RES)

    integer, parameter :: CAP = 8
    integer :: idep

    ! set simple members
    this%id       =  id
    this%status   =  status
    this%visited  =  (status == STATUS_MAIN) ! main vars are implicitly visited
    this%analyzed =  (status == STATUS_MAIN) ! main vars are implicitly analyzed
    this%const    =  .false.
    this%v        => v

    ! main var index
    if (status == STATUS_MAIN) then
      this%iimvar = iimvar
    else
      this%iimvar = 0
    end if

    ! set rest of members only for non-dependency node
    if (status == STATUS_DEP) return

    ! save depgraph equation pointer
    this%iequ = e%id

    ! init vectors
    call this%parents%init(       0, c = CAP)
    call this%parents_t%init(     0, c = CAP)
    call this%partial_jaco%init(  0, c = CAP)
    call this%partial_jaco_t%init(0, c = CAP)
    call this%partial_jaco_p%init(0, c = CAP)

    ! get parents, parents_t
    if      (status == STATUS_PROV) then
      ! collect dependencies
      do idep = 1, e%e%vdep%n
        if (associated(e%e%jaco(iprov,idep)%p)) then
          ! add node as dependency to graph and save as parent
          call this%parents%push(g%add_node(e, e%e%vdep%d(idep)%p, 0, STATUS_DEP))

          ! save jacobian matrix pointer for partial derivatives
          call this%partial_jaco%push(jacobian_matrix_ptr(e%e%jaco(iprov,idep)%p%matr))

          ! if precondioner should be computed: use precondionter if avail, otherwise use exact jacobian
          if (e%e%has_precon(iprov, idep)) then
            call this%partial_jaco_p%push(jacobian_matrix_ptr(e%e%jaco_p(iprov,idep)%p%matr))
          else
            call this%partial_jaco_p%push(jacobian_matrix_ptr(e%e%jaco(  iprov,idep)%p%matr))
          end if
        end if
      end do

    else if (status == STATUS_RES) then
      ! cast equation to residual equation
      select type (r => e%e)
        class is (res_equation)
          do idep = 1, r%vdep%n
            if (associated(r%jaco_f(idep)%p)) then
              ! add node as dependency to graph and save as parent
              call this%parents%push(g%add_node(e, r%vdep%d(idep)%p, 0, STATUS_DEP))

              ! save jacobian matrix pointer for partial derivatives
              call this%partial_jaco%push(jacobian_matrix_ptr(r%jaco_f(idep)%p%matr))

              ! if precondioner should be computed: use precondionter if avail, otherwise use exact jacobian
              if (r%has_precon_f(idep)) then
                call this%partial_jaco_p%push(jacobian_matrix_ptr(r%jaco_fp(idep)%p%matr))
              else
                call this%partial_jaco_p%push(jacobian_matrix_ptr(r%jaco_f( idep)%p%matr))
              end if
            end if
            if (associated(r%jaco_ft(idep)%p)) then
              ! add node as dependency to graph and save as parent
              call this%parents_t%push(g%add_node(e, r%vdep%d(idep)%p, 0, STATUS_DEP))

              ! save jacobian matrix pointer for partial derivatives
              call this%partial_jaco_t%push(jacobian_matrix_ptr(r%jaco_ft(idep)%p%matr))
            end if
          end do

        class default
          call program_error("status = STATUS_RES but equation is not a residual equation")

      end select
    end if
  end subroutine

  subroutine node_destruct(this)
    !! free memory
    class(node), intent(inout) :: this

    integer :: i, j

    if (allocated(this%jchain)) then
      do i = 1, size(this%jchain)
        do j = 1, this%jchain(i)%n
          if (associated(this%jchain(i)%d(j)%p)) deallocate (this%jchain(i)%d(j)%p)
        end do
      end do
    end if

    if (allocated(this%jchain_p)) then
      do i = 1, size(this%jchain_p)
        do j = 1, this%jchain_p(i)%n
          if (associated(this%jchain_p(i)%d(j)%p)) deallocate (this%jchain_p(i)%d(j)%p)
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

  recursive subroutine node_analyze(this, g)
    !! initialize jacobian chains and total derivatives
    class(node),    intent(inout) :: this
    type(depgraph), intent(inout) :: g

    integer, parameter :: CAP = 8
    integer :: i

    ! check if this node was already visited
    if (this%visited) then
      call this%v%print()
      call program_error("circular dependency detected")
    else
      this%visited = .true.
    end if

    ! set constant flag, reset later if any parent is not constant (main vars are not constant and already analyzed)
    this%const = .true.

    ! analyze parents which are not already analyzed
    do i = 1, this%parents%n
      associate (parent => g%nodes%d(this%parents%d(i))%p)
        if (.not. parent%analyzed) call parent%analyze(g)

        ! reset const flag if parent is not constant
        if (.not. parent%const) this%const = .false.
      end associate
    end do

    ! analyze parents_t which are not already analyzed
    do i = 1, this%parents_t%n
      associate (parent => g%nodes%d(this%parents_t%d(i))%p)
        if (.not. parent%analyzed) call parent%analyze(g)

        ! reset const flag if parent is not constant
        if (.not. parent%const) this%const = .false.
      end associate
    end do

    ! init jacobian chains
    if (this%const) then
      ! evaluate provider equation if constant node (no derivatives needed, just provided values)
      call g%equs%d(this%iequ)%e%eval()
    else
      ! allocate total derivative jacobian matrices and chain vectors
      allocate (this%total_jaco(  g%imvar%n), &
        &       this%total_jaco_p(g%imvar%n), &
        &       this%total_jaco_t(g%imvar%n), &
        &       this%jchain(      g%imvar%n), &
        &       this%jchain_p(    g%imvar%n), &
        &       this%jchain_t(    g%imvar%n)  )

      ! init jacobian chain vectors
      do i = 1, g%imvar%n
        call this%jchain(  i)%init(0, c=CAP)
        call this%jchain_t(i)%init(0, c=CAP)
        call this%jchain_p(i)%init(0, c=CAP)
      end do

      ! add jacobian chains for parents
      call init_jchains(            this%parents,   this%partial_jaco,   this%total_jaco,   this%jchain,   .false., .false.)
      call init_jchains(            this%parents_t, this%partial_jaco_t, this%total_jaco_t, this%jchain_t, .true.,  .false.)
      if (g%prec) call init_jchains(this%parents,   this%partial_jaco_p, this%total_jaco_p, this%jchain_p, .false., .true. )
    end if

    ! set analyzed flag
    this%analyzed = .true.

    associate (e => g%equs%d(this%iequ))
      ! analyze other nodes of equation
      if (e%ires > 0) then
        if (.not. g%nodes%d(e%ires)%p%analyzed) call g%nodes%d(e%ires)%p%analyze(g)
      end if
      do i = 1, size(e%iprov)
        if (.not. g%nodes%d(e%iprov(i))%p%analyzed) call g%nodes%d(e%iprov(i))%p%analyze(g)
      end do

      ! add equation to evaluation list
      if (.not. e%evaluated) then
        call g%ieval%push(this%iequ)
        e%evaluated = .true.
      end if
    end associate

  contains

    subroutine init_jchains(parents, partial_jaco, total_jaco, jchain, time, prec)
      type(vector_int),                 intent(in)    :: parents
        !! parents or parents_t
      type(vector_jacobian_matrix_ptr), intent(in)    :: partial_jaco
        !! partial_jaco or partial_jaco_t
      type(jacobian_matrix_ptr),        intent(inout) :: total_jaco(:)
        !! total_jaco or total_jaco_t
      type(vector_jacobian_chain_ptr),  intent(inout) :: jchain(:)
        !! jchain or jchain_t
      logical,                          intent(in)    :: time
        !! time derivative flag
      logical,                          intent(in)    :: prec
        !! preconditioner flag

      integer                                :: i, j, k
      type(jacobian_matrix),     pointer     :: n_total_jaco
      type(jacobian_matrix_ptr)              :: jmul(2)
      type(jacobian_matrix_ptr), allocatable :: jadd(:)
      type(jacobian_chain_ptr)               :: chain

      ! add multiplication chains for parents
      do i = 1, parents%n
        associate (n => g%nodes%d(parents%d(i))%p)
          if      (n%const) then
            ! do nothing for const node
            cycle
          else if (n%status == STATUS_MAIN) then
            ! set total derivatives
            total_jaco(n%iimvar)%p => partial_jaco%d(i)%p
          else
            ! chain rule for all main vars
            do j = 1, g%imvar%n
              if      (time) then
                n_total_jaco => n%total_jaco_t(j)%p
              else if (prec) then
                if (associated(n%total_jaco_p(j)%p)) then
                  n_total_jaco => n%total_jaco_p(j)%p
                else
                  n_total_jaco => n%total_jaco(j)%p
                end if
              else
                n_total_jaco => n%total_jaco(j)%p
              end if

              ! do nothing if no dependency on this main var
              if (.not. associated(n_total_jaco)) cycle

              ! new multiplication chain
              nullify (chain%p)
              allocate (jacobian_mul_chain :: chain%p)

              ! init multiplication chain
              jmul(1)%p => partial_jaco%d(i)%p
              jmul(2)%p => n_total_jaco
              call chain%p%init(jmul)

              ! save chain
              call jchain(j)%push(chain)
            end do
          endif
        end associate
      end do

      ! add addition chains for parents (only after all multiplication chains have been added)
      do i = 1, g%imvar%n
        if      (jchain(i)%n == 0) then ! no chain for this main var
          cycle
        else if ((.not. associated(total_jaco(i)%p)) .and. (jchain(i)%n == 1)) then ! only one chain, set total derivatives
          total_jaco(i)%p => jchain(i)%d(1)%p%result
        else ! more than one chain, add results
          ! allocate temporary jacobian matrix pointers for addition
          if (associated(total_jaco(i)%p)) then
            ! partial derivative wrt main var exists
            allocate (jadd(1 + jchain(i)%n))
            jadd(1)%p => total_jaco(i)%p
            k = 1
          else
            ! no partial derivative wrt main var exists
            allocate (jadd(jchain(i)%n))
            k = 0
          end if

          ! save rest of jacobian matrix pointers
          do j = 1, jchain(i)%n
            jadd(k + j)%p => jchain(i)%d(j)%p%result
          end do

          ! init addition chain
          nullify (chain%p)
          allocate (jacobian_add_chain :: chain%p)
          call chain%p%init(jadd)

          ! deallocate temporary jacobian matrix pointers
          deallocate (jadd)

          ! save chain
          call jchain(i)%push(chain)

          ! set/update total derivative jacobian matrix
          total_jaco(i)%p => jchain(i)%d(jchain(i)%n)%p%result
        end if
      end do
    end subroutine

  end subroutine

  subroutine node_eval(this)
    !! evaluate non-constant parts of jacobian chains
    class(node), intent(inout) :: this ! DAG node

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

    if (allocated(this%jchain_p)) then
      do i = 1, size(this%jchain_p)
        ! evaluate non-constant parts of jacobian chains
        do j = 1, this%jchain_p(i)%n
          call this%jchain_p(i)%d(j)%p%eval(const = .false.)
        end do
      end do
    end if
  end subroutine

  subroutine depgraph_equ_init(this, g, id, e)
    !! initialize dependency graph equation
    class(depgraph_equ),     intent(out)   :: this
    type(depgraph),          intent(inout) :: g
      !! dependency graph
    integer,                 intent(in)    :: id
      !! index this depgraph_equ is saved at in g%equs
    class(equation), target, intent(in)    :: e
      !! underlying equation

    integer :: iprov, iprov_f

    ! set simple members
    this%id        =  id
    this%e         => e
    this%evaluated =  .false.

    ! main variable + residuals
    iprov_f = -1
    select type (e)
      class is (res_equation)
        this%imain = g%add_node(this, e%mvar, 0, STATUS_MAIN)
        this%ires  = g%add_node(this, e%f,    0, STATUS_RES )
        iprov_f    = e%iprov_f
      class default
        this%imain = -1
        this%ires  = -1
    end select

    ! provide variable selectors (dependent variable selectors are automatically added by node_init)
    allocate (this%iprov(e%vprov%n))
    do iprov = 1, e%vprov%n
      if (iprov == iprov_f) then
        this%iprov(iprov) = this%ires ! already added
      else
        this%iprov(iprov) = g%add_node(this, e%vprov%d(iprov)%p, iprov, STATUS_PROV)
      end if
    end do
  end subroutine

  subroutine depgraph_init(this, prec)
    !! initialize dependency graph
    class(depgraph), intent(out) :: this
    logical,         intent(in)  :: prec
      !! should a preconditioner be computed?

    integer, parameter :: CAP = 32

    this%prec = prec

    call this%equs%init( 0, c = CAP)
    call this%nodes%init(0, c = CAP)
    call this%hnodes%init(  c = CAP)
    call this%ires%init( 0, c = CAP)
    call this%imvar%init(0, c = CAP)
    call this%ieval%init(0, c = CAP)
  end subroutine

  subroutine depgraph_destruct(this)
    !! release memory
    class(depgraph), intent(inout) :: this

    integer :: i

    call this%equs%destruct()
    do i = 1, this%nodes%n
      call this%nodes%d(i)%p%destruct()
      deallocate (this%nodes%d(i)%p)
    end do
    call this%nodes%destruct()
    call this%ires%destruct()
    call this%imvar%destruct()
    call this%ieval%destruct()
  end subroutine

  subroutine depgraph_add_equ(this, e)
    !! add new equation to dependency graph
    class(depgraph), intent(inout) :: this
    class(equation), intent(in)    :: e
      !! equation to add

    type(depgraph_equ) :: tmp

    ! init new depgraph equation
    call this%equs%push(tmp)
    call this%equs%d(this%equs%n)%init(this, this%equs%n, e)

    ! check for residual equation
    select type (e)
      class is (res_equation)
        call this%ires%push(this%equs%n)
    end select
  end subroutine

  recursive function depgraph_add_node(this, e, v, iprov, status) result(inode)
    !! add new node to dependency graph
    class(depgraph),    intent(inout) :: this
    type(depgraph_equ), intent(in)    :: e
      !! dependency graph equation
    type(vselector),    intent(in)    :: v
      !! variable selector to add
    integer,            intent(in)    :: iprov
      !! index in e%e%vprov (only used for STATUS_PROV)
    integer,            intent(in)    :: status
      !! node status
    integer                           :: inode
      !! return node index

    integer :: hkey(v%hashkey_size())
    logical :: hstat

    ! get hashmap key
    hkey = v%hashkey()

    ! check if vselector already exists
    call this%hnodes%get(hkey, inode, hstat)
    if (hstat) then
      ! node was found
      associate (m => this%nodes%d(inode)%p)
        ! do nothing if dependency
        if (status == STATUS_DEP) return

        ! check if provided multiple times
        if (m%status /= STATUS_DEP) then
          print *, "vselector:"
          call v%print()
          print *, "equations:"
          print *, e%e%name
          print *, this%equs%d(m%iequ)%e%name
          call program_error("vselector provided multiple times")
        end if
      end associate
    else
      ! node not found, create a new one
      block
        type(node_ptr) :: tmp
        allocate (tmp%p)
        call this%nodes%push(tmp)
      end block
      inode = this%nodes%n
      call this%hnodes%set(hkey, inode)
    end if

    ! init/update node
    if (status == STATUS_MAIN) then
      call this%imvar%push(inode)
      call this%nodes%d(inode)%p%init(this, inode, e, v, 0,     this%imvar%n, status)
    else
      call this%nodes%d(inode)%p%init(this, inode, e, v, iprov, 0,            status)
    end if
  end function

  subroutine depgraph_analyze(this)
    !! analyze dependency graph: get evaluation list etc.
    class(depgraph), intent(inout) :: this

    integer :: i

    ! analyze nodes, start with residuals
    do i = 1, this%ires%n
      associate (e => this%equs%d(this%ires%d(i)))
        associate (n => this%nodes%d(e%ires)%p)
          if (.not. n%analyzed) call n%analyze(this)
        end associate
      end associate
    end do
  end subroutine

  subroutine depgraph_output(this, fname, pdf, png, svg)
    !! create output of dependency graph.
    !!
    !! always outputs graphviz file: `fname//'.gv'`
    !! pdf/png/svg files are optionally created.
    !! create other formats by: `$ dot -T<EXT> fname.gv > fname.<EXT>`
    !!
    !! graphviz examples: https://graphviz.org/gallery/
    class(depgraph),   intent(in) :: this
    character(*),      intent(in) :: fname
      !! file name without extension, e.g. 'output/depgraph'
    logical, optional, intent(in) :: pdf
      !! create pdf file (default: true)
    logical, optional, intent(in) :: png
      !! create png file (default: false)
    logical, optional, intent(in) :: svg
      !! create svg file (default: false)

    character(len=*), parameter :: COLOR(0:3) = ['red  ', 'black', 'blue ', 'green']
    integer :: iounit, i, j
    logical :: pdf_, png_, svg_

    ! optional args
    pdf_ = .true.
    if (present(pdf)) pdf_ = pdf
    png_ = .false.
    if (present(png)) png_ = png
    svg_ = .false.
    if (present(svg)) svg_ = svg

    ! create file + write header
    open (newunit=iounit, file=fname//'.gv', action='WRITE')
    write (iounit, '(A)') 'digraph mygraph {'
    write (iounit, '(A)') '  rankdir = "LR"'
    write (iounit, '(A)') '  node [shape=box]'

    ! write nodes
    do i = 1, this%nodes%n
      associate (n => this%nodes%d(i)%p)
        write (iounit, '(A,I0,A)', advance='no') '  ', i, ' ['
        write (iounit, '(2A)', advance='no') 'color=', trim(COLOR(n%status))
        write (iounit, '(2A)', advance='no') ' label="', n%v%name
        if (n%iequ > 0) then
          write (iounit, '(2A)', advance='no') '\n', this%equs%d(n%iequ)%e%name
        end if
        write (iounit, '(A)', advance='no') '"'
        write (iounit, '(A)', advance='no') ' penwidth=2.0'
        if (n%const) then
          write (iounit, '(A)', advance='no') ' style=filled'
          write (iounit, '(A)', advance='no') ' color=gray'
        end if
        write (iounit, '(A)') ']'
      end associate
    end do

    ! write edges
    do i = 1, this%nodes%n
      associate (n => this%nodes%d(i)%p)
        do j = 1, n%parents%n
          write (iounit, '(A,I0,A,I0)') '  ', n%parents%d(j), ' -> ', i
        end do
      end associate
    end do

    ! finish
    write (iounit, '(A)') '}'
    close (unit=iounit)

    ! create pdf/png/svg files
    if (pdf_) call execute_command_line('dot -Tpdf '//fname//'.gv > '//fname//'.pdf')
    if (png_) call execute_command_line('dot -Tpng '//fname//'.gv > '//fname//'.png')
    if (svg_) call execute_command_line('dot -Tsvg '//fname//'.gv > '//fname//'.svg')
  end subroutine

end module
