#include "../util/macro.f90.inc"

module equation_m

  use error_m
  use grid_table_m, only: grid_table, grid_table_ptr
  use jacobian_m,   only: jacobian, jacobian_ptr
  use stencil_m,    only: stencil_ptr
  use variable_m,   only: variable
  use vector_m,     only: vector_int
  use vselector_m,  only: vselector, vselector_ptr, vector_vselector_ptr

  implicit none

  private
  public equation, equation_ptr
  public equation_realloc_jaco
  public equation_destruct
  public equation_reset
  public equation_set_jaco_matr

  type, abstract :: equation
    !! abstract equation base

    character(:), allocatable  :: name
      !! equation name

    type(vector_vselector_ptr) :: vprov
      !! provided variable selectors
    type(vector_vselector_ptr) :: vdep
      !! dependency variable selectors
    type(vector_int)           :: vprov_alc
      !! which vprov%d(i)%p were allocated in this%provide(var, tab)?
    type(vector_int)           :: vdep_alc
      !! which vdep%d(i)%p were allocated in this%depend(var, tab)?

    type(jacobian_ptr), allocatable :: jaco(:,:)
      !! derivatives of vprov wrt vdep (vprov%n x vdep%n)

    logical :: finished_init
      !! indicates whether init_final was called
  contains
    procedure :: equation_init
    procedure :: destruct      => equation_destruct
    procedure :: reset         => equation_reset
    procedure :: realloc_jaco  => equation_realloc_jaco
    procedure :: init_jaco     => equation_init_jaco
    procedure :: init_final    => equation_init_final
    procedure :: set_jaco_matr => equation_set_jaco_matr
    procedure :: test          => equation_test

    procedure, private :: provide_vselector     => equation_provide_vselector
    procedure, private :: provide_variable_ntab => equation_provide_variable_ntab
    procedure, private :: provide_variable      => equation_provide_variable
    generic            :: provide => provide_vselector, provide_variable_ntab, provide_variable

    procedure, private :: depend_vselector     => equation_depend_vselector
    procedure, private :: depend_variable_ntab => equation_depend_variable_ntab
    procedure, private :: depend_variable      => equation_depend_variable
    generic            :: depend => depend_vselector, depend_variable_ntab, depend_variable

    procedure(equation_eval), deferred :: eval
  end type

  abstract interface
    subroutine equation_eval(this)
      !! evaluate equation (implement in child classes)
      import equation
      class(equation), intent(inout) :: this
    end subroutine
  end interface

  type equation_ptr
    class(equation), pointer :: p => null()
  end type

#define T equation_ptr
#define TT type(equation_ptr)
#include "../util/vector_def.f90.inc"

contains

#define T equation_ptr
#define TT type(equation_ptr)
#include "../util/vector_imp.f90.inc"

  subroutine equation_init(this, name)
    !! initialize equation base
    class(equation), intent(out) :: this
    character(*),    intent(in)  :: name
      !! name of equation

    integer, parameter :: cap = 16

    ! save name
    this%name = name

    ! init vprov and vdep vectors
    call this%vprov%init(    0, c = cap)
    call this%vdep%init(     0, c = cap)
    call this%vprov_alc%init(0, c = cap)
    call this%vdep_alc%init( 0, c = cap)

    ! allocate jaco
    allocate (this%jaco(cap,cap))

    ! init_final has not been called yet
    this%finished_init = .false.
  end subroutine

  subroutine equation_destruct(this)
    !! destruct equation
    class(equation), intent(inout) :: this

    integer :: i, j

    if (allocated(this%name)) deallocate (this%name)

    do i = 1, this%vprov_alc%n
      deallocate (this%vprov%d(this%vprov_alc%d(i))%p)
    end do
    do i = 1, this%vdep_alc%n
      deallocate (this%vdep%d(this%vdep_alc%d(i))%p)
    end do

    call this%vprov%destruct()
    call this%vdep%destruct()
    call this%vprov_alc%destruct()
    call this%vdep_alc%destruct()

    ! destruct jacobians
    do j = 1, size(this%jaco,2); do i = 1, size(this%jaco,1)
      if (associated(this%jaco(i,j)%p)) then
        call this%jaco(i,j)%p%destruct()
        deallocate (this%jaco(i,j)%p)
      end if
    end do; end do
    deallocate (this%jaco)
  end subroutine

  subroutine equation_reset(this)
    !! reset provided var selectors and non-const jacobians
    class(equation), intent(inout) :: this

    integer :: i, j

    ! reset provided var selectors
    do i = 1, this%vprov%n
      call this%vprov%d(i)%p%reset()
    end do

    ! reset non-const parts of jacobians
    do j = 1, size(this%jaco,2); do i = 1, size(this%jaco,1)
      if (.not. associated(this%jaco(i,j)%p)) cycle
      call this%jaco(i,j)%p%reset(const = .false.)
    end do; end do
  end subroutine

  function equation_provide_vselector(this, vsel) result(iprov)
    !! provide var selector
    class(equation),         intent(inout) :: this
    type(vselector), target, intent(in)    :: vsel
      !! new provided var selector
    integer                                :: iprov
      !! return provided index

    ! reallocate jaco if necessary
    if (this%vprov%n >= size(this%vprov%d)) then
      call this%realloc_jaco((this%vprov%n + 1) * 2, size(this%vdep%d))
    end if

    ! add var selector to provided variables
    call this%vprov%push(vsel%get_ptr())

    ! return index
    iprov = this%vprov%n
  end function

  function equation_provide_variable_ntab(this, var, tab, name) result(iprov)
    !! provide variable for multiple grid tables, creates var selector internally
    class(equation),        intent(inout) :: this
    class(variable),        intent(in)    :: var
      !! new provided variable
    type(grid_table_ptr),   intent(in)    :: tab(:)
      !! grid table pointers
    character(*), optional, intent(in)    :: name
      !! name of new var selector (default: var%name)
    integer                               :: iprov
      !! return provided index

    character(:), allocatable :: name_
    type(vselector), pointer  :: vsel

    ! optional argument
    if (present(name)) then
      allocate (name_, source = name)
    else
      allocate (name_, source = var%name)
    end if

    ! create vselector from variable and keep track of memory
    allocate (vsel)
    call vsel%init(var, tab, name=name_)
    call this%vprov_alc%push(this%vprov%n)

    ! add provided vselector
    iprov = this%provide(vsel)
  end function

  function equation_provide_variable(this, var, tab, name) result(iprov)
    !! provide variable for single grid table, creates var selector internally
    class(equation),        intent(inout) :: this
    class(variable),        intent(in)    :: var
      !! new provided variable
    type(grid_table),       intent(in)    :: tab
      !! grid table pointers
    character(*), optional, intent(in)    :: name
      !! name of new var selector (default: var%name)
    integer                               :: iprov
      !! return provided index

    iprov = this%provide(var, [tab%get_ptr()], name=name)
  end function

  function equation_depend_vselector(this, vsel) result(idep)
    !! depend on var selector
    class(equation),          intent(inout) :: this
    class(vselector), target, intent(in)    :: vsel
      !! new dependency var selector
    integer                                 :: idep
      !! return dependency index

    type(vselector_ptr) :: vptr

    ! reallocate g if necessary
    if (this%vdep%n >= size(this%vdep%d)) then
      call this%realloc_jaco(size(this%vprov%d), (this%vdep%n + 1) * 2)
    end if

    ! add var selector to dependent variables
    vptr%p => vsel
    call this%vdep%push(vptr)

    ! return dependency index
    idep = this%vdep%n
  end function

  function equation_depend_variable_ntab(this, var, tab, name) result(idep)
    !! add new dependency variable, creates var selector internally
    class(equation),        intent(inout) :: this
    class(variable),        intent(in)    :: var
      !! new dependency variable
    type(grid_table_ptr),   intent(in)    :: tab(:)
      !! grid table pointers
    character(*), optional, intent(in)    :: name
      !! name of new var selector (default: var%name)
    integer                               :: idep
      !! return dependency index

    character(:), allocatable :: name_
    type(vselector), pointer  :: vsel

    ! optional argument
    if (present(name)) then
      allocate (name_, source = name)
    else
      allocate (name_, source = var%name)
    end if

    ! create vselector from variable and keep track of memory
    allocate (vsel)
    call vsel%init(var, tab, name=name_)
    call this%vdep_alc%push(this%vdep%n)

    ! add dependency vselector
    idep = this%depend(vsel)
  end function

  function equation_depend_variable(this, var, tab, name) result(idep)
    !! add new dependency variable for single grid table, creates var selector internally
    class(equation),        intent(inout) :: this
    class(variable),        intent(in)    :: var
      !! new dependency variable
    type(grid_table),       intent(in)    :: tab
      !! grid table pointers
    character(*), optional, intent(in)    :: name
      !! name of new var selector (default: var%name)
    integer                               :: idep
      !! return dependency index

    idep = this%depend(var, [tab%get_ptr()], name=name)
  end function

  subroutine equation_realloc_jaco(this, cprov, cdep)
    !! reallocate this%jaco, if initial capacity was not big enough
    class(equation), intent(inout) :: this
    integer,         intent(in)    :: cprov
      !! new prov capacity
    integer,         intent(in)    :: cdep
      !! new dep capacity

    type(jacobian_ptr), allocatable :: jaco_tmp(:,:)

    ! reallocate this%jaco
    allocate (jaco_tmp(cprov, cdep))
    jaco_tmp(1:this%vprov%n,1:this%vdep%n) = this%jaco(1:this%vprov%n,1:this%vdep%n)
    call move_alloc(jaco_tmp, this%jaco)
  end subroutine

  function equation_init_jaco(this, iprov, idep, st, const, zero, valmsk) result(jaco)
    !! allocate and initialize jacobian
    class(equation),   intent(inout) :: this
    integer,           intent(in)    :: iprov
      !! provided var selector index
    integer,           intent(in)    :: idep
      !! dependency var selector index
    type(stencil_ptr), intent(in)    :: st(:)
      !! stencils (v1%ntab)
    logical, optional, intent(in)    :: const
      !! const flag; default = false; the same for all blocks
    logical, optional, intent(in)    :: zero(:,:)
      !! zero flags (v1%ntab x v2%ntab); default: set automatically by checking stencils
    logical, optional, intent(in)    :: valmsk(:,:)
      !! value mask (v1%nval x v2%nval); default = true; the same for all blocks
    type(jacobian), pointer          :: jaco
      !! return pointer to newly created jacobian

    associate (vprov => this%vprov%d(iprov)%p, vdep => this%vdep%d(idep)%p)
      ! allocate and init jacobian
      allocate (this%jaco(iprov,idep)%p)
      call this%jaco(iprov,idep)%p%init(vprov, vdep, st, const = const, zero = zero, valmsk = valmsk)

      ! return pointer to jacobian
      jaco => this%jaco(iprov,idep)%p
    end associate
  end function

  subroutine equation_init_final(this)
    !! set constant parts of jacobian matrices
    class(equation), intent(inout) :: this

    if (this%finished_init) call program_error("init_final called multiple times")
    this%finished_init = .true.

    call this%set_jaco_matr(const = .true., nonconst = .false.)
  end subroutine

  subroutine equation_set_jaco_matr(this, const, nonconst)
    !! set jacobian matrices
    class(equation),   intent(inout) :: this
    logical, optional, intent(in)    :: const
      !! enable processing of const blocks (default: true)
    logical, optional, intent(in)    :: nonconst
      !! enable processing of non-const blocks (default: true)

    integer :: i, j

    do i = 1, this%vprov%n
      do j = 1, this%vdep%n
        if (associated(this%jaco(i,j)%p)) call this%jaco(i,j)%p%set_matr(const = const, nonconst = nonconst)
      end do
    end do
  end subroutine

  subroutine equation_test(this)
    !! test jacobians with finite differences
    class(equation), intent(inout) :: this

    ! FIXME
  end subroutine

end module
