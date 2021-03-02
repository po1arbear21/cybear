#include "../util/macro.f90.inc"

module res_equation_m

  use equation_m,   only: equation, equation_realloc_jaco, equation_set_jaco_matr, equation_reset, equation_destruct
  use error_m
  use grid_table_m, only: grid_table, grid_table_ptr
  use jacobian_m,   only: jacobian, jacobian_ptr
  use stencil_m,    only: stencil_ptr
  use variable_m,   only: variable, variable_ptr
  use vselector_m,  only: vselector

  implicit none

  private
  public res_equation, res_equation_ptr, vector_res_equation_ptr

  type, abstract, extends(equation) :: res_equation
    !! residual equation (is solved, not eliminated)

    type(vselector), pointer :: mvar => null()
      !! main var selector
    logical                  :: mvar_alc = .false.
      !! has mvar been allocated in this object?

    type(variable_ptr), allocatable :: fvar(:)
      !! residual data. size: mvar%nval (mvar supplied at init).
    type(vselector)                 :: f
      !! residual var selector

    type(jacobian_ptr), allocatable :: jaco_f(:)
      !! derivatives of f wrt vdep. size: this%vdep%d
    type(jacobian_ptr), allocatable :: jaco_ft(:)
      !! derivatives of f wrt d/dt(vdep), must be const. size: this%vdep%d


  contains
    procedure :: destruct      => res_equation_destruct
    procedure :: reset         => res_equation_reset
    procedure :: realloc_jaco  => res_equation_realloc_jaco
    procedure :: init_jaco_f   => res_equation_init_jaco_f
    procedure :: set_jaco_matr => res_equation_set_jaco_matr

    procedure, private :: res_equation_init_f_vsel
    procedure, private :: res_equation_init_f_var_ntab
    procedure, private :: res_equation_init_f_var
    generic            :: init_f => res_equation_init_f_vsel, res_equation_init_f_var_ntab, res_equation_init_f_var
  end type

  type res_equation_ptr
    class(res_equation), pointer :: p => null()
  end type

#define T res_equation_ptr
#define TT type(res_equation_ptr)
#include "../util/vector_def.f90.inc"

contains

#define T res_equation_ptr
#define TT type(res_equation_ptr)
#include "../util/vector_imp.f90.inc"

  subroutine res_equation_init_f_vsel(this, mvar)
    !! set main variable and initialize residual data
    class(res_equation),      intent(inout) :: this
    class(vselector), target, intent(in)    :: mvar
      !! main variable

    integer :: ival

    ! save mvar
    this%mvar => mvar

    ! allocate residual data
    allocate (this%fvar(mvar%nval))
    do ival = 1, mvar%nval
      allocate (this%fvar(ival)%p, source = mvar%v(ival)%p)
      this%fvar(ival)%p%name = this%fvar(ival)%p%name//"_res"
    end do

    ! init residual var selector based on mvar
    call this%f%init(this%fvar, mvar%tab, mvar%name//"_res")

    ! allocate f jacobians
    allocate (this%jaco_f( size(this%vdep%d)))
    allocate (this%jaco_ft(size(this%vdep%d)))
  end subroutine

  subroutine res_equation_init_f_var_ntab(this, mvar, tab, name)
    !! set main variable and initialize residual data
    class(res_equation),     intent(inout) :: this
    class(variable), target, intent(in)    :: mvar
      !! main variable
    type(grid_table_ptr),    intent(in)    :: tab(:)
      !! grid table pointers
    character(*), optional,  intent(in)    :: name
      !! name of new var selector (default: var%name)

    character(:), allocatable :: name_
    type(vselector), pointer  :: vsel

    ! optional argument
    if (present(name)) then
      allocate (name_, source = name)
    else
      allocate (name_, source = mvar%name)
    end if

    ! create vselector from variable and keep track of memory
    allocate (vsel)
    call vsel%init(mvar, tab, name=name_)
    this%mvar_alc = .true.

    ! init by vselector
    call this%init_f(vsel)
  end subroutine

  subroutine res_equation_init_f_var(this, mvar, tab, name)
    !! set main variable and initialize residual data
    class(res_equation),     intent(inout) :: this
    class(variable), target, intent(in)    :: mvar
      !! main variable
    type(grid_table),        intent(in)    :: tab
      !! grid table
    character(*), optional,  intent(in)    :: name
      !! name of new var selector (default: var%name)

    call this%init_f(mvar, [tab%get_ptr()], name=name)
  end subroutine

  subroutine res_equation_destruct(this)
    !! destruct residual equation
    class(res_equation), intent(inout) :: this

    integer :: i

    if (this%mvar_alc) then
      if (associated(this%mvar)) deallocate (this%mvar)
    else
      nullify (this%mvar)
    end if

    do i = 1, size(this%fvar)
      if (associated(this%fvar(i)%p)) deallocate (this%fvar(i)%p)
    end do
    deallocate (this%fvar)
    do i = 1, size(this%jaco_f)
      if (associated(this%jaco_f(i)%p)) then
        call this%jaco_f(i)%p%destruct()
        deallocate (this%jaco_f(i)%p)
      end if
    end do
    deallocate (this%jaco_f)
    do i = 1, size(this%jaco_ft)
      if (associated(this%jaco_ft(i)%p)) then
        call this%jaco_ft(i)%p%destruct()
        deallocate (this%jaco_ft(i)%p)
      end if
    end do
    deallocate (this%jaco_ft)

    ! destruct base
    call equation_destruct(this)
  end subroutine

  subroutine res_equation_reset(this)
    !! reset provided vars, f and non-const jacobians
    class(res_equation), intent(inout) :: this

    integer :: i

    ! reset vprov and this%jaco
    call equation_reset(this)

    ! reset f
    call this%f%reset()

    ! reset non-const parts of jaco_f
    do i = 1, size(this%jaco_f)
      if (.not. associated(this%jaco_f(i)%p)) cycle
      call this%jaco_f(i)%p%reset(const = .false.)
    end do
  end subroutine

  subroutine res_equation_realloc_jaco(this, cprov, cdep)
    !! reallocate this%jaco, this%jaco_f and this%jaco_ft if initial capacity was not big enough
    class(res_equation), intent(inout) :: this
    integer,             intent(in)    :: cprov
      !! new prov capacity
    integer,             intent(in)    :: cdep
      !! new dep capacity

    type(jacobian_ptr), allocatable :: jaco_f_tmp(:), jaco_ft_tmp(:)

    ! call base (reallocate this%jaco)
    call equation_realloc_jaco(this, cprov, cdep)

    ! reallocate this%jaco_f and this%jaco_ft
    if (cdep > this%vdep%n) then
      allocate (jaco_f_tmp( cdep))
      allocate (jaco_ft_tmp(cdep))
      jaco_f_tmp( 1:this%vdep%n) = this%jaco_f
      jaco_ft_tmp(1:this%vdep%n) = this%jaco_ft
      call move_alloc(jaco_f_tmp,  this%jaco_f )
      call move_alloc(jaco_ft_tmp, this%jaco_ft)
    end if
  end subroutine

  function res_equation_init_jaco_f(this, idep, st, const, zero, valmsk, dtime) result(jaco)
    !! allocate and initialize f/ft jacobian
    class(res_equation), intent(inout) :: this
    integer,             intent(in)    :: idep
      !! dependency var index
    type(stencil_ptr),   intent(in)    :: st(:)
      !! stencils (v1%ntab)
    logical, optional,   intent(in)    :: const
      !! const flag; default = false for f, true for ft; the same for all blocks
    logical, optional,   intent(in)    :: zero(:,:)
      !! zero flags (v1%ntab x v2%ntab); default: set automatically by checking stencils
    logical, optional,   intent(in)    :: valmsk(:,:)
      !! value mask (v1%nval x v2%nval); default = true; the same for all blocks
    logical, optional,   intent(in)    :: dtime
      !! false: init f; true: init ft (default: false)
    type(jacobian), pointer            :: jaco
      !! return pointer to newly created jacobian

    logical :: const_, dtime_

    dtime_ = .false.
    if (present(dtime)) dtime_ = dtime
    const_ = dtime_
    if (present(const)) const_ = const

    associate (vdep => this%vdep%d(idep)%p)
      if (dtime_) then
        ASSERT(const_)
        allocate (this%jaco_ft(idep)%p)
        call this%jaco_ft(idep)%p%init(this%f, vdep, st, const = const_, zero = zero, valmsk = valmsk)
        jaco => this%jaco_ft(idep)%p
      else
        allocate (this%jaco_f(idep)%p)
        call this%jaco_f(idep)%p%init(this%f, vdep, st, const = const_, zero = zero, valmsk = valmsk)
        jaco => this%jaco_f(idep)%p
      end if
    end associate
  end function

  subroutine res_equation_set_jaco_matr(this, const, nonconst)
    !! set jacobian matrices
    class(res_equation), intent(inout) :: this
    logical, optional,   intent(in)    :: const
      !! enable processing of const blocks (default: true)
    logical, optional,   intent(in)    :: nonconst
      !! enable processing of non-const blocks (default: true)

    integer :: i

    ! base
    call equation_set_jaco_matr(this, const = const, nonconst = nonconst)

    ! jaco_f and jaco_ft
    do i = 1, this%vdep%n
      if (associated(this%jaco_f( i)%p)) call this%jaco_f( i)%p%set_matr(const = const, nonconst = nonconst)
      if (associated(this%jaco_ft(i)%p)) call this%jaco_ft(i)%p%set_matr(const = const, nonconst = nonconst)
    end do
  end subroutine

end module
