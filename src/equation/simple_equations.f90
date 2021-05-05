#include "../util/macro.f90.inc"

module simple_equations_m

  use equation_m,     only: equation
  use grid_m,         only: grid_table, grid_table_ptr
  use jacobian_m,     only: jacobian, jacobian_ptr
  use qsort_m,        only: qsort
  use res_equation_m, only: res_equation
  use stencil_m,      only: stencil_ptr, dirichlet_stencil, empty_stencil
  use variable_m,     only: variable, variable_ptr
  use vselector_m,    only: vselector, vselector_ptr

  implicit none

  private
  public dummy_equation, selector_equation, input_equation

  type, extends(equation) :: dummy_equation
    !! dummy equation which provides values but does not change anything
  contains
    generic   :: init => dummy_equation_init_vsel,      &
      &                  dummy_equation_init_nvar_ntab, &
      &                  dummy_equation_init_var_ntab,  &
      &                  dummy_equation_init_nvar_tab,  &
      &                  dummy_equation_init_var_tab
    procedure :: eval => dummy_equation_eval

    procedure, private :: dummy_equation_init_vsel,      &
      &                   dummy_equation_init_nvar_ntab, &
      &                   dummy_equation_init_var_ntab,  &
      &                   dummy_equation_init_nvar_tab,  &
      &                   dummy_equation_init_var_tab
  end type

  type, extends(equation) :: selector_equation
    !! equation that selects variables from one or multiple var selectors

    type(dirichlet_stencil) :: dir_st
    type(empty_stencil)     :: emp_st
  contains
    procedure :: init => selector_equation_init
    procedure :: eval => selector_equation_eval
  end type

  type, extends(res_equation) :: input_equation
    !! provide value as input parameter

    real, allocatable :: appl(:)
      !! applied values

    type(dirichlet_stencil)  :: st
  contains
    generic   :: init  => input_equation_init_vsel,      &
      &                   input_equation_init_nvar_ntab, &
      &                   input_equation_init_var_ntab,  &
      &                   input_equation_init_nvar_tab,  &
      &                   input_equation_init_var_tab
    procedure :: apply => input_equation_apply
    procedure :: eval  => input_equation_eval

    procedure, private :: input_equation_init_vsel,      &
      &                   input_equation_init_nvar_ntab, &
      &                   input_equation_init_var_ntab,  &
      &                   input_equation_init_nvar_tab,  &
      &                   input_equation_init_var_tab,   &
      &                   input_equation_init_body
  end type

contains

  subroutine dummy_equation_init_vsel(this, v)
    !! initialize dummy equation
    class(dummy_equation), intent(out) :: this
    class(vselector),      intent(in)  :: v
      !! dummy variable

    integer :: iprov

    ! init base
    call this%equation_init("Provide"//v%name)

    ! add provided var
    iprov = this%provide(v)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine dummy_equation_init_nvar_ntab(this, v, tab, name)
    !! initialize dummy equation
    class(dummy_equation), intent(out) :: this
    type(variable_ptr),    intent(in)  :: v(:)
      !! variable pointers
    type(grid_table_ptr),  intent(in)  :: tab(:)
      !! grid table pointers
    character(*),          intent(in)  :: name
      !! selector name

    integer :: iprov

    ! init base
    call this%equation_init("Provide"//name)

    ! add provided var
    iprov = this%provide(v, tab, name)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine dummy_equation_init_var_ntab(this, v, tab, name)
    !! initialize dummy equation
    class(dummy_equation),  intent(out) :: this
    class(variable),        intent(in)  :: v
      !! provided variable
    type(grid_table_ptr),   intent(in)  :: tab(:)
      !! grid table pointers
    character(*), optional, intent(in)  :: name
      !! name of new var selector (default: v%name)

    integer :: iprov

    ! init base
    if (present(name)) then
      call this%equation_init("Provide"//name)
    else
      call this%equation_init("Provide"//v%name)
    end if

    ! add provided var
    iprov = this%provide(v, tab, name=name)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine dummy_equation_init_nvar_tab(this, v, name, tab)
    !! initialize dummy equation
    class(dummy_equation),      intent(out) :: this
    type(variable_ptr),         intent(in)  :: v(:)
      !! variable pointers
    character(*),               intent(in)  :: name
      !! selector name
    type(grid_table), optional, intent(in)  :: tab
      !! grid table (default: variables' whole grids via v%g%tab_all)

    integer :: iprov

    ! init base
    call this%equation_init("Provide"//name)

    ! add provided var
    iprov = this%provide(v, name, tab=tab)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine dummy_equation_init_var_tab(this, v, tab, name)
    !! initialize dummy equation
    class(dummy_equation),      intent(out) :: this
    class(variable),            intent(in)  :: v
      !! new provided variable
    type(grid_table), optional, intent(in)  :: tab
      !! grid table (default: variable's whole grid via v%g%tab_all)
    character(*),     optional, intent(in)  :: name
      !! name of new var selector (default: var%name)

    integer :: iprov

    ! init base
    if (present(name)) then
      call this%equation_init("Provide"//name)
    else
      call this%equation_init("Provide"//v%name)
    end if

    ! add provided var
    iprov = this%provide(v, tab=tab, name=name)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine dummy_equation_eval(this)
    !! evaluate dummy equation
    class(dummy_equation), intent(inout) :: this

    ! do nothing
    IGNORE(this)
  end subroutine

  subroutine selector_equation_init(this, v1, v2, status)
    !! initialize selector equation
    class(selector_equation), target, intent(out) :: this
    type(vselector),          target, intent(in)  :: v1
      !! result var selector (select v1 ...)
    type(vselector_ptr),              intent(in)  :: v2(:)
      !! dependency var selectors (... from v2(:))
    logical,                          intent(out) :: status
      !! return success/fail: success=true, fail=false

    integer              :: i, i0, i1, j, iprov, ival1, ival2tmp, itab1
    integer              :: ivsel2(v1%ntab), ival2(v1%ntab), itab2(v1%ntab), idx2(v1%g%idx_dim), perm(v1%ntab)
    logical, allocatable :: valmsk(:,:)
    type(jacobian_ptr)   :: jaco
    type(stencil_ptr)    :: st(v1%ntab)

    ! init base
    call this%equation_init("Select"//v1%name)

    ! provide v1
    iprov = this%provide(v1)

    ! stencils
    call this%dir_st%init(v1%g)
    call this%emp_st%init()

    ! loop over v1 variables
    do ival1 = 1, v1%nval
      itab2 = -1
      ! loop over source var selectors
      do i = 1, size(v2)
        ! find variable in v2(i)%p%v which corresponds to ival1
        do ival2tmp = 1, v2(i)%p%nval
          if (associated(v1%v(ival1)%p, target = v2(i)%p%v(ival2tmp)%p)) exit
        end do
        if (ival2tmp > v2(i)%p%nval) cycle

        ! check if any tables from ival, ival2tmp match
        do itab1 = 1, v1%ntab
          if (itab2(itab1) > 0) cycle
          do j = 1, v2(i)%p%ntab
            if (associated(v1%tab(itab1)%p, target = v2(i)%p%tab(j)%p)) exit
          end do
          if (j > v2(i)%p%ntab) cycle
          ivsel2(itab1) = i
          ival2( itab1) = ival2tmp
          itab2( itab1) = j
        end do
        if (all(itab2 > 0)) exit
      end do

      ! exit if some itab1 from ival1 were not found in v2
      if (any(itab2 <= 0)) then
        status = .false.
        call this%destruct() ! clean up
        return
      end if

      ! sort ivsel2 to find out how many different ivsel2 to use
      call qsort(ivsel2, perm=perm)

      i1 = 0
      do while (i1 < v1%ntab)
        i0 = i1 + 1
        ! i=i0..i1 denote same ivsel2
        do i1 = i0, v1%ntab-1
          if (ivsel2(i1+1) /= ivsel2(i0)) exit
        end do

        allocate (valmsk(v1%nval,v2(ivsel2(i0))%p%nval), source = .false.)
        st = this%emp_st%get_ptr()
        do i = i0, i1
          valmsk(ival1,ival2(perm(i))) = .true.
          st(perm(i)) = this%dir_st%get_ptr()
        end do

        ! init jacobian
        jaco%p => this%init_jaco(iprov, this%depend(v2(ivsel2(i0))%p), st, const = .true., valmsk = valmsk)

        ! set jacobian entries
        do i = i0, i1
          itab1 = perm(i)
          do j = 1, v1%tab(itab1)%p%n
            idx2 = v1%tab(itab1)%p%get_idx(j)
            ! set derivative to 1
            call jaco%p%set(itab1, j, idx2, ival1, ival2(itab1), 1.0)
          end do
        end do

        deallocate (valmsk)
      end do
    end do

    ! finish initialization
    status = .true.
    call this%init_final()
  end subroutine

  subroutine selector_equation_eval(this)
    !! evaluate selector equation
    class(selector_equation), intent(inout) :: this

    ! do nothing
    IGNORE(this)
  end subroutine

  subroutine input_equation_init_vsel(this, v)
    !! initialize input equation (f = v - appl == 0)
    class(input_equation),    intent(out) :: this
    class(vselector), target, intent(in)  :: v
      !! input parameter variable

    ! init base
    call this%equation_init("Input"//v%name)

    ! init residual data
    call this%init_f(v)

    ! init rest
    call this%input_equation_init_body()
  end subroutine

  subroutine input_equation_init_nvar_ntab(this, v, tab, name)
    !! initialize input equation
    class(input_equation), intent(out) :: this
    type(variable_ptr),    intent(in)  :: v(:)
      !! variable pointers
    type(grid_table_ptr),  intent(in)  :: tab(:)
      !! grid table pointers
    character(*),          intent(in)  :: name
      !! selector name

    ! init base
    call this%equation_init("Input"//name)

    ! init residual data
    call this%init_f(v, tab, name)

    ! init rest
    call this%input_equation_init_body()
  end subroutine

  subroutine input_equation_init_var_ntab(this, v, tab, name)
    !! initialize input equation
    class(input_equation),  intent(out) :: this
    class(variable),        intent(in)  :: v
      !! provided variable
    type(grid_table_ptr),   intent(in)  :: tab(:)
      !! grid table pointers
    character(*), optional, intent(in)  :: name
      !! name of new var selector (default: v%name)

    ! init base
    call this%equation_init("Input"//v%name)

    ! init residual data
    call this%init_f(v, tab, name=name)

    ! init rest
    call this%input_equation_init_body()
  end subroutine

  subroutine input_equation_init_nvar_tab(this, v, name, tab)
    !! initialize input equation
    class(input_equation),      intent(out) :: this
    type(variable_ptr),         intent(in)  :: v(:)
      !! variable pointers
    character(*),               intent(in)  :: name
      !! selector name
    type(grid_table), optional, intent(in)  :: tab
      !! grid table

    ! init base
    call this%equation_init("Input"//name)

    ! init residual data
    call this%init_f(v, name, tab=tab)

    ! init rest
    call this%input_equation_init_body()
  end subroutine

  subroutine input_equation_init_var_tab(this, v, tab, name)
    !! initialize input equation
    class(input_equation),      intent(out) :: this
    class(variable),            intent(in)  :: v
      !! new provided variable
    type(grid_table), optional, intent(in)  :: tab
      !! grid table
    character(*),     optional, intent(in)  :: name
      !! name of new var selector (default: var%name)

    ! init base
    call this%equation_init("Input"//v%name)

    ! init residual data
    call this%init_f(v, tab=tab, name=name)

    ! init rest
    call this%input_equation_init_body()
  end subroutine

  subroutine input_equation_init_body(this)
    !! initialize input equation
    class(input_equation), target, intent(inout) :: this

    integer                 :: i, idx(this%mvar%g%idx_dim), idep, itab, ival
    logical                 :: valmsk(this%mvar%nval,this%mvar%nval)
    type(stencil_ptr)       :: st(this%mvar%ntab)
    type(jacobian), pointer :: jaco

    ! depend on main variable
    idep = this%depend(this%mvar)

    ! dirichlet stencil
    call this%st%init(this%mvar%g)
    do itab = 1, this%mvar%ntab
      st(itab)%p => this%st
    end do

    ! value mask
    valmsk = .false.
    do ival = 1, this%mvar%nval
      valmsk(ival,ival) = .true.
    end do

    ! init jacobian
    jaco => this%init_jaco_f(idep, st, const=.true., valmsk=valmsk)
    do ival = 1, this%mvar%nval
      do itab = 1, this%mvar%ntab
        do i = 1, this%mvar%tab(itab)%p%n
          idx = this%mvar%tab(itab)%p%get_idx(i)
          call jaco%set(itab, i, idx, ival, ival, 1.0)
        end do
      end do
    end do

    ! allocate input parameters
    allocate (this%appl(this%mvar%n), source = 0.0)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine input_equation_apply(this, appl)
    !! apply input parameters
    class(input_equation), intent(inout) :: this
    real,                  intent(in)    :: appl(:)

    ! save input parameters
    call this%mvar%set(appl)
    this%appl = appl
  end subroutine

  subroutine input_equation_eval(this)
    !! evaluate input equation (f = v - appl == 0)
    class(input_equation), intent(inout) :: this

    call this%f%set(this%mvar%get() - this%appl)
  end subroutine

end module
