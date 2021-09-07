#include "../util/macro.f90.inc"

module simple_equations_m

  use array_m,        only: array_int, array2_log
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
    call this%equation_init("provide_"//v%name)

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
    call this%equation_init("provide_"//name)

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
      call this%equation_init("provide_"//name)
    else
      call this%equation_init("provide_"//v%name)
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
    call this%equation_init("provide_"//name)

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
      call this%equation_init("provide_"//name)
    else
      call this%equation_init("provide_"//v%name)
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

    integer                 :: i, j, ival1, ival2, itab1, itab2, idx(v1%g%idx_dim), iprov, idep
    logical                 :: v2_used(size(v2))
    type(array_int)         :: get_v2(v1%ntab,v1%nval), get_ival2(v1%ntab,v1%nval)
    type(array2_log)        :: valmsk(size(v2))
    type(jacobian), pointer :: jaco
    type(stencil_ptr)       :: st(v1%ntab)

    ! init base
    call this%equation_init("select_"//v1%name)

    ! provide v1
    iprov = this%provide(v1)

    ! stencils
    call this%dir_st%init(v1%g)
    call this%emp_st%init()

    ! loop over v1 values
    v2_used = .false.
    do ival1 = 1, v1%nval
      do itab1 = 1, v1%ntab
        allocate (get_v2(   itab1,ival1)%d(v1%tab(itab1)%p%n), source = 0)
        allocate (get_ival2(itab1,ival1)%d(v1%tab(itab1)%p%n), source = 0)
      end do

      ! loop over source var selectors
      do j = 1, size(v2)
        ! check if grids match
        if (.not. associated(v2(j)%p%g, target = v1%g)) cycle

        ! check if idx_type, idx_dir match
        if ((v1%idx_type /= v2(j)%p%idx_type) .or. (v1%idx_dir /= v2(j)%p%idx_dir)) cycle

        ! find value in v2(j)%p%v which corresponds to ival1
        do ival2 = 1, v2(j)%p%nval
          if (associated(v1%v(ival1)%p, target = v2(j)%p%v(ival2)%p)) exit
        end do
        if (ival2 > v2(j)%p%nval) cycle

        ! allocate value mask
        if (.not. allocated(valmsk(j)%d)) allocate (valmsk(j)%d(v1%nval,v2(j)%p%nval), source = .false.)

        ! find points
        do itab1 = 1, v1%ntab
          do i = 1, v1%tab(itab1)%p%n
            if (get_v2(itab1,ival1)%d(i) > 0) cycle

            idx   = v1%tab(itab1)%p%get_idx(i)
            itab2 = v2(j)%p%itab%get(idx)

            if (itab2 > 0) then
              get_v2(   itab1,ival1)%d(i) = j
              get_ival2(itab1,ival1)%d(i) = ival2
              v2_used(j) = .true.
              valmsk(j)%d(ival1,ival2) = .true.
            end if
          end do
        end do
      end do

      ! make sure this value was found for all points
      do itab1 = 1, v1%ntab
        if (any(get_v2(itab1,ival1)%d <= 0)) then
          ! clean up and return fail
          call this%destruct()
          status = .false.
          return
        end if
      end do
    end do

    ! init dependencies
    do j = 1, size(v2)
      if (.not. v2_used(j)) cycle
      idep = this%depend(v2(j)%p)

      ! get stencils
      st = this%emp_st%get_ptr()
      do itab1 = 1, v1%ntab
        do ival1 = 1, v1%nval
          if (any(get_v2(itab1,ival1)%d == j)) then
            st(itab1) = this%dir_st%get_ptr()
            exit
          end if
        end do
      end do

      ! init jacobian
      jaco => this%init_jaco(iprov, idep, st, const = .true., valmsk = valmsk(j)%d)

      ! set jacobian entries
      do itab1 = 1, v1%ntab
        do ival1 = 1, v1%nval
          if (all(get_v2(itab1,ival1)%d /= j)) cycle
          do i = 1, v1%tab(itab1)%p%n
            idx = v1%tab(itab1)%p%get_idx(i)
            call jaco%set(itab1, i, idx, ival1, get_ival2(itab1,ival1)%d(i), 1.0)
          end do
        end do
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
    call this%equation_init("input_"//v%name)

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
    call this%equation_init("input_"//name)

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
    call this%equation_init("input_"//v%name)

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
    call this%equation_init("input_"//name)

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
    call this%equation_init("input_"//v%name)

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
