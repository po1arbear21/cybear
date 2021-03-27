#include "../util/macro.f90.inc"

module esystem_m
  use array_m,            only: array_int
  use error_m
  use equation_m,         only: equation
  use esystem_depgraph_m, only: depgraph, STATUS_DEP
  use grid_m,             only: grid_table, grid_table_ptr
  use matrix_m,           only: block_real, matrix_real, sparse_real, sparse_cmplx, spbuild_real, matrix_convert
  use newton_m,           only: newton_opt, newton
  use res_equation_m,     only: res_equation
  use simple_equations_m, only: vector_dummy_equation_ptr, dummy_equation_ptr, selector_equation_ptr, vector_selector_equation_ptr
  use variable_m,         only: variable, variable_ptr
  use vector_m,           only: vector_int
  use vselector_m,        only: vselector, vselector_ptr, vector_vselector_ptr

  implicit none

  private
  public esystem

  type esystem
    !! equation system

    character(:), allocatable :: name
      !! system name

    type(depgraph) :: g
      !! dependency graph

    type(vector_dummy_equation_ptr)    :: edum
      !! dummy equations
    type(vector_selector_equation_ptr) :: eselect
      !! selector equations

    integer                      :: nbl
      !! number of blocks
    integer                      :: n
      !! total number of values
    type(array_int), allocatable :: res2block(:)
      !! (iimvar, main var table index) -> block index;  (ntab) x (size(requs))
    integer,         allocatable :: block2res(:,:)
      !! block index -> (iimvar, main var table index); (2, nbl)
    integer,         allocatable :: i0(:)
      !! start rows for flat residuals (= cols of main vars)
    integer,         allocatable :: i1(:)
      !! end   rows for flat residuals (= cols of main vars)

    type(block_real)     :: df
      !! derivatives of f wrt x
    logical, allocatable :: dfconst(:,:)
      !! constant flags for df
    type(block_real)     :: dft
      !! derivatives of f wrt d/dt(x); constant

    logical :: finished_init
      !! indicates whether init_final was called
  contains
    procedure :: init            => esystem_init
    procedure :: destruct        => esystem_destruct
    procedure :: add_equation    => esystem_add_equation
    generic   :: provide         => esystem_provide_vsel,      &
      &                             esystem_provide_nvar_ntab, &
      &                             esystem_provide_var_ntab,  &
      &                             esystem_provide_nvar_tab,  &
      &                             esystem_provide_var_tab
    procedure :: init_final      => esystem_init_final
    procedure :: eval            => esystem_eval
    procedure :: get_main_var    => esystem_get_main_var
    procedure :: get_res_equ     => esystem_get_res_equ
    procedure :: search_main_var => esystem_search_main_var
    procedure :: solve           => esystem_solve
    generic   :: get_x           => esystem_get_x,    esystem_get_x_block
    generic   :: set_x           => esystem_set_x,    esystem_set_x_block
    generic   :: update_x        => esystem_update_x, esystem_update_x_block
    generic   :: get_df          => esystem_get_df,   esystem_get_df_cmplx
    generic   :: get_dft         => esystem_get_dft,  esystem_get_dft_cmplx
    procedure :: print           => esystem_print

    procedure, private :: esystem_provide_vsel,      &
      &                   esystem_provide_nvar_ntab, &
      &                   esystem_provide_var_ntab,  &
      &                   esystem_provide_nvar_tab,  &
      &                   esystem_provide_var_tab
    procedure, private :: esystem_set_x,    esystem_set_x_block
    procedure, private :: esystem_update_x, esystem_update_x_block
    procedure, private :: esystem_get_df,   esystem_get_df_cmplx
    procedure, private :: esystem_get_dft,  esystem_get_dft_cmplx
    procedure, private :: esystem_get_x,    esystem_get_x_block
  end type

contains

  subroutine esystem_init(this, name)
    !! initialize equation system
    class(esystem), intent(out) :: this
    character(*),   intent(in)  :: name

    this%name = name

    ! init dependency graph
    call this%g%init()

    ! init vectors
    call this%edum%init(   0, 8)
    call this%eselect%init(0, 8)

    ! init_final has not been called yet
    this%finished_init = .false.
  end subroutine

  subroutine esystem_destruct(this)
    !! release memory
    class(esystem), intent(inout) :: this

    integer :: i

    ! destruct dependency graph
    call this%g%destruct()

    ! destruct dummy and selection equations
    do i = 1, this%edum%n
      if (associated(this%edum%d(i)%p)) then
        call this%edum%d(i)%p%destruct()
        deallocate (this%edum%d(i)%p)
      end if
    end do
    call this%edum%destruct()
    do i = 1, this%eselect%n
      if (associated(this%eselect%d(i)%p)) then
        call this%eselect%d(i)%p%destruct()
        deallocate (this%eselect%d(i)%p)
      end if
    end do
    call this%eselect%destruct()

    ! destruct jacobians
    call this%df%destruct()
    call this%dft%destruct()
  end subroutine

  subroutine esystem_provide_vsel(this, v)
    !! provide vselector (create dummy equation)
    class(esystem),          intent(inout) :: this
    type(vselector), target, intent(in)    :: v

    type(dummy_equation_ptr) :: ep

    allocate (ep%p)
    call ep%p%init(v)
    call this%edum%push(ep)

    call this%add_equation(ep%p)
  end subroutine

  subroutine esystem_provide_nvar_ntab(this, v, tab, name)
    !! initialize dummy equation
    class(esystem),       intent(inout) :: this
    type(variable_ptr),   intent(in)    :: v(:)
      !! variable pointers
    type(grid_table_ptr), intent(in)    :: tab(:)
      !! grid table pointers
    character(*),         intent(in)    :: name
      !! selector name

    type(dummy_equation_ptr) :: ep

    allocate (ep%p)
    call ep%p%init(v, tab, name)
    call this%edum%push(ep)

    call this%add_equation(ep%p)
  end subroutine

  subroutine esystem_provide_var_ntab(this, v, tab, name)
    !! initialize dummy equation
    class(esystem),         intent(inout) :: this
    class(variable),        intent(in)    :: v
      !! provided variable
    type(grid_table_ptr),   intent(in)    :: tab(:)
      !! grid table pointers
    character(*), optional, intent(in)    :: name
      !! name of new var selector (default: v%name)

    type(dummy_equation_ptr) :: ep

    allocate (ep%p)
    call ep%p%init(v, tab, name=name)
    call this%edum%push(ep)

    call this%add_equation(ep%p)
  end subroutine

  subroutine esystem_provide_nvar_tab(this, v, name, tab)
    !! initialize dummy equation
    class(esystem),             intent(inout) :: this
    type(variable_ptr),         intent(in)    :: v(:)
      !! variable pointers
    character(*),               intent(in)    :: name
      !! selector name
    type(grid_table), optional, intent(in)    :: tab
      !! grid table (default: variables' whole grids via v%g%tab_all)

    type(dummy_equation_ptr) :: ep

    allocate (ep%p)
    call ep%p%init(v, name, tab=tab)
    call this%edum%push(ep)

    call this%add_equation(ep%p)
  end subroutine

  subroutine esystem_provide_var_tab(this, v, tab, name)
    !! initialize dummy equation
    class(esystem),             intent(inout) :: this
    class(variable),            intent(in)    :: v
      !! new provided variable
    type(grid_table), optional, intent(in)    :: tab
      !! grid table (default: variable's whole grid via v%g%tab_all)
    character(*),     optional, intent(in)    :: name
      !! name of new var selector (default: var%name)

    type(dummy_equation_ptr) :: ep

    allocate (ep%p)
    call ep%p%init(v, tab=tab, name=name)
    call this%edum%push(ep)

    call this%add_equation(ep%p)
  end subroutine

  subroutine esystem_add_equation(this, e)
    !! add equation to system
    class(esystem),          intent(inout) :: this
    class(equation), target, intent(in)    :: e
      !! equation to add

    ! make sure initialization of equation is finished
    if (.not. e%finished_init) call program_error("init_final was not called for "//e%name)

    ! add equation to dependency graph
    call this%g%add_equ(e)
  end subroutine

  subroutine esystem_init_final(this)
    !! finish initialization
    class(esystem), intent(inout) :: this

    if (this%finished_init) call program_error("init_final called multiple times")
    this%finished_init = .true.

    ! try to generate selector equations for missing variable selectors
    call try_fix()

    call check_provided()

    ! analyze dependency graph
    call this%g%analyze()

    call init_blocks()
    call init_jacobians()

  contains

    subroutine try_fix()
      !! try to generate missing equations

      integer                     :: i
      logical                     :: status
      type(vselector_ptr)         :: vp
      type(vector_vselector_ptr)  :: prov
      type(vector_int)            :: fix_list
      type(selector_equation_ptr) :: ep

      ! get list of provided vars, init fix list (unprovided)
      call prov%init(0, 32)
      call fix_list%init(0, 32)
      do i = 1, this%g%nodes%n
        associate (n => this%g%nodes%d(i))
          if (n%status == STATUS_DEP) then
            call fix_list%push(i)
          else
            vp%p => n%v
            call prov%push(vp)
          end if
        end associate
      end do

      ! fix nodes
      do i = 1, fix_list%n
        associate (n => this%g%nodes%d(fix_list%d(i)))
          ! new selector equation
          allocate (ep%p)
          call ep%p%init(n%v, prov%d(1:prov%n), status)

          ! check status
          if (.not. status) then
            ! can not be resolved by selector equation
            deallocate (ep%p)

            print "(A)", "Missing variable selector:"
            call n%v%print()
            call program_error("Cannot fix equation system")
          end if

          ! save equation and add to system
          call this%eselect%push(ep)
          call this%add_equation(ep%p)

          ! add var to prov
          vp%p => n%v
          call prov%push(vp)
        end associate
      end do
    end subroutine

    subroutine check_provided()
      integer :: i
      logical :: fail

      ! check if all vars are provided/main vars
      fail = .false.
      do i = 1, this%g%nodes%n
        associate (n => this%g%nodes%d(i))
          if (n%status == STATUS_DEP) then
            fail = .true.
            print "(A)", "Not provided:"
            call n%v%print()
            print *
          end if
        end associate
      end do
      if (fail) call program_error("missing one or more variable selectors")
    end subroutine

    subroutine init_blocks()
      !! initialize block index tables etc.

      integer :: iimvar, ibl, itab

      ! count blocks
      this%nbl = 0
      do iimvar = 1, this%g%imvar%n
        associate (v => this%g%nodes%d(this%g%imvar%d(iimvar))%v)
          this%nbl = this%nbl + v%ntab
        end associate
      end do

      ! set block indices
      allocate (this%res2block(this%g%imvar%n), this%block2res(2,this%nbl), this%i0(this%nbl), this%i1(this%nbl))
      ibl    = 0
      this%n = 0
      do iimvar = 1, this%g%imvar%n
        associate (v => this%g%nodes%d(this%g%imvar%d(iimvar))%v)
          allocate (this%res2block(iimvar)%d(v%ntab))

          ! mvar is split into v%ntab blocks
          do itab = 1, v%ntab
            ibl = ibl + 1

            ! set residual equation <-> block translation tables
            this%res2block(iimvar)%d(itab) = ibl
            this%block2res(1,ibl)          = iimvar
            this%block2res(2,ibl)          = itab

            ! set flat start/end indices for block
            this%i0(ibl) = this%n + 1
            this%n       = this%n + v%nvals(itab)
            this%i1(ibl) = this%n
          end do
        end associate
      end do
    end subroutine

    subroutine init_jacobians()
      !! initialize jacobians

      integer :: iimvar, jimvar, itab1, itab2, ibl1, ibl2

      ! allocate jacobians
      call this%df%init( this%i1 - this%i0 + 1)
      call this%dft%init(this%i1 - this%i0 + 1)
      allocate (this%dfconst(this%nbl,this%nbl))

      ! set jacobians (loop over residuals and main variables)
      do iimvar = 1, this%g%ires%n
        do jimvar = 1, this%g%imvar%n
          associate (fn => this%g%nodes%d(this%g%equs%d(this%g%ires%d(iimvar))%ires), &
            &        vn => this%g%nodes%d(this%g%imvar%d(jimvar))                     )
            ! loop over blocks from both main variables
            do itab1 = 1, fn%v%ntab
              do itab2 = 1, vn%v%ntab
                ! get block indices
                ibl1 = this%res2block(iimvar)%d(itab1)
                ibl2 = this%res2block(jimvar)%d(itab2)

                ! matrix size zero ?
                if ((fn%v%nvals(itab1) <= 0) .or. (vn%v%nvals(itab2) <= 0)) cycle

                ! init df matrix
                if (associated(fn%total_jaco(jimvar)%p)) then
                  ! set pointer
                  call this%df%set_ptr(ibl1, ibl2, fn%total_jaco(jimvar)%p%b(itab1,itab2)%p)
                  this%dfconst(ibl1,ibl2) = fn%total_jaco(jimvar)%p%const(itab1,itab2)
                end if

                ! init dft matrix
                if (associated(fn%total_jaco_t(jimvar)%p)) then
                  if (.not. fn%total_jaco_t(jimvar)%p%const(itab1,itab2)) then
                    call program_error("time derivative matrix is not constant")
                  end if

                  ! set pointer
                  call this%dft%set_ptr(ibl1, ibl2, fn%total_jaco_t(jimvar)%p%b(itab1,itab2)%p)
                end if
              end do
            end do
          end associate
        end do
      end do
    end subroutine

  end subroutine

  subroutine esystem_eval(this, f, df)
    !! evaluate equations, get residuals and jacobians
    class(esystem),              intent(inout) :: this ! equation system
    real,              optional, intent(out)   :: f(:)
      !! output residuals
    type(sparse_real), optional, intent(out)   :: df
      !! output jacobian

    integer :: i, j, k0, k1, ibl1, ibl2

    if (.not. this%finished_init) call program_error("init_final was not called")

    ! loop over evaluation list
    do i = 1, this%g%ieval%n
      ! get dependency graph equation
      associate (e => this%g%equs%d(this%g%ieval%d(i)))
        ! evaluate equation
        call e%e%eval()
        call e%e%set_jaco_matr(const = .false., nonconst = .true.)

        ! perform non-const jacobian chain operations for provided vars
        do j = 1, size(e%iprov)
          call this%g%nodes%d(e%iprov(j))%eval()
        end do

        ! perform non-const jacobian chain operations for residuals
        if (e%ires > 0) then
          call this%g%nodes%d(e%ires)%eval()
        end if
      end associate
    end do

    ! set residuals flat array
    if (present(f)) then
      do i = 1, this%g%ires%n
        associate (fn => this%g%nodes%d(this%g%equs%d(this%g%ires%d(i))%ires))
          ! get block indices
          ibl1 = this%res2block(i)%d(1)
          ibl2 = this%res2block(i)%d(fn%v%ntab)

          ! get flat indices
          k0 = this%i0(ibl1)
          k1 = this%i1(ibl2)

          ! set residuals
          f(k0:k1) = fn%v%get()
        end associate
      end do
    end if

    ! output jacobian
    if (present(df)) then
      call this%get_df(df)
    end if
  end subroutine

  function esystem_get_main_var(this, iimvar) result(mv)
    !! get main var selector
    class(esystem), intent(in) :: this
    integer,        intent(in) :: iimvar
      !! main var index
    class(vselector), pointer  :: mv
      !! return pointer to main var

    mv => this%g%nodes%d(this%g%imvar%d(iimvar))%v
  end function

  function esystem_get_res_equ(this, iires) result(re)
    !! get residual equation
    class(esystem),   intent(in) :: this
    integer,          intent(in) :: iires
      !! residual equation index
    class(res_equation), pointer :: re
      !! return pointer to residual equation

    select type (e => this%g%equs%d(this%g%ires%d(iires))%e)
      class is (res_equation)
        re => e
      class default
        re => null()
    end select
  end function

  function esystem_search_main_var(this, name) result(iimvar)
    !! search for main var selector by name
    class(esystem), intent(in) :: this
    character(*),   intent(in) :: name
      !! main var name
    integer                    :: iimvar
      !! return main var index

    integer :: jimvar
    logical :: found

    found = .false.
    do jimvar = 1, this%g%imvar%n
      associate (n => this%g%nodes%d(this%g%imvar%d(jimvar)))
        if (n%v%name == name) then
          if (found) call program_error("multiple main variables with name "//name//" found")
          found = .true.
          iimvar = jimvar
        end if
      end associate
    end do

    if (.not. found) then
      call program_error("main variable with name "//name//" not found")
    end if
  end function

  subroutine esystem_solve(this, opt)
    !! solves equations system by newton-raphson method.
    class(esystem),             intent(inout) :: this
    type(newton_opt), optional, intent(in)    :: opt

    real                      :: p(0)
    real, allocatable         :: x(:)
    type(newton_opt)          :: opt_
    type(sparse_real), target :: df

    if (present(opt)) then
      opt_ = opt
    else
      call opt_%init(this%n)
    end if

    allocate (x(this%n))

    call newton(fun, p, opt_, this%get_x(), x)

    ! save newton's result in esystem's variables
    call this%set_x(x)

  contains
    subroutine fun(x, p, f, dfdx, dfdp)
      real,                        intent(in)  :: x(:)
        !! arguments
      real,                        intent(in)  :: p(:)
        !! parameters
      real,                        intent(out) :: f(:)
        !! output function values
      class(matrix_real), pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
      real, optional,              intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p

      ! params not needed
      ASSERT(size(p) == 0)
      ASSERT(.not. present(dfdp))
      IGNORE(p)
      IGNORE(dfdp)

      ! save input variable
      call this%set_x(x)

      ! compute residue and jacobian
      call this%eval(f, df)
      dfdx => df
    end subroutine
  end subroutine

  function esystem_get_x(this) result(x)
    !! get main variables in flat array
    class(esystem), intent(in) :: this
    real                       :: x(this%n)
      !! return values in flat array

    integer :: ibl

    ! get values for all blocks
    do ibl = 1, this%nbl
      x(this%i0(ibl):this%i1(ibl)) = this%get_x(ibl)
    end do
  end function

  function esystem_get_x_block(this, ibl) result(x)
    !! get main variables from one block in flat array
    class(esystem), intent(in) :: this
    integer,        intent(in) :: ibl
      !! block index
    real                       :: x(this%i1(ibl)-this%i0(ibl)+1)
      !! return values in flat array

    integer :: itab, imvar

    ! get values
    imvar = this%g%imvar%d(this%block2res(1,ibl))
    associate (n => this%g%nodes%d(imvar))
      ! get table index
      itab = this%block2res(2,ibl)

      ! return values
      x = n%v%get(itab)
    end associate
  end function

  subroutine esystem_set_x(this, x)
    !! set main variables in flat array
    class(esystem), intent(inout) :: this
    real,           intent(in)    :: x(:)
      !! set values from flat array

    integer :: ibl

    ! set values for all blocks
    do ibl = 1, this%nbl
      call this%set_x(ibl, x(this%i0(ibl):this%i1(ibl)))
    end do
  end subroutine

  subroutine esystem_set_x_block(this, ibl, x)
    !! set main variables from one block in flat array
    class(esystem), intent(inout) :: this
    integer,        intent(in)    :: ibl
      !! block index
    real,           intent(in)    :: x(:)
      !! set values from flat array

    integer :: itab, imvar

    ! set values
    imvar = this%g%imvar%d(this%block2res(1,ibl))
    associate (n => this%g%nodes%d(imvar))
      ! get table index
      itab = this%block2res(2,ibl)

      ! set values
      call n%v%set(itab, x)
    end associate
  end subroutine

  subroutine esystem_update_x(this, dx)
    !! update main variables in flat array
    class(esystem), intent(inout) :: this
    real,           intent(in)    :: dx(:)
      !! delta values in flat array

    integer :: ibl

    ! update values for all blocks
    do ibl = 1, this%nbl
      call this%update_x(ibl, dx(this%i0(ibl):this%i1(ibl)))
    end do
  end subroutine

  subroutine esystem_update_x_block(this, ibl, dx)
    !! update main variables from one block in flat array
    class(esystem), intent(inout) :: this
    integer,        intent(in)    :: ibl
      !! block index
    real,           intent(in)    :: dx(:)
      !! delta values in flat array

    integer :: itab, imvar

    ! update values
    imvar = this%g%imvar%d(this%block2res(1,ibl))
    associate (n => this%g%nodes%d(imvar))
      ! get table index
      itab = this%block2res(2,ibl)

      ! update values
      call n%v%update(itab, dx)
    end associate
  end subroutine

  subroutine esystem_get_df(this, df)
    !! get df as sparse matrix
    class(esystem),    intent(in)  :: this
    type(sparse_real), intent(out) :: df
      !! output sparse matrix

    ! local variables
    type(spbuild_real) :: sb

    call df%init(this%n)
    call sb%init(df)
    call matrix_convert(this%df, sb)
    call sb%save()
  end subroutine

  subroutine esystem_get_df_cmplx(this, df)
    !! get df as complex sparse matrix
    class(esystem),     intent(in)  :: this
    type(sparse_cmplx), intent(out) :: df
      !! output sparse matrix

    ! local variables
    type(sparse_real)  :: df_real
    type(spbuild_real) :: sb

    call df_real%init(this%n)
    call sb%init(df_real)
    call matrix_convert(this%df, sb)
    call sb%save()
    call matrix_convert(df_real, df)
  end subroutine

  subroutine esystem_get_dft(this, dft)
    !! get df as sparse matrix
    class(esystem),    intent(in)  :: this
    type(sparse_real), intent(out) :: dft
      !! output sparse matrix

    ! local variables
    type(spbuild_real) :: sb

    call dft%init(this%n)
    call sb%init(dft)
    call matrix_convert(this%dft, sb)
    call sb%save()
  end subroutine

  subroutine esystem_get_dft_cmplx(this, dft)
    !! get df as complex sparse matrix
    class(esystem),     intent(in)  :: this
    type(sparse_cmplx), intent(out) :: dft
      !! output sparse matrix

    ! local variables
    type(sparse_real)  :: dft_real
    type(spbuild_real) :: sb

    call dft_real%init(this%n)
    call sb%init(dft_real)
    call matrix_convert(this%dft, sb)
    call sb%save()
    call matrix_convert(dft_real, dft)
  end subroutine

  subroutine esystem_print(this)
    !! print main variables with bounds
    class(esystem), intent(in) :: this

    ! local variables
    integer :: i, j, ibl

    ibl = 0
    do i = 1, this%g%imvar%n
      associate (v => this%g%nodes%d(this%g%imvar%d(i))%v)
        call v%print()
        do j = 1, v%ntab
          ibl = ibl + 1
          print "(I0,A,I0)", this%i0(ibl), " ", this%i1(ibl)
        end do
        print *
      end associate
    end do
  end subroutine
end module
