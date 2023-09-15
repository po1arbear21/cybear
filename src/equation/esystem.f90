m4_include(../util/macro.f90.inc)

module esystem_m

  use array_m,            only: array1_int
  use color_m,            only: COL_DEFAULT, COL_MAGENTA
  use error_m,            only: assert_failed, program_error
  use equation_m,         only: equation, equation_ptr, vector_equation_ptr
  use esystem_depgraph_m, only: depgraph, depgraph_equ, STATUS_DEP
  use gmres_m,            only: gmres_options
  use grid_table_m,       only: grid_table, grid_table_ptr
  use hashmap_m,          only: hashmap_int
  use json_m,             only: json_object
  use matrix_m,           only: block_real, dense_real, dense_cmplx, matrix_real, matrix_cmplx, sparse_real, sparse_cmplx, matrix_convert
  use newton_m,           only: newton_opt, newton
  use output_file_m,      only: output_file
  use res_equation_m,     only: res_equation
  use simple_equations_m, only: dummy_equation, selector_equation, input_equation
  use variable_m,         only: variable, variable_ptr, vector_variable_ptr
  use vector_m,           only: vector_int
  use vselector_m,        only: vselector, vselector_ptr, vector_vselector_ptr

  implicit none

  private
  public esystem

  type esystem
    !! equation system

    character(:), allocatable :: name
      !! system name

    type(vector_variable_ptr) :: vars
      !! list of all variables
    type(hashmap_int)         :: hvars
      !! hashmap for variable => vars index

    type(depgraph) :: g
      !! dependency graph

    type(vector_equation_ptr) :: ealloc
      !! automatically allocated equations

    integer                       :: nbl
      !! number of blocks
    integer                       :: n
      !! total number of values
    type(array1_int), allocatable :: res2block(:)
      !! (iimvar, main var table index) -> block index;  (ntab) x (size(requs))
    integer,          allocatable :: block2res(:,:)
      !! block index -> (iimvar, main var table index); (2, nbl)
    integer,          allocatable :: i0(:)
      !! start rows for flat residuals (= cols of main vars)
    integer,          allocatable :: i1(:)
      !! end   rows for flat residuals (= cols of main vars)

    type(vector_int)     :: input_equs
      !! input equation indices
    integer, allocatable :: input_i0(:)
      !! start rows for input variables
    integer, allocatable :: input_i1(:)
      !! end   rows for input variables
    integer              :: ninput
      !! total number of input values

    logical                       :: dense
      !! use dense df, dft
    type(block_real)              :: df
      !! derivatives of f wrt x
    logical,          allocatable :: dfconst(:,:)
      !! constant flags for df
    type(block_real), allocatable :: dfp
      !! preconditioner derivatives of f wrt x
      !! allocatable: indicates if preconditioner is being computed
    type(block_real)              :: dft
      !! derivatives of f wrt d/dt(x); constant

    logical :: finished_init
      !! indicates whether init_final was called
    logical :: parallel_eval
      !! Use parallel evaluation of equations in eval()
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
    generic   :: set_input       => esystem_set_input_all, esystem_set_input_single
    procedure :: eval            => esystem_eval
    procedure :: get_main_var    => esystem_get_main_var
    procedure :: get_res_equ     => esystem_get_res_equ
    procedure :: search_main_var => esystem_search_main_var
    procedure :: solve           => esystem_solve
    generic   :: get_x           => esystem_get_x,    esystem_get_x_block
    generic   :: set_x           => esystem_set_x,    esystem_set_x_block
    generic   :: update_x        => esystem_update_x, esystem_update_x_block
    generic   :: get_df          => esystem_get_df,   esystem_get_df_cmplx
    generic   :: get_dfp         => esystem_get_dfp,  esystem_get_dfp_cmplx
    generic   :: get_dft         => esystem_get_dft,  esystem_get_dft_cmplx
    procedure :: print           => esystem_print
    procedure :: output_info     => esystem_output_info
    procedure :: output_data     => esystem_output_data

    procedure, private :: esystem_provide_vsel,      &
      &                   esystem_provide_nvar_ntab, &
      &                   esystem_provide_var_ntab,  &
      &                   esystem_provide_nvar_tab,  &
      &                   esystem_provide_var_tab
    procedure, private :: esystem_set_input_all, esystem_set_input_single
    procedure, private :: esystem_set_x,    esystem_set_x_block
    procedure, private :: esystem_update_x, esystem_update_x_block
    procedure, private :: esystem_get_df,   esystem_get_df_cmplx
    procedure, private :: esystem_get_dfp,  esystem_get_dfp_cmplx
    procedure, private :: esystem_get_dft,  esystem_get_dft_cmplx
    procedure, private :: esystem_get_x,    esystem_get_x_block
  end type

contains

  subroutine esystem_init(this, name, precon, dense, parallel_eval)
    !! initialize equation system
    class(esystem),    intent(out) :: this
    character(*),      intent(in)  :: name
    logical, optional, intent(in)  :: precon
      !! should a preconditioner matrix be created? (default: false)
    logical, optional, intent(in)  :: dense
      !! use dense jacobian matrix instead of sparse
    logical, optional, intent(in)  :: parallel_eval
      !! use parallel evaluation of equations

    logical :: precon_

    this%name = name

    ! preconditoner
    precon_ = .false.
    if (present(precon)) precon_ = precon
    if (precon_) allocate(this%dfp)

    ! dense jacobian?
    this%dense = .false.
    if (present(dense)) this%dense = dense

    ! parallel eval?
    this%parallel_eval = .false.
    if (present(parallel_eval)) this%parallel_eval = parallel_eval

    ! init variable list
    call this%vars%init(0, c = 8)
    call this%hvars%init()

    ! init dependency graph
    call this%g%init(precon_)

    ! init allocated equation vector
    call this%ealloc%init(0, c = 8)

    ! input variable indices
    call this%input_equs%init(0, c = 8)

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
    do i = 1, this%ealloc%n
      if (associated(this%ealloc%d(i)%p)) then
        call this%ealloc%d(i)%p%destruct()
        deallocate (this%ealloc%d(i)%p)
      end if
    end do
    call this%ealloc%destruct()

    ! destruct jacobians
    call this%df%destruct()
    call this%dft%destruct()
    if (allocated(this%dfp)) call this%dfp%destruct()
  end subroutine

  subroutine esystem_provide_vsel(this, v, input)
    !! provide vselector (create dummy equation)
    class(esystem),          intent(inout) :: this
    type(vselector), target, intent(in)    :: v
    logical, optional,       intent(in)    :: input
      !! input variable (default: false)

    logical            :: input_
    type(equation_ptr) :: ptr

    input_ = .false.
    if (present(input)) input_ = input

    if (input_) then
      allocate (input_equation :: ptr%p)
      call this%ealloc%push(ptr)
      select type (i => ptr%p)
      type is (input_equation)
        call i%init(v)
        call this%add_equation(i)
        call this%input_equs%push(this%g%equs%n)
      end select
    else
      allocate (dummy_equation :: ptr%p)
      call this%ealloc%push(ptr)
      select type (d => ptr%p)
      type is (dummy_equation)
        call d%init(v)
        call this%add_equation(d)
      end select
    end if
  end subroutine

  subroutine esystem_provide_nvar_ntab(this, v, tab, name, input)
    !! initialize dummy equation
    class(esystem),       intent(inout) :: this
    type(variable_ptr),   intent(in)    :: v(:)
      !! variable pointers
    type(grid_table_ptr), intent(in)    :: tab(:)
      !! grid table pointers
    character(*),         intent(in)    :: name
      !! selector name
    logical, optional,    intent(in)    :: input
      !! input variable (default: false)

    logical            :: input_
    type(equation_ptr) :: ptr

    input_ = .false.
    if (present(input)) input_ = input

    if (input_) then
      allocate (input_equation :: ptr%p)
      call this%ealloc%push(ptr)
      select type (i => ptr%p)
      type is (input_equation)
        call i%init(v, tab, name)
        call this%add_equation(i)
        call this%input_equs%push(this%g%equs%n)
      end select
    else
      allocate (dummy_equation :: ptr%p)
      call this%ealloc%push(ptr)
      select type (d => ptr%p)
      type is (dummy_equation)
        call d%init(v, tab, name)
        call this%add_equation(d)
      end select
    end if
  end subroutine

  subroutine esystem_provide_var_ntab(this, v, tab, name, input)
    !! initialize dummy equation
    class(esystem),         intent(inout) :: this
    class(variable),        intent(in)    :: v
      !! provided variable
    type(grid_table_ptr),   intent(in)    :: tab(:)
      !! grid table pointers
    character(*), optional, intent(in)    :: name
      !! name of new var selector (default: v%name)
    logical,      optional, intent(in)    :: input
      !! input variable (default: false)

    logical            :: input_
    type(equation_ptr) :: ptr

    input_ = .false.
    if (present(input)) input_ = input

    if (input_) then
      allocate (input_equation :: ptr%p)
      call this%ealloc%push(ptr)
      select type (i => ptr%p)
      type is (input_equation)
        call i%init(v, tab, name=name)
        call this%add_equation(i)
        call this%input_equs%push(this%g%equs%n)
      end select
    else
      allocate (dummy_equation :: ptr%p)
      call this%ealloc%push(ptr)
      select type (d => ptr%p)
      type is (dummy_equation)
        call d%init(v, tab, name=name)
        call this%add_equation(d)
      end select
    end if
  end subroutine

  subroutine esystem_provide_nvar_tab(this, v, name, tab, input)
    !! initialize dummy equation
    class(esystem),             intent(inout) :: this
    type(variable_ptr),         intent(in)    :: v(:)
      !! variable pointers
    character(*),               intent(in)    :: name
      !! selector name
    type(grid_table), optional, intent(in)    :: tab
      !! grid table (default: whole grid via v%g%tab_all)
    logical,          optional, intent(in)    :: input
      !! input variable (default: false)

    logical            :: input_
    type(equation_ptr) :: ptr

    input_ = .false.
    if (present(input)) input_ = input

    if (input_) then
      allocate (input_equation :: ptr%p)
      call this%ealloc%push(ptr)
      select type (i => ptr%p)
      type is (input_equation)
        call i%init(v, name, tab=tab)
        call this%add_equation(i)
        call this%input_equs%push(this%g%equs%n)
      end select
    else
      allocate (dummy_equation :: ptr%p)
      call this%ealloc%push(ptr)
      select type (d => ptr%p)
      type is (dummy_equation)
        call d%init(v, name, tab=tab)
        call this%add_equation(d)
      end select
    end if
  end subroutine

  subroutine esystem_provide_var_tab(this, v, tab, name, input)
    !! initialize dummy equation
    class(esystem),             intent(inout) :: this
    class(variable),            intent(in)    :: v
      !! new provided variable
    type(grid_table), optional, intent(in)    :: tab
      !! grid table (default: whole grid via v%g%tab_all)
    character(*),     optional, intent(in)    :: name
      !! name of new var selector (default: var%name)
    logical,          optional, intent(in)    :: input
      !! input variable (default: false)

    logical            :: input_
    type(equation_ptr) :: ptr

    input_ = .false.
    if (present(input)) input_ = input

    if (input_) then
      allocate (input_equation :: ptr%p)
      call this%ealloc%push(ptr)
      select type (i => ptr%p)
      type is (input_equation)
        call i%init(v, tab=tab, name=name)
        call this%add_equation(i)
        call this%input_equs%push(this%g%equs%n)
      end select
    else
      allocate (dummy_equation :: ptr%p)
      call this%ealloc%push(ptr)
      select type (d => ptr%p)
      type is (dummy_equation)
        call d%init(v, tab=tab, name=name)
        call this%add_equation(d)
      end select
    end if
  end subroutine

  subroutine esystem_add_equation(this, e)
    !! add equation to system
    class(esystem),          intent(inout) :: this
    class(equation), target, intent(in)    :: e
      !! equation to add

    integer :: i

    ! make sure initialization of equation is finished
    if (.not. e%finished_init) call program_error("init_final was not called for "//e%name)

    ! add equation to dependency graph
    call this%g%add_equ(e)

    ! collect variables
    do i = 1, e%vprov%n
      call add_vselector(e%vprov%d(i)%p)
    end do
    do i = 1, e%vdep%n
      call add_vselector(e%vdep%d(i)%p)
    end do
    select type (e)
    class is (res_equation)
      call add_vselector(e%mvar)
    end select

  contains

    subroutine add_vselector(v)
      !! add variables from vselector to list
      type(vselector), target, intent(in) :: v

      integer :: i

      do i = 1, size(v%v)
        call add_variable(v%v(i)%p)
      end do
    end subroutine

    subroutine add_variable(v)
      !! add variable to list
      class(variable), target, intent(in) :: v

      integer :: ivar, hkey(m4_ptrsize)
      logical :: status

      hkey = v%hashkey()

      call this%hvars%get(hkey, ivar, status = status)
      if (.not. status) then
        call this%vars%push(v%get_ptr())
        call this%hvars%set(hkey, this%vars%n)
      end if
    end subroutine

  end subroutine

  subroutine esystem_init_final(this)
    !! finish initialization
    class(esystem), intent(inout) :: this

    m4_assert(.not. this%finished_init)
    this%finished_init = .true.

    ! try to generate selector equations for missing variable selectors
    call try_fix()

    call check_provided()

    ! analyze dependency graph
    call this%g%analyze()

    call init_blocks()
    call init_jacobians()

    call init_input()

  contains

    subroutine try_fix()
      !! try to generate missing equations

      integer, parameter :: CAP = 32
      integer                          :: i
      logical                          :: status
      type(vselector_ptr)              :: vp
      type(vector_vselector_ptr)       :: prov
      type(vector_int)                 :: fix_list
      type(equation_ptr)               :: ptr

      ! get list of provided vars, init fix list (unprovided)
      call prov%init(    0, c = CAP)
      call fix_list%init(0, c = CAP)
      do i = 1, this%g%nodes%n
        associate (n => this%g%nodes%d(i)%p)
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
        associate (n => this%g%nodes%d(fix_list%d(i))%p)
          ! new selector equation
          allocate (selector_equation :: ptr%p)
          select type (e => ptr%p)
          type is (selector_equation)
            call e%init(n%v, prov%d(1:prov%n), status)

            ! check status
            if (.not. status) then
              ! can not be resolved by selector equation
              deallocate (ptr%p)

              print "(A)", "Missing variable selector:"
              call n%v%print()
              call program_error("Cannot fix equation system")
            end if

            ! save equation and add to system
            call this%ealloc%push(ptr)
            call this%add_equation(e)
          end select

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
        associate (n => this%g%nodes%d(i)%p)
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
        associate (v => this%g%nodes%d(this%g%imvar%d(iimvar))%p%v)
          this%nbl = this%nbl + v%ntab
        end associate
      end do

      ! set block indices
      allocate (this%res2block(this%g%imvar%n), this%block2res(2,this%nbl), this%i0(this%nbl), this%i1(this%nbl))
      ibl    = 0
      this%n = 0
      do iimvar = 1, this%g%imvar%n
        associate (v => this%g%nodes%d(this%g%imvar%d(iimvar))%p%v)
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
      if (allocated(this%dfp)) call this%dfp%init(this%i1 - this%i0 + 1)
      allocate (this%dfconst(this%nbl,this%nbl))

      ! set jacobians (loop over residuals and main variables)
      do iimvar = 1, this%g%ires%n
        do jimvar = 1, this%g%imvar%n
          associate (fn => this%g%nodes%d(this%g%equs%d(this%g%ires%d(iimvar))%ires)%p, &
            &        vn => this%g%nodes%d(this%g%imvar%d(jimvar))%p                     )
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

                ! init dfp matrix
                if (associated(fn%total_jaco_p(jimvar)%p)) then
                  ! set pointer
                  call this%dfp%set_ptr(ibl1, ibl2, fn%total_jaco_p(jimvar)%p%b(itab1,itab2)%p)
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

    subroutine init_input()
      !! initialize input indices

      integer :: i, iimvar, ibl0, ibl1

      allocate (this%input_i0(this%input_equs%n), this%input_i1(this%input_equs%n))

      this%ninput = 0
      do i = 1, this%input_equs%n
        ! get mvar index
        iimvar = this%g%nodes%d(this%g%equs%d(this%input_equs%d(i))%imain)%p%iimvar

        ! get first and last block
        ibl0 = this%res2block(iimvar)%d(1)
        ibl1 = this%res2block(iimvar)%d(size(this%res2block(iimvar)%d))

        ! set start and end indices (covering all blocks)
        this%input_i0(i) = this%i0(ibl0)
        this%input_i1(i) = this%i1(ibl1)

        ! count total number of input variables
        this%ninput = this%ninput + this%input_i1(i) - this%input_i0(i) + 1
      end do
    end subroutine

  end subroutine

  subroutine esystem_set_input_all(this, x, only_appl)
    !! set input values for all input equation
    class(esystem),    intent(inout) :: this
    real,              intent(in)    :: x(:)
      !! new values for all input variables
    logical, optional, intent(in)    :: only_appl
      !! only set the applied value, leave actual value unchanged (default: false)

    integer :: i, i0, i1

    m4_assert(size(x) == this%ninput)

    i1 = 0
    do i = 1, this%input_equs%n
      i0 = i1 + 1
      i1 = i0 + this%input_i1(i) - this%input_i0(i)
      call this%esystem_set_input_single(i, x(i0:i1), only_appl)
    end do
  end subroutine

  subroutine esystem_set_input_single(this, i, x, only_appl)
    !! set input values for one input equation
    class(esystem),    intent(inout) :: this
    integer,           intent(in)    :: i
      !! input variable index
    real,              intent(in)    :: x(:)
      !! new values for i-th input variable
    logical, optional, intent(in)    :: only_appl
      !! only set the applied value, leave actual value unchanged (default: false)

    m4_assert(size(x) == this%input_i1(i)-this%input_i0(i)+1)

    select type (e => this%g%equs%d(this%input_equs%d(i))%e)
      class is (input_equation)
        call e%apply(x, only_appl)
    end select
  end subroutine

  subroutine esystem_eval(this, f, df, dfp)
    !! evaluate equations, get residuals, jacobians, and preconditioner
    class(esystem),     target,   intent(inout) :: this ! equation system
    real,               optional, intent(out)   :: f(:)
      !! output residuals
    class(matrix_real), optional, intent(out)   :: df
      !! output jacobian
    type(sparse_real),  optional, intent(out)   :: dfp
      !! output preconditioner jacobian

    integer :: i, j, k, k0, k1, ibl1, ibl2
    type(depgraph_equ), pointer :: e

    m4_assert(this%finished_init)

    ! Parallelized evaluation of equations in dependency blocks
    if (this%parallel_eval) then
      do i = 1, this%g%neblks
        !$omp parallel do default(none) schedule(dynamic) private(j,e,k) shared(i,this)
        do j = this%g%ieblks(i), this%g%ieblks(i + 1) - 1
          ! get dependency graph equation
          e => this%g%equs%d(this%g%ieval_final(j))

          ! reset non-constant parts of provided var selectors and jacobians
          call e%e%reset(const = .false., nonconst = .true.)

          ! evaluate equation
          call e%e%eval()
          call e%e%set_jaco_matr(const = .false., nonconst = .true.)

          ! perform non-const jacobian chain operations for provided vars
          do k = 1, size(e%iprov)
            call this%g%nodes%d(e%iprov(k))%p%eval()
          end do

          ! perform non-const jacobian chain operations for residuals
          if (e%ires > 0) call this%g%nodes%d(e%ires)%p%eval()
        end do
        !$omp end parallel do
      end do

    ! loop over original evaluation list
    else
      do i = 1, this%g%ieval%n
        ! get dependency graph equation
        e => this%g%equs%d(this%g%ieval%d(i))

        ! reset non-constant parts of provided var selectors and jacobians
        call e%e%reset(const = .false., nonconst = .true.)

        ! evaluate equation
        call e%e%eval()
        call e%e%set_jaco_matr(const = .false., nonconst = .true.)

        ! perform non-const jacobian chain operations for provided vars
        do j = 1, size(e%iprov)
          call this%g%nodes%d(e%iprov(j))%p%eval()
        end do

        ! perform non-const jacobian chain operations for residuals
        if (e%ires > 0) call this%g%nodes%d(e%ires)%p%eval()
      end do
    end if

    ! set residuals flat array
    if (present(f)) then
      m4_assert(size(f) == this%n)
      do i = 1, this%g%ires%n
        associate (fn => this%g%nodes%d(this%g%equs%d(this%g%ires%d(i))%ires)%p)
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
    if (present(df)) call this%get_df(df)

    ! output preconditioner jacobian
    if (present(dfp)) call this%get_dfp(dfp)
  end subroutine

  function esystem_get_main_var(this, iimvar) result(mv)
    !! get main var selector
    class(esystem), intent(in) :: this
    integer,        intent(in) :: iimvar
      !! main var index
    class(vselector), pointer  :: mv
      !! return pointer to main var

    mv => this%g%nodes%d(this%g%imvar%d(iimvar))%p%v
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
      associate (n => this%g%nodes%d(this%g%imvar%d(jimvar))%p)
        if (n%v%name == name) then
          if (found) call program_error("multiple main variables with name "//name//" found")
          found = .true.
          iimvar = jimvar
        end if
      end associate
    end do

    if (.not. found) call program_error("main variable with name "//name//" not found")
  end function

  subroutine esystem_solve(this, nopt, gopt)
    !! solves steady-state by newton-raphson method (deprecated, use analysis/steady_state instead)
    class(esystem),                intent(inout) :: this
    type(newton_opt),    optional, intent(in)    :: nopt
      !! options for the newton solver
    type(gmres_options), optional, intent(in)    :: gopt
      !! options for the iterative solver in each newton iteration (only used when nopt%it_solver == true)

    real                      :: p(0)
    real, allocatable         :: x(:)
    type(newton_opt)          :: nopt_
    type(sparse_real), target :: df, dfp

    print "(A)", COL_MAGENTA//"esystem_solve is deprecated, use analysis/steady_state instead"//COL_DEFAULT

    m4_assert(.not. this%dense)

    ! parse optional newton options
    if (present(nopt)) then
      m4_assert(size(nopt%atol) == this%n)
      nopt_ = nopt
    else
      call nopt_%init(this%n)
    end if

    allocate (x(this%n))

    call newton(fun, p, nopt_, this%get_x(), x, gmres_opt = gopt)

    ! save newton's result in esystem's variables
    call this%set_x(x)

    ! release factorizations
    if (nopt_%it_solver) then
      call dfp%destruct()
    else
      call df%destruct()
    end if

  contains

    subroutine fun(x, p, f, dfdx, dfdx_prec, dfdp)
      real,                                  intent(in)  :: x(:)
        !! arguments
      real,                                  intent(in)  :: p(:)
        !! parameters
      real,               optional,          intent(out) :: f(:)
        !! output function values
      class(matrix_real), optional, pointer, intent(out) :: dfdx
        !! output pointer to jacobian of f wrt x
      class(matrix_real), optional, pointer, intent(out) :: dfdx_prec
        !! optional output pointer to preconditioner jacobian of f wrt x
      real,               optional,          intent(out) :: dfdp(:,:)
        !! optional output jacobian of f wrt p

      ! params not needed
      m4_assert(size(p) == 0)
      m4_ignore(p)
      m4_assert(present(dfdx))
      m4_assert(.not. present(dfdp))
      m4_ignore(dfdp)

      ! save input variable
      call this%set_x(x)

      ! compute residue, jacobian, and possible preconditioner
      if (nopt_%it_solver) then
        m4_assert(present(dfdx_prec))
        call this%eval(f = f, df = df, dfp = dfp)
        dfdx_prec => dfp
      else
        m4_assert(.not. present(dfdx_prec))
        call this%eval(f = f, df = df)
      end if
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
    associate (n => this%g%nodes%d(imvar)%p)
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

    m4_assert(size(x) == this%n)

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

    m4_assert(size(x) == this%i1(ibl)-this%i0(ibl)+1)

    ! set values
    imvar = this%g%imvar%d(this%block2res(1,ibl))
    associate (n => this%g%nodes%d(imvar)%p)
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

    m4_assert(size(dx) == this%n)

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

    m4_assert(size(dx) == this%i1(ibl)-this%i0(ibl)+1)

    ! update values
    imvar = this%g%imvar%d(this%block2res(1,ibl))
    associate (n => this%g%nodes%d(imvar)%p)
      ! get table index
      itab = this%block2res(2,ibl)

      ! update values
      call n%v%update(itab, dx)
    end associate
  end subroutine

  subroutine esystem_get_df(this, df)
    !! get jacobian
    class(esystem),     intent(in)  :: this
    class(matrix_real), intent(out) :: df
      !! output sparse or dense matrix

    select type (df)
    class is (sparse_real)
      call matrix_convert(this%df, df)          ! sparse_real <- block_real
    class is (dense_real)
      call matrix_convert(this%df, df)          ! dense_real <- block_real
    end select
  end subroutine

  subroutine esystem_get_df_cmplx(this, df)
    !! get jacobian as complex matrix
    class(esystem),      intent(in)  :: this
    class(matrix_cmplx), intent(out) :: df
      !! output sparse or dense matrix

    type(dense_real)  :: d
    type(sparse_real) :: s

    select type (df)
    class is (sparse_cmplx)
      call matrix_convert(this%df, s) ! sparse_real  <- block_real
      call matrix_convert(s, df)      ! sparse_cmplx <- sparse_real
    class is (dense_cmplx)
      call matrix_convert(this%df, d) ! dense_real  <- block_real
      call matrix_convert(d, df)      ! dense_cmplx <- dense_real
    end select
  end subroutine

  subroutine esystem_get_dfp(this, dfp)
    !! get dfp as sparse matrix
    class(esystem),    intent(in)  :: this
    type(sparse_real), intent(out) :: dfp
      !! output sparse matrix

    m4_assert(allocated(this%dfp))

    call matrix_convert(this%dfp, dfp)        ! sparse_real <- block_real
  end subroutine

  subroutine esystem_get_dfp_cmplx(this, dfp)
    !! get dfp as complex sparse matrix
    class(esystem),     intent(in)  :: this
    type(sparse_cmplx), intent(out) :: dfp
      !! output sparse matrix

    type(sparse_real)  :: dfp_real

    m4_assert(allocated(this%dfp))

    call matrix_convert(this%dfp, dfp_real)   ! sparse_real  <- block_real
    call matrix_convert(dfp_real, dfp     )   ! sparse_cmplx <- sparse_real
  end subroutine

  subroutine esystem_get_dft(this, dft)
    !! get time derivative jacobian
    class(esystem),     intent(in)  :: this
    class(matrix_real), intent(out) :: dft
      !! output sparse or dense matrix

    select type (dft)
    class is (sparse_real)
      call matrix_convert(this%dft, dft)          ! sparse_real <- block_real
    class is (dense_real)
      call matrix_convert(this%dft, dft)          ! dense_real <- block_real
    end select
  end subroutine

  subroutine esystem_get_dft_cmplx(this, dft)
    !! get time derivative jacobian as complex matrix
    class(esystem),      intent(in)  :: this
    class(matrix_cmplx), intent(out) :: dft
      !! output sparse or dense matrix

    type(dense_real)  :: d
    type(sparse_real) :: s

    select type (dft)
    class is (sparse_cmplx)
      call matrix_convert(this%dft, s) ! sparse_real  <- block_real
      call matrix_convert(s, dft)      ! sparse_cmplx <- sparse_real
    class is (dense_cmplx)
      call matrix_convert(this%dft, d) ! dense_real  <- block_real
      call matrix_convert(d, dft)      ! dense_cmplx <- dense_real
    end select
  end subroutine

  subroutine esystem_print(this)
    !! print main variables with bounds
    class(esystem), intent(in) :: this

    ! local variables
    integer :: i, j, ibl

    ibl = 0
    do i = 1, this%g%imvar%n
      associate (v => this%g%nodes%d(this%g%imvar%d(i))%p%v)
        call v%print()
        do j = 1, v%ntab
          ibl = ibl + 1
          print "(I0,A,I0)", this%i0(ibl), " ", this%i1(ibl)
        end do
        print *
      end associate
    end do
  end subroutine

  subroutine esystem_output_info(this, of)
    !! output info of all variables
    class(esystem),    intent(in)    :: this
    type(output_file), intent(inout) :: of
      !! output file handle

    integer :: i

    do i = 1, this%vars%n
      call this%vars%d(i)%p%output_info(of)
    end do
  end subroutine

  subroutine esystem_output_data(this, of, obj)
    !! output data of all variables
    class(esystem),             intent(in)    :: this
    type(output_file),          intent(inout) :: of
      !! output file handle
    type(json_object), pointer, intent(inout) :: obj
      !! parent object in output file

    integer :: i

    do i = 1, this%vars%n
      call this%vars%d(i)%p%output_data(of, obj)
    end do
  end subroutine

end module
