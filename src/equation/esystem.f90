module esystem_m
  use array_m
  use esystem_dag_m
  use matrix_m
  use simple_equations_m
  implicit none

  type esystem
    !! equation system

    character(:), allocatable :: name
      !! system name

    type(dag) :: d
      !! directed acyclic graph

    type(vector_dummy_equation_ptr)    :: edum
      !! dummy equations
    type(vector_selector_equation_ptr) :: eselect
      !! selector equations

    integer                      :: nbl
      !! number of blocks
    integer                      :: n
      !! total number of values
    type(array_int), allocatable :: res2block(:)
      !! (res equation index, main var table index) -> block index;  (ntab) x (size(requs))
    integer,         allocatable :: block2res(:,:)
      !! block index -> (res equation index, main var table index); (nbl x 2)
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
    procedure :: provide         => esystem_provide
    procedure :: add_equation    => esystem_add_equation
    procedure :: try_fix         => esystem_try_fix
    procedure :: init_final      => esystem_init_final
    procedure :: eval            => esystem_eval
    procedure :: get_main_var    => esystem_get_main_var
    procedure :: get_res_equ     => esystem_get_res_equ
    procedure :: search_main_var => esystem_search_main_var

    procedure :: esystem_get_x
    procedure :: esystem_get_x_block
    generic   :: get_x => esystem_get_x, esystem_get_x_block

    procedure :: esystem_set_x
    procedure :: esystem_set_x_block
    generic   :: set_x => esystem_set_x, esystem_set_x_block

    procedure :: esystem_update_x
    procedure :: esystem_update_x_block
    generic   :: update_x => esystem_update_x, esystem_update_x_block

    procedure :: esystem_get_df
    procedure :: esystem_get_df_cmplx
    generic   :: get_df => esystem_get_df, esystem_get_df_cmplx

    procedure :: esystem_get_dft
    procedure :: esystem_get_dft_cmplx
    generic   :: get_dft => esystem_get_dft, esystem_get_dft_cmplx

    procedure :: print => esystem_print
  end type

contains

  subroutine esystem_init(this, name)
    !! initialize equation system
    class(esystem), intent(out) :: this
    character(*),   intent(in)  :: name

    this%name = name

    ! init dag
    call this%d%init()

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

    ! destruct dag
    call this%d%destruct()

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

  subroutine esystem_provide(this, v)
    !! provide vselector (create dummy equation)
    class(esystem),          intent(inout) :: this
    type(vselector), target, intent(in)    :: v

    type(dummy_equation_ptr) :: ep

    allocate (ep%p)
    call ep%p%init(v)
    call this%edum%push(ep)

    call this%add_equation(ep%p)
  end subroutine

  subroutine esystem_add_equation(this, e)
    !! add equation to system
    class(esystem),       intent(inout) :: this
    class(equation), target, intent(in)    :: e
      !! equation to add

    ! make sure initialization of equation is finished
    if (.not. e%finished_init) call program_error("init_final was not called for "//e%name)

    ! add equation to dag
    call this%d%add_equ(e)
  end subroutine

  subroutine esystem_try_fix(this)
    !! try to generate missing equations
    class(esystem), intent(inout) :: this

    integer                     :: i
    logical                     :: status
    type(vselector_ptr)         :: vp
    type(vector_vselector_ptr)  :: prov
    type(vector_dag_node_ptr)   :: fix_list
    type(selector_equation_ptr) :: ep

    ! get list of provided vars, init fix list (unprovided)
    call prov%init(0, 32)
    call fix_list%init(0, 32)
    do i = 1, this%d%nodes%n
      associate (np => this%d%nodes%d(i))
        if (np%p%status == NDSTATUS_DEP) then
          call fix_list%push(np)
        else
          vp%p => np%p%v
          call prov%push(vp)
        end if
      end associate
    end do

    ! fix nodes
    do i = 1, fix_list%n
      associate (n => fix_list%d(i)%p)
        ! new selector equation
        allocate (ep%p)
        call ep%p%init(n%v, prov%d(1:prov%n), status)

        ! check status
        if (.not. status) then
          ! can not be resolved by selector equation
          deallocate (ep%p)

          print "(A)", "Missing variable selector:"
          call n%v%print()
          call program_error("Can not fix equation system")
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

  subroutine esystem_init_final(this)
    !! finish initialization
    class(esystem), intent(inout) :: this

    integer :: i, j, k, mvari, resi, itab1, itab2, ibl1, ibl2
    logical :: fail

    if (this%finished_init) call program_error("init_final called multiple times")
    this%finished_init = .true.

    ! try to generate selector equations for missing variable selectors
    call this%try_fix()

    ! check if all vars are provided/main vars
    fail = .false.
    do i = 1, this%d%nodes%n
      associate (n => this%d%nodes%d(i)%p)
        if (n%status == NDSTATUS_DEP) then
          fail = .true.
          print "(A)", "Not provided:"
          call n%v%print()
          print *
        end if
      end associate
    end do
    if (fail) call program_error("missing one or more variable selectors")

    ! analyze dag
    call this%d%analyze()

    ! count blocks
    this%nbl = 0
    do i = 1, this%d%mvari%n
      mvari = this%d%mvari%d(i)
      associate (v => this%d%nodes%d(mvari)%p%v)
        this%nbl = this%nbl + v%ntab
      end associate
    end do

    ! set block indices
    allocate (this%res2block(this%d%mvari%n))
    allocate (this%block2res(this%nbl,2), source = 0)
    allocate (this%i0(this%nbl), source = 0)
    allocate (this%i1(this%nbl), source = 0)
    ibl1 = 0
    k    = 0
    do i = 1, this%d%mvari%n
      mvari = this%d%mvari%d(i)
      associate (v => this%d%nodes%d(mvari)%p%v)
        allocate (this%res2block(i)%d(v%ntab), source = 0)

        ! mvar is split into v%ntab blocks
        do j = 1, v%ntab
          ibl1 = ibl1 + 1

          ! set residual equation <-> block translation tables
          this%res2block(i)%d(j) = ibl1
          this%block2res(ibl1,1) = i
          this%block2res(ibl1,2) = j

          ! set flat start/end indices for block
          this%i0(ibl1) = k + 1
          this%i1(ibl1) = k + v%nvals(j)
          k             = this%i1(ibl1)
        end do
      end associate
    end do
    this%n = this%i1(this%nbl)

    ! allocate jacobians
    call this%df%init( this%i1 - this%i0 + 1)
    call this%dft%init(this%i1 - this%i0 + 1)
    allocate (this%dfconst(this%nbl,this%nbl))

    ! set jacobians (loop over residuals and main variables)
    do i = 1, this%d%resi%n
      resi = this%d%resi%d(i)
      do j = 1, this%d%mvari%n
        mvari = this%d%mvari%d(j)
        associate (fn => this%d%equs%d(resi)%p%res, &
          &        vn => this%d%nodes%d(mvari)%p)
          ! loop over blocks from both main variables
          do itab1 = 1, fn%v%ntab
            do itab2 = 1, vn%v%ntab
              ! get block indices
              ibl1 = this%res2block(i)%d(itab1)
              ibl2 = this%res2block(j)%d(itab2)

              ! matrix size zero ?
              if ((fn%v%nvals(itab1) <= 0) .or. (vn%v%nvals(itab2) <= 0)) cycle

              ! init df matrix
              if (associated(fn%total_jaco(j)%p)) then
                ! set pointer
                call this%df%set_ptr(ibl1, ibl2, fn%total_jaco(j)%p%b(itab1,itab2)%p)
                this%dfconst(ibl1,ibl2) = fn%total_jaco(j)%p%const(itab1,itab2)
              end if

              ! init dft matrix
              if (associated(fn%total_jaco_t(j)%p)) then
                if (.not. fn%total_jaco_t(j)%p%const(itab1,itab2)) then
                  call program_error("time derivative matrix is not constant")
                end if

                ! set pointer
                call this%dft%set_ptr(ibl1, ibl2, fn%total_jaco_t(j)%p%b(itab1,itab2)%p)
              end if
            end do
          end do
        end associate
      end do
    end do
  end subroutine

  subroutine esystem_eval(this, f)
    !! evaluate equations, get residuals and jacobians
    class(esystem), intent(inout) :: this ! equation system
    real,           intent(out)   :: f(:)
      !! output residuals

    integer :: i, j, k0, k1, ibl1, ibl2

    if (.not. this%finished_init) call program_error("init_final was not called")

    ! loop over evaluation list
    do i = 1, this%d%ev%n
      ! get DAG equation
      associate (dag_e => this%d%equs%d(this%d%ev%d(i))%p)
        ! evaluate equation
        call dag_e%e%eval()
        call dag_e%e%set_jaco_matr(const = .false., nonconst = .true.)

        ! perform non-const graph operations for provided vars
        do j = 1, size(dag_e%prov)
          call dag_e%prov(j)%p%eval()
        end do

        ! perform non-const graph operations for residual
        if (associated(dag_e%res)) then
          call dag_e%res%eval()
        end if
      end associate
    end do

    ! set residuals flat array
    do i = 1, this%d%resi%n
      associate (fn => this%d%equs%d(this%d%resi%d(i))%p%res)
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
  end subroutine

  function esystem_get_main_var(this, i) result(mv)
    !! get main var selector
    class(esystem), intent(in) :: this
    integer,        intent(in) :: i
      !! main var index
    class(vselector), pointer  :: mv
      !! return pointer to main var

    mv => this%d%nodes%d(this%d%mvari%d(i))%p%v
  end function

  function esystem_get_res_equ(this, i) result(re)
    !! get residual equation
    class(esystem),   intent(in) :: this
    integer,          intent(in) :: i
      !! residual equation index
    class(res_equation), pointer :: re
      !! return pointer to residual equation

    select type (e => this%d%equs%d(this%d%resi%d(i))%p%e)
      class is (res_equation)
        re => e
      class default
        re => null()
    end select
  end function

  function esystem_search_main_var(this, name) result(i)
    !! search for main var selector by name
    class(esystem), intent(in) :: this
    character(*),   intent(in) :: name
      !! main var name
    integer                    :: i
      !! return main var index

    integer :: j
    logical :: found

    found = .false.
    do j = 1, this%d%mvari%n
      associate (nd => this%d%nodes%d(this%d%mvari%d(j))%p)
        if (nd%v%name == name) then
          if (found) call program_error("multiple main variables with name "//name//" found")
          found = .true.
          i = j
        end if
      end associate
    end do

    if (.not. found) then
      call program_error("main variable with name "//name//" not found")
    end if
  end function

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

    integer :: itab, mvari

    ! get values
    mvari = this%d%mvari%d(this%block2res(ibl,1))
    associate (nd => this%d%nodes%d(mvari)%p)
      ! get table index
      itab = this%block2res(ibl,2)

      ! return values
      x = nd%v%get(itab)
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

    integer :: itab, mvari

    ! set values
    mvari = this%d%mvari%d(this%block2res(ibl,1))
    associate (nd => this%d%nodes%d(mvari)%p)
      ! get table index
      itab = this%block2res(ibl,2)

      ! set values
      call nd%v%set(itab, x)
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

    integer :: itab, mvari

    ! update values
    mvari = this%d%mvari%d(this%block2res(ibl,1))
    associate (nd => this%d%nodes%d(mvari)%p)
      ! get table index
      itab = this%block2res(ibl,2)

      ! update values
      call nd%v%update(itab, dx)
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
    call this%df%to_sparse(sb)
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
    call this%df%to_sparse(sb)
    call sb%save()
    call df_real%to_cmplx(df)
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
    call this%dft%to_sparse(sb)
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
    call this%dft%to_sparse(sb)
    call sb%save()
    call dft_real%to_cmplx(dft)
  end subroutine

  subroutine esystem_print(this)
    !! print main variables with bounds
    class(esystem), intent(in) :: this

    ! local variables
    integer :: i, j, ibl

    ibl = 0
    do i = 1, this%d%mvari%n
      associate (v => this%d%nodes%d(this%d%mvari%d(i))%p%v)
        call v%print()
        do j = 1, v%ntab
          ibl = ibl + 1
          print "(2I6)", this%i0(ibl), this%i1(ibl)
        end do
        print *
      end associate
    end do
  end subroutine

end module