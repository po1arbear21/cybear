#include "../util/macro.f90.inc"

module esystem_m
  use array_m,            only: array_int
  use error_m
  use equation_m,         only: equation
  use esystem_depgraph_m, only: depgraph, STATUS_DEP
  use matrix_m,           only: block_real, matrix_real, sparse_real, sparse_cmplx, spbuild_real, matrix_convert
  use newton_m,           only: newton_opt, newton
  use res_equation_m,     only: res_equation
  use simple_equations_m, only: vector_dummy_equation_ptr, dummy_equation_ptr, selector_equation_ptr, vector_selector_equation_ptr
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
    procedure :: solve           => esystem_solve

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
    class(esystem),          intent(inout) :: this
    class(equation), target, intent(in)    :: e
      !! equation to add

    ! make sure initialization of equation is finished
    if (.not. e%finished_init) call program_error("init_final was not called for "//e%name)

    ! add equation to dependency graph
    call this%g%add_equ(e)
  end subroutine

  subroutine esystem_try_fix(this)
    !! try to generate missing equations
    class(esystem), intent(inout) :: this

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

  subroutine esystem_init_final(this)
    !! finish initialization
    class(esystem), intent(inout) :: this

    integer :: i, j, k, imvar, ires, itab1, itab2, ibl1, ibl2
    logical :: fail

    if (this%finished_init) call program_error("init_final called multiple times")
    this%finished_init = .true.

    ! try to generate selector equations for missing variable selectors
    call this%try_fix()

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

    ! analyze dependency graph
    call this%g%analyze()

    ! count blocks
    this%nbl = 0
    do i = 1, this%g%imvar%n
      imvar = this%g%imvar%d(i)
      associate (v => this%g%nodes%d(imvar)%v)
        this%nbl = this%nbl + v%ntab
      end associate
    end do

    ! set block indices
    allocate (this%res2block(this%g%imvar%n))
    allocate (this%block2res(this%nbl,2), source = 0)
    allocate (this%i0(this%nbl), source = 0)
    allocate (this%i1(this%nbl), source = 0)
    ibl1 = 0
    k    = 0
    do i = 1, this%g%imvar%n
      imvar = this%g%imvar%d(i)
      associate (v => this%g%nodes%d(imvar)%v)
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

    if (this%nbl > 0) then
      this%n = this%i1(this%nbl)
    else
      this%n = 0
    end if

    ! allocate jacobians
    call this%df%init( this%i1 - this%i0 + 1)
    call this%dft%init(this%i1 - this%i0 + 1)
    allocate (this%dfconst(this%nbl,this%nbl))

    ! set jacobians (loop over residuals and main variables)
    do i = 1, this%g%ires%n
      ires = this%g%ires%d(i)
      do j = 1, this%g%imvar%n
        imvar = this%g%imvar%d(j)
        associate (fn => this%g%nodes%d(this%g%equs%d(ires)%ires), &
          &        vn => this%g%nodes%d(imvar))
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

  subroutine esystem_eval(this, f, df)
    !! evaluate equations, get residuals and jacobians
    class(esystem),    intent(inout) :: this ! equation system
    real,              intent(out)   :: f(:)
      !! output residuals
    type(sparse_real), intent(out)   :: df
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

    ! output jacobian
    call this%get_df(df)
  end subroutine

  function esystem_get_main_var(this, i) result(mv)
    !! get main var selector
    class(esystem), intent(in) :: this
    integer,        intent(in) :: i
      !! main var index
    class(vselector), pointer  :: mv
      !! return pointer to main var

    mv => this%g%nodes%d(this%g%imvar%d(i))%v
  end function

  function esystem_get_res_equ(this, i) result(re)
    !! get residual equation
    class(esystem),   intent(in) :: this
    integer,          intent(in) :: i
      !! residual equation index
    class(res_equation), pointer :: re
      !! return pointer to residual equation

    select type (e => this%g%equs%d(this%g%ires%d(i))%e)
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
    do j = 1, this%g%imvar%n
      associate (n => this%g%nodes%d(this%g%imvar%d(j)))
        if (n%v%name == name) then
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
    imvar = this%g%imvar%d(this%block2res(ibl,1))
    associate (n => this%g%nodes%d(imvar))
      ! get table index
      itab = this%block2res(ibl,2)

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
    imvar = this%g%imvar%d(this%block2res(ibl,1))
    associate (n => this%g%nodes%d(imvar))
      ! get table index
      itab = this%block2res(ibl,2)

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
    imvar = this%g%imvar%d(this%block2res(ibl,1))
    associate (n => this%g%nodes%d(imvar))
      ! get table index
      itab = this%block2res(ibl,2)

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
