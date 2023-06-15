m4_include(../util/macro.f90.inc)

module equation_m

  use array_m,         only: array1_real
  use color_m
  use error_m
  use grid_table_m,    only: grid_table, grid_table_ptr
  use ieee_arithmetic, only: ieee_is_nan
  use jacobian_m,      only: jacobian, jacobian_ptr
  use stencil_m,       only: stencil_ptr
  use variable_m,      only: variable, variable_ptr
  use vector_m,        only: vector_int
  use vselector_m,     only: vselector, vselector_ptr, vector_vselector_ptr

  implicit none

  private
  public equation, equation_ptr, vector_equation_ptr
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
    type(jacobian_ptr), allocatable :: jaco_p(:,:)
      !! preconditioner derivatives of vprov wrt vdep (vprov%n x vdep%n)

    logical :: finished_init
      !! indicates whether init_final was called
  contains
    procedure :: equation_init
    procedure :: get_ptr       => equation_get_ptr
    procedure :: destruct      => equation_destruct
    procedure :: reset         => equation_reset
    generic   :: provide       => equation_provide_vsel,      &
      &                           equation_provide_nvar_ntab, &
      &                           equation_provide_var_ntab,  &
      &                           equation_provide_nvar_tab,  &
      &                           equation_provide_var_tab
    generic   :: depend        => equation_depend_vsel,       &
      &                           equation_depend_nvar_ntab,  &
      &                           equation_depend_var_ntab,   &
      &                           equation_depend_nvar_tab,   &
      &                           equation_depend_var_tab
    procedure :: realloc_jaco  => equation_realloc_jaco
    procedure :: init_jaco     => equation_init_jaco
    procedure :: set_jaco_matr => equation_set_jaco_matr
    procedure :: init_final    => equation_init_final
    procedure :: has_precon    => equation_has_precon
    procedure :: test          => equation_test

    procedure(equation_eval), deferred :: eval

    procedure, private :: equation_provide_vsel
    procedure, private :: equation_provide_nvar_ntab
    procedure, private :: equation_provide_var_ntab
    procedure, private :: equation_provide_nvar_tab
    procedure, private :: equation_provide_var_tab
    procedure, private :: equation_depend_vsel
    procedure, private :: equation_depend_nvar_ntab
    procedure, private :: equation_depend_var_ntab
    procedure, private :: equation_depend_nvar_tab
    procedure, private :: equation_depend_var_tab
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

  m4_define({T},{equation_ptr})
  m4_include(../util/vector_def.f90.inc)

contains

  m4_define({T},{equation_ptr})
  m4_include(../util/vector_imp.f90.inc)

  subroutine equation_init(this, name, precon)
    !! initialize equation base
    class(equation),   intent(out) :: this
    character(*),      intent(in)  :: name
      !! name of equation
    logical, optional, intent(in)  :: precon
      !! preconditioner flag (will preconditioner be used?). (default: .false.)

    integer, parameter :: CAP = 16
    logical            :: precon_

    ! optional argument
    precon_ = .false.
    if (present(precon)) precon_ = precon

    ! save name
    this%name = name

    ! init vprov and vdep vectors
    call this%vprov%init(    0, c = CAP)
    call this%vdep%init(     0, c = CAP)
    call this%vprov_alc%init(0, c = CAP)
    call this%vdep_alc%init( 0, c = CAP)

    ! allocate jaco + preconditioner
    allocate (this%jaco(CAP,CAP))
    if (precon_) allocate (this%jaco_p(CAP,CAP))

    ! init_final has not been called yet
    this%finished_init = .false.
  end subroutine

  function equation_get_ptr(this) result(ptr)
    !! return pointer to this equation
    class(equation), target, intent(in) :: this
    type(equation_ptr)                  :: ptr

    ptr%p => this
  end function

  subroutine equation_destruct(this)
    !! destruct equation
    class(equation), intent(inout) :: this

    integer :: i, j

    if (allocated(this%name)) deallocate (this%name)

    ! destruct jacobians
    if (allocated(this%jaco)) then
      do j = 1, size(this%jaco, dim = 2)
        do i = 1, size(this%jaco, dim = 1)
          if (associated(this%jaco(i,j)%p)) then
            call this%jaco(i,j)%p%destruct()
            deallocate (this%jaco(i,j)%p)
          end if
        end do
      end do
      deallocate (this%jaco)
    end if
    if (allocated(this%jaco_p)) then
      do j = 1, size(this%jaco_p, dim = 2)
        do i = 1, size(this%jaco_p, dim = 1)
          if (associated(this%jaco_p(i,j)%p)) then
            call this%jaco_p(i,j)%p%destruct()
            deallocate (this%jaco_p(i,j)%p)
          end if
        end do
      end do
      deallocate (this%jaco_p)
    end if

    ! destruct vselectors
    if (allocated(this%vprov_alc%d)) then
      do i = 1, this%vprov_alc%n
        deallocate (this%vprov%d(this%vprov_alc%d(i))%p)
      end do
    end if
    if (allocated(this%vdep_alc%d)) then
      do i = 1, this%vdep_alc%n
        deallocate (this%vdep%d(this%vdep_alc%d(i))%p)
      end do
    end if
    call this%vprov%destruct()
    call this%vdep%destruct()
    call this%vprov_alc%destruct()
    call this%vdep_alc%destruct()
  end subroutine

  subroutine equation_reset(this, const, nonconst)
    !! reset provided var selectors and jacobians
    class(equation),   intent(inout) :: this
    logical, optional, intent(in)    :: const
      !! reset constant parts (default: true)
    logical, optional, intent(in)    :: nonconst
      !! reset non-constant parts (default: true)

    integer :: i, j, itab
    logical :: const_, nonconst_, vprov_const

    ! optional arguments
    const_ = .true.
    if (present(const)) const_ = const
    nonconst_ = .true.
    if (present(nonconst)) nonconst_ = nonconst

    ! reset provided var selectors
    do i = 1, this%vprov%n
      do itab = 1, this%vprov%d(i)%p%ntab
        vprov_const = .true.
        do j = 1, this%vdep%n
          if (associated(this%jaco(i,j)%p)) vprov_const = all(this%jaco(i,j)%p%matr%const(itab,:))
          if (.not. vprov_const) exit
        end do
        if ((const_ .and. vprov_const) .or. (nonconst_ .and. .not. vprov_const)) call this%vprov%d(i)%p%reset(itab = itab)
      end do
    end do

    ! reset jacobians
    do j = 1, this%vdep%n
      do i = 1, this%vprov%n
        if (associated(this%jaco(i,j)%p)) call this%jaco(i,j)%p%reset(const = const, nonconst = nonconst)
        if (allocated(this%jaco_p)) then
          if (associated(this%jaco_p(i,j)%p)) call this%jaco_p(i,j)%p%reset(const = const, nonconst = nonconst)
        end if
      end do
    end do
  end subroutine

  subroutine equation_realloc_jaco(this, cprov, cdep)
    !! reallocate this%jaco + this%jaco_p, if initial capacity was not big enough
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

    ! reallocate this%jaco_p
    if (allocated(this%jaco_p)) then
      allocate (jaco_tmp(cprov, cdep))
      jaco_tmp(1:this%vprov%n,1:this%vdep%n) = this%jaco_p(1:this%vprov%n,1:this%vdep%n)
      call move_alloc(jaco_tmp, this%jaco_p)
    end if
  end subroutine

  function equation_init_jaco(this, iprov, idep, st, const, zero, valmsk, valmsk_tab, precon) result(jaco)
    !! allocate and initialize jacobian
    class(equation),             intent(inout) :: this
    integer,                     intent(in)    :: iprov
      !! provided var selector index
    integer,                     intent(in)    :: idep
      !! dependency var selector index
    type(stencil_ptr), optional, intent(in)    :: st(:)
      !! stencils (v1%ntab); default = dirichtlet_stencil (works only if provided and dependend variables are defined on same grid)
    logical,           optional, intent(in)    :: const
      !! const flag; default = false; the same for all blocks
    logical,           optional, intent(in)    :: zero(:,:)
      !! zero flags (v1%ntab x v2%ntab); default: set automatically by checking stencils
    logical,           optional, intent(in)    :: valmsk(:,:)
      !! value mask (v1%nval x v2%nval); default = true; the same for all blocks
    logical,           optional, intent(in)    :: valmsk_tab(:,:,:,:)
      !! value mask per block (v1%nval x v2%nval x v1%ntab x v2%ntab)
    logical,           optional, intent(in)    :: precon
      !! create preconditioner jacobian; default = false
    type(jacobian),    pointer                 :: jaco
      !! return pointer to newly created jacobian

    logical :: precon_

    associate (vprov => this%vprov%d(iprov)%p, vdep => this%vdep%d(idep)%p)
      ! preconditioner flag
      precon_ = .false.
      if (present(precon)) precon_ = precon
      m4_assert((.not. precon_) .or. allocated(this%jaco_p))

      ! allocate and init jacobian
      allocate (jaco)
      call jaco%init(vprov, vdep, st = st, const = const, zero = zero, valmsk = valmsk, valmsk_tab = valmsk_tab)

      ! save pointer to jacobian
      if (precon_) then
        m4_assert(.not. associated(this%jaco_p(iprov,idep)%p))
        this%jaco_p(iprov,idep)%p => jaco
      else
        m4_assert(.not. associated(this%jaco(iprov,idep)%p))
        this%jaco(iprov,idep)%p => jaco
      end if
    end associate
  end function

  subroutine equation_init_final(this)
    !! set constant parts of jacobian matrices
    class(equation), intent(inout) :: this

    m4_assert(.not. this%finished_init)
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
        if (allocated(this%jaco_p)) then
          if (associated(this%jaco_p(i,j)%p)) call this%jaco_p(i,j)%p%set_matr(const = const, nonconst = nonconst)
        end if
      end do
    end do
  end subroutine

  function equation_has_precon(this, iprov, idep) result(precon)
    !! determines if block (iprov, idep) has a preconditioner
    class(equation), intent(in) :: this
    integer,         intent(in) :: iprov
      !! provided var selector index
    integer,         intent(in) :: idep
      !! dependency var selector index
    logical                     :: precon

    precon = .false.
    if (allocated(this%jaco_p)) precon = associated(this%jaco_p(iprov,idep)%p)
  end function

  subroutine equation_test(this, ntest, rx, ax, rtol, atol, no_sign_change)
    !! test jacobians with finite differences
    class(equation),   intent(inout) :: this
    integer,           intent(in)    :: ntest
      !! number of tests to perform
    real,    optional, intent(in)    :: rx
      !! relative x perturbation (default: 1e-4)
    real,    optional, intent(in)    :: ax
      !! absolute x perturbation (default: 1e-8)
    real,    optional, intent(in)    :: rtol
      !! relative tolerance (default: 1e-3)
    real,    optional, intent(in)    :: atol
      !! absolute tolerance (default: 1e-6)
    logical, optional, intent(in)    :: no_sign_change
      !! forbid sign change when calculating finite differences (default: false)

    integer                        :: i, j, k, kk, l
    logical                        :: nan, nan1, nan2, nsgn
    real                           :: rx_, ax_, rtol_, atol_, r
    real                           :: dx1, dydx, dydx1, dydx2
    real,              allocatable :: x0(:), xm(:), xp(:), dx(:)
    type(array1_real), allocatable :: y0(:), ym(:), yp(:), dy(:)

    ! optional arguments
    rx_ = 1e-4
    if (present(rx)) rx_ = rx
    ax_ = 1e-8
    if (present(ax)) ax_ = ax
    rtol_ = 1e-3
    if (present(rtol)) rtol_ = rtol
    atol_ = 1e-6
    if (present(atol)) atol_ = atol
    nsgn = .false.
    if (present(no_sign_change)) nsgn = .true.

    ! evaluate equation and set jacobians in matrix form
    call this%eval()
    call this%set_jaco_matr()

    ! get provided variables
    allocate (y0(this%vprov%n), ym(this%vprov%n), yp(this%vprov%n), dy(this%vprov%n))
    do i = 1, this%vprov%n
      allocate (y0(i)%d(this%vprov%d(i)%p%n), ym(i)%d(this%vprov%d(i)%p%n), yp(i)%d(this%vprov%d(i)%p%n), dy(i)%d(this%vprov%d(i)%p%n))
      y0(i)%d = this%vprov%d(i)%p%get()
    end do

    do j = 1, this%vdep%n
      associate (dep => this%vdep%d(j)%p)
        ! allocate memory for dependencies
        allocate (x0(dep%n), dx(dep%n), xm(dep%n), xp(dep%n))

        ! get original x values
        x0 = dep%get()
        xm = x0
        xp = x0
        dx = 0

        ! test all derivatives
        do kk = 1, ntest
          if (mod(kk-1, 100) == 0) print *, kk
          call random_number(r)
          k = ceiling(size(x0) * r)

          ! get delta x
          dx1 = max(abs(x0(k)) * rx_, ax_)
          if (nsgn) dx1 = min(dx1, abs(x0(k) * (1.0 - 1e-16)))

          ! set xp, xm, dx
          xm(k) = x0(k) - dx1
          xp(k) = x0(k) + dx1
          dx(k) = dx1

          ! get ym
          call dep%set(xm)
          call this%eval()
          do i = 1, this%vprov%n
            ym(i)%d = this%vprov%d(i)%p%get()
          end do

          ! get yp
          call dep%set(xp)
          call this%eval()
          do i = 1, this%vprov%n
            yp(i)%d = this%vprov%d(i)%p%get()
          end do

          ! check derivatives by finite differences
          do i = 1, this%vprov%n
            ! does d(f_i)/d(x_j) exist? zero by default
            if (.not. associated(this%jaco(i,j)%p)) cycle

            ! get dy = jaco * dx
            call this%jaco(i,j)%p%matr%mul_vec(dx, dy(i)%d)

            associate (prov => this%vprov%d(i)%p)
              do l = 1, prov%n
                ! check if both derivatives are zero
                if ((yp(i)%d(l) == y0(i)%d(l)) .and. (dy(i)%d(l) == 0)) cycle

                ! entry from matrix
                dydx = dy(i)%d(l) / dx1

                ! finite differences (forward and centered)
                dydx1 = (yp(i)%d(l) - y0(i)%d(l)) / dx1
                dydx2 = (yp(i)%d(l) - ym(i)%d(l)) / (2 * dx1)

                ! check for nan
                nan  = ieee_is_nan(dydx)
                nan1 = ieee_is_nan(dydx1)
                nan2 = ieee_is_nan(dydx2)
                if (nan .or. nan1 .or. nan2) then
                  print *
                  print *
                  print "(A)", "provided:"
                  call prov%print()
                  print *
                  print "(A)", "dependency:"
                  call dep%print()

                  if (nan ) call program_error("Matrix entry is NaN")
                  if (nan1) call program_error("Forward finite difference is NaN")
                  if (nan2) call program_error("Central finite difference is NaN")
                end if

                if (abs(dydx - dydx2) > max(max(2 * abs(dydx2 - dydx1), atol_), rtol_ * abs(dydx2))) then
                  print *
                  print *
                  print "(A)", COL_MAGENTA//"Possible Error detected:"//COL_DEFAULT
                  print "(A)", "provided:"
                  call prov%print()
                  print *
                  print "(A)", "dependency:"
                  call dep%print()
                  print *
                  print "(A, I5,      A, I5     )", "k(dep) = ", k,     ";                    l(prov) = ", l
                  print "(A, ES24.16, A, ES24.16)", "xm     = ", xm(k), "; ym      = ", ym(i)%d(l)
                  print "(A, ES24.16, A, ES24.16)", "x0     = ", x0(k), "; y0      = ", y0(i)%d(l)
                  print "(A, ES24.16, A, ES24.16)", "xp     = ", xp(k), "; yp      = ", yp(i)%d(l)
                  print "(A, ES24.16, A         )", "dydx   = ", dydx,  " (from matrix       )"
                  print "(A, ES24.16, A         )", "dydx1  = ", dydx1, " (fin.diff. forward )"
                  print "(A, ES24.16, A         )", "dydx2  = ", dydx2, " (fin.diff. centered)"
                end if
              end do
            end associate
          end do

          ! reset xp, xm, dx
          xp(k) = x0(k)
          xm(k) = x0(k)
          dx(k) = 0
        end do

        ! restore x
        call dep%set(x0)

        ! free memory
        deallocate (x0, xm, xp, dx)
      end associate
    end do
  end subroutine

  function equation_provide_vsel(this, vsel) result(iprov)
    !! provide var selector
    class(equation),         intent(inout) :: this
    type(vselector), target, intent(in)    :: vsel
      !! new provided var selector
    integer                                :: iprov
      !! return provided index

    ! reallocate jaco if necessary
    if (this%vprov%n >= size(this%vprov%d)) call this%realloc_jaco((this%vprov%n + 1) * 2, size(this%vdep%d))

    ! add var selector to provided variables
    call this%vprov%push(vsel%get_ptr())

    ! return index
    iprov = this%vprov%n
  end function

  function equation_provide_nvar_ntab(this, v, tab, name) result(iprov)
    !! provide variables for multiple grid tables, creates var selector internally
    class(equation),      intent(inout) :: this
    type(variable_ptr),   intent(in)    :: v(:)
      !! variable pointers
    type(grid_table_ptr), intent(in)    :: tab(:)
      !! grid table pointers
    character(*),         intent(in)    :: name
      !! selector name
    integer                             :: iprov

    type(vselector), pointer :: vsel

    ! create vselector from variable and keep track of memory
    allocate (vsel)
    call vsel%init(v, tab, name)
    iprov = this%provide(vsel)
    call this%vprov_alc%push(iprov)
  end function

  function equation_provide_var_ntab(this, v, tab, name) result(iprov)
    !! provide variable for multiple grid tables, creates var selector internally
    class(equation),        intent(inout) :: this
    class(variable),        intent(in)    :: v
      !! new provided variable
    type(grid_table_ptr),   intent(in)    :: tab(:)
      !! grid table pointers
    character(*), optional, intent(in)    :: name
      !! name of new var selector (default: v%name)
    integer                               :: iprov
      !! return provided index

    type(vselector), pointer :: vsel

    ! create vselector from variable and keep track of memory
    allocate (vsel)
    call vsel%init(v, tab, name=name)
    iprov = this%provide(vsel)
    call this%vprov_alc%push(iprov)
  end function

  function equation_provide_nvar_tab(this, v, name, tab) result(iprov)
    !! provide variables for multiple grid tables, creates var selector internally
    class(equation),            intent(out) :: this
    type(variable_ptr),         intent(in)  :: v(:)
      !! variable pointers
    character(*),               intent(in)  :: name
      !! selector name
    type(grid_table), optional, intent(in)  :: tab
      !! grid table (default: variable's whole grid via v%g%tab_all)
    integer                                 :: iprov
      !! return provided index

    type(vselector), pointer :: vsel

    ! create vselector from variable and keep track of memory
    allocate (vsel)
    call vsel%init(v, name, tab=tab)
    iprov = this%provide(vsel)
    call this%vprov_alc%push(iprov)
  end function

  function equation_provide_var_tab(this, v, tab, name) result(iprov)
    !! provide variable for single grid table, creates var selector internally
    class(equation),            intent(inout) :: this
    class(variable),            intent(in)    :: v
      !! new provided variable
    type(grid_table), optional, intent(in)    :: tab
      !! grid table (default: variable's whole grid via v%g%tab_all)
    character(*),     optional, intent(in)    :: name
      !! name of new var selector (default: v%name)
    integer                                   :: iprov
      !! return provided index

    type(vselector), pointer :: vsel

    ! create vselector from variable and keep track of memory
    allocate (vsel)
    call vsel%init(v, tab=tab, name=name)
    iprov = this%provide(vsel)
    call this%vprov_alc%push(iprov)
  end function

  function equation_depend_vsel(this, vsel) result(idep)
    !! depend on var selector
    class(equation),          intent(inout) :: this
    class(vselector), target, intent(in)    :: vsel
      !! new dependency var selector
    integer                                 :: idep
      !! return dependency index

    ! reallocate g if necessary
    if (this%vdep%n >= size(this%vdep%d)) call this%realloc_jaco(size(this%vprov%d), (this%vdep%n + 1) * 2)

    ! add var selector to dependent variables
    call this%vdep%push(vsel%get_ptr())

    ! return dependency index
    idep = this%vdep%n
  end function

  function equation_depend_nvar_ntab(this, v, tab, name) result(idep)
    !! depend on variables for multiple grid tables, creates var selector internally
    class(equation),      intent(inout) :: this
    type(variable_ptr),   intent(in)    :: v(:)
      !! variable pointers
    type(grid_table_ptr), intent(in)    :: tab(:)
      !! grid table pointers
    character(*),         intent(in)    :: name
      !! selector name
    integer                             :: idep

    type(vselector), pointer :: vsel

    ! create vselector from variable and keep track of memory
    allocate (vsel)
    call vsel%init(v, tab, name)
    idep = this%depend(vsel)
    call this%vdep_alc%push(idep)
  end function

  function equation_depend_var_ntab(this, v, tab, name) result(idep)
    !! depend on variable for multiple grid tables, creates var selector internally
    class(equation),        intent(inout) :: this
    class(variable),        intent(in)    :: v
      !! new dependent variable
    type(grid_table_ptr),   intent(in)    :: tab(:)
      !! grid table pointers
    character(*), optional, intent(in)    :: name
      !! name of new var selector (default: var%name)
    integer                               :: idep
      !! return dependency index

    type(vselector), pointer :: vsel

    ! create vselector from variable and keep track of memory
    allocate (vsel)
    call vsel%init(v, tab, name=name)
    idep = this%depend(vsel)
    call this%vdep_alc%push(idep)
  end function

  function equation_depend_nvar_tab(this, v, name, tab) result(idep)
    !! depend on variables for multiple grid tables, creates var selector internally
    class(equation),            intent(out) :: this
    type(variable_ptr),         intent(in)  :: v(:)
      !! variable pointers
    character(*),               intent(in)  :: name
      !! selector name
    type(grid_table), optional, intent(in)  :: tab
      !! grid table (default: variable's whole grid via v%g%tab_all)
    integer                                 :: idep
      !! return dependency index

    type(vselector), pointer :: vsel

    ! create vselector from variable and keep track of memory
    allocate (vsel)
    call vsel%init(v, name, tab=tab)
    idep = this%depend(vsel)
    call this%vdep_alc%push(idep)
  end function

  function equation_depend_var_tab(this, v, tab, name) result(idep)
    !! depend variable for single grid table, creates var selector internally
    class(equation),            intent(inout) :: this
    class(variable),            intent(in)    :: v
      !! new dependent variable
    type(grid_table), optional, intent(in)    :: tab
      !! grid table (default: variable's whole grid via v%g%tab_all)
    character(*),     optional, intent(in)    :: name
      !! name of new var selector (default: var%name)
    integer                                   :: idep
      !! return dependency index

    type(vselector), pointer :: vsel

    ! create vselector from variable and keep track of memory
    allocate (vsel)
    call vsel%init(v, tab=tab, name=name)
    idep = this%depend(vsel)
    call this%vdep_alc%push(idep)
  end function

end module
