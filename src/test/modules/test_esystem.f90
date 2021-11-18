module test_esystem_m

  use equation_m,     only: equation
  use error_m,        only: program_error
  use esystem_m,      only: esystem
  use jacobian_m,     only: jacobian
  use grid_m,         only: grid, IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL, grid_data1_real, grid_table
  use grid1D_m,       only: grid1D
  use res_equation_m, only: res_equation
  use stencil_m,      only: dirichlet_stencil
  use test_case_m,    only: test_case
  use variable_m,     only: variable
  use vselector_m,    only: vselector

  implicit none

  private
  public test_esystem

  type, extends(variable) :: var
    type(grid_data1_real), pointer :: d => null()
  contains
    procedure :: init => var_init
  end type

  type, extends(equation) :: equation1
    !! y = 3x+1

    type(vselector) :: x
      !! dependent
    type(vselector) :: y
      !! provided

    type(jacobian), pointer :: dydx => null()

    real, allocatable :: y_tmp(:)

  contains
    procedure :: init => equation1_init
    procedure :: eval => equation1_eval
  end type

  type, extends(res_equation) :: res_equation1
    !! residuum: f = y-x == 0

    type(vselector) :: x
      !! main variable
    type(vselector) :: y
      !! dependent

    type(jacobian), pointer :: dfdx => null()
    type(jacobian), pointer :: dfdy => null()

    real, allocatable :: f_tmp(:)

  contains
    procedure :: init => res_equation1_init
    procedure :: eval => res_equation1_eval
  end type

  type, extends(res_equation) :: res_equation2
    !! residuum: f = z*x+y == 0

    type(var), pointer :: z => null()
      !! main variable
    type(var), pointer :: x => null()
      !! dependent
    type(var), pointer :: y => null()
      !! dependent

    type(jacobian), pointer :: dfdx => null()
    type(jacobian), pointer :: dfdy => null()
    type(jacobian), pointer :: dfdz => null()
  contains
    procedure :: init => res_equation2_init
    procedure :: eval => res_equation2_eval
  end type

  type(grid1D),     target :: g
  type(grid_table), target :: gtab

contains

  subroutine test_esystem()
    integer         :: i
    type(test_case) :: tc

    type(var)           :: x, y, z
    type(equation1)     :: eq1
    type(esystem)       :: es
    type(res_equation1) :: req1
    type(res_equation2) :: req2


    call tc%init("esystem: dir solver")

    ! init grid
    call g%init(real([0, 1, 2, 4]))

    ! init vars
    call x%init("x", IDX_VERTEX)
    call y%init("y", IDX_VERTEX)
    call z%init("z", IDX_VERTEX)

    ! init by random data. otherwise: newton doesnt converge
    x%d%data = g%x+100
    y%d%data = -10*(g%x+1)
    z%d%data = 10.0+5*g%x

    ! setup grid table + ptr
    call gtab%init('all vertices', g, IDX_VERTEX, 0, initial_flags = .true.)
    call gtab%init_final()

    ! equation system
    call es%init('test eqs system')

    ! equation1
    call eq1%init(y, x)
    call es%add_equation(eq1)

    ! res eq
    call req1%init(x, y)
    call es%add_equation(req1)

    ! res eq
    call req2%init(z, x, y)
    call es%add_equation(req2)

    ! finish equation system
    call es%init_final()

    ! solve esystem
    call es%solve()

    call tc%assert_eq([(-0.5, i=1, 4)], x%d%data, 1e-14, "es: solve: result x" )
    call tc%assert_eq([(-1.0, i=1, 4)], z%d%data, 1e-14, "es: solve: result z" )

    call tc%finish()
  end subroutine

  subroutine var_init(this, name, idx_type)
    class(var),   intent(out) :: this
    character(*), intent(in)  :: name
    integer,      intent(in)  :: idx_type

    select case (idx_type)
      case (IDX_VERTEX)
        call this%variable_init(name//'_vert', 'V', g, idx_type, 0)
      case (IDX_CELL)
        call this%variable_init(name//'_cell', 'V', g, idx_type, 0)
      case default
        call program_error('idx_type not defined for 1d var')
    end select

    this%d => this%data%get_ptr1()
  end subroutine

  subroutine equation1_init(this, y, x)
    class(equation1), target, intent(out) :: this
    type(var),        target, intent(in)  :: y
      !! prov var
    type(var),        target, intent(in)  :: x
      !! dependent var

    integer :: i, idx1(1), idx2(1), ix, iy

    ! init base
    call this%equation_init('eq1: y=x+1')

    ! variable selectors
    call this%x%init(x, gtab)
    call this%y%init(y, gtab)

    ! temp data
    allocate (this%y_tmp(this%y%n))

    ! provide y
    iy = this%provide(this%y)

    ! add dependency
    ix = this%depend(this%x)

    ! init jacobian
    this%dydx => this%init_jaco(iy, ix, const = .true.)
    do i = 1, gtab%n
      idx1 = gtab%get_idx(i)
      idx2 = idx1
      call this%dydx%set(idx1, idx2, 3.0)
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine equation1_eval(this)
    !! evaluate equation: y=3x+1
    class(equation1), intent(inout) :: this

    ! y_tmp <- jacobian * x
    call this%dydx%matr%mul_vec(this%x%get(), this%y_tmp)

    ! y: flat array to grid
    call this%y%set(this%y_tmp+1)
  end subroutine

  subroutine res_equation1_init(this, x, y)
    class(res_equation1), target, intent(out) :: this
    type(var),            target, intent(in)  :: x
      !! main var
    type(var),            target, intent(in)  :: y
      !! dependent var

    integer :: i, idx1(1), idx2(1), idep

    ! init base
    call this%equation_init('req1: f=y-x')

    ! variable selectors
    call this%x%init(x, gtab)
    call this%y%init(y, gtab)

    call this%init_f(this%x)

    allocate (this%f_tmp(this%f%n))

    ! add dependencies
    idep = this%depend(this%x)

    ! init jacobian for x: dfdx
    this%dfdx => this%init_jaco_f(idep, const = .true.)
    do i = 1, gtab%n
      idx1 = gtab%get_idx(i)
      idx2 = idx1
      call this%dfdx%set(idx1, idx2, -1.0)
    end do

    idep = this%depend(this%y)

    ! init jacobian for y: dfdy
    this%dfdy => this%init_jaco_f(idep, const = .true.)
    do i = 1, gtab%n
      idx1 = gtab%get_idx(i)
      idx2 = idx1
      call this%dfdy%set(idx1, idx2, 1.0)
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine res_equation1_eval(this)
    !! evaluate equation: f=y-x
    class(res_equation1), intent(inout) :: this

    ! y_tmp <- jacobian * x
    call this%dfdx%matr%mul_vec(this%x%get(), this%f_tmp)
    call this%dfdy%matr%mul_vec(this%y%get(), this%f_tmp, fact_y=1.0)

    ! y: flat array to grid
    call this%f%set(this%f_tmp)
  end subroutine

  subroutine res_equation2_init(this, z, x, y)
    class(res_equation2), target, intent(out) :: this
    type(var),            target, intent(in)  :: z
      !! main var
    type(var),            target, intent(in)  :: x
      !! dependent var
    type(var),            target, intent(in)  :: y
      !! dependent var

    integer :: i, idx(1), idep

    ! init base
    call this%equation_init('req1: f=z*x+y')

    ! save variables
    this%x => x
    this%y => y
    this%z => z

    ! setting main var
    call this%init_f(z, gtab)

    ! add main: z
    idep = this%depend(z, gtab)
    this%dfdz => this%init_jaco_f(idep)

    ! add dep: x
    idep = this%depend(x, gtab)
    this%dfdx => this%init_jaco_f(idep)

    ! add dep: y
    idep = this%depend(y, gtab)
    this%dfdy => this%init_jaco_f(idep, const = .true.)
    do i = 1, gtab%n
      idx = gtab%get_idx(i)
      call this%dfdy%set(idx, idx, 1.0)
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine res_equation2_eval(this)
    !! evaluate equation: f=z*x+y
    class(res_equation2), intent(inout) :: this

    integer :: i, idx(1)
    real    :: x, y, z

    do i = 1, gtab%n
      idx = gtab%get_idx(i)

      x = this%x%get(idx)
      y = this%y%get(idx)
      z = this%z%get(idx)

      call this%f%set(idx, [z*x+y])

      call this%dfdx%set(idx, idx, z)
      call this%dfdz%set(idx, idx, x)
    end do
  end subroutine

end module
