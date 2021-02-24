#include "../../util/macro.f90.inc"

module test_esystem_m

  use error_m
  use test_case_m

  use equation_m,       only: equation
  use esystem_m,        only: esystem
  use jacobian_m,       only: jacobian
  use grid_data_m,      only: grid_data1_real
  use grid_m,           only: IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL
  use grid_table_m,     only: grid_table
  use grid1D_m,         only: grid1D
  use res_equation_m,   only: res_equation
  use stencil_m,        only: dirichlet_stencil
  use variable_m,       only: variable
  use vselector_m,      only: vselector

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

    type(dirichlet_stencil) :: st
      !! y=y(x_i) only depends on depending variabls from same index (certex index i)

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

    type(dirichlet_stencil) :: st
      !! y=y(x_i) only depends on depending variabls from same index (certex index i) fixme

    real, allocatable :: f_tmp(:)

  contains
    procedure :: init => res_equation1_init
    procedure :: eval => res_equation1_eval
  end type

  type, extends(res_equation) :: res_equation2
    !! residuum: f = z*x+y == 0

    type(vselector) :: z
      !! main variable
    type(vselector) :: x
      !! dependent
    type(vselector) :: y
      !! dependent

    type(jacobian), pointer :: dfdx => null()
    type(jacobian), pointer :: dfdy => null()
    type(jacobian), pointer :: dfdz => null()

    type(dirichlet_stencil) :: st
      !! y=y(x_i) only depends on depending variabls from same index (certex index i) fixme
  contains
    procedure :: init => res_equation2_init
    procedure :: eval => res_equation2_eval
  end type

  type(grid1D),     target :: g
  type(grid_table), target :: gtab

contains

  subroutine test_esystem()
    type(test_case) :: tc

    type(var) :: x, y, z

    type(equation1)     :: eq1
    type(esystem)       :: es
    type(res_equation1) :: req1
    type(res_equation2) :: req2

    integer :: i

    print "(A)", "test_esystem"
    call tc%init("esystem")

    ! init grid
    call g%init(real([0, 1, 2, 4]))

    ! init vars
    call x%init(IDX_VERTEX)
    call y%init(IDX_VERTEX)
    call z%init(IDX_VERTEX)

    x%d%data = g%x+1
    y%d%data = 0.1*(g%x+1)
    z%d%data = 1.0+0.1*g%x

    ! setup grid table + ptr
    call gtab%init('all vertices', g, IDX_VERTEX, 0)
    call gtab%flags%set([(.true., i=1, 4)])
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

    print *, x%d%data
    print *, z%d%data

    call tc%assert_eq([(-0.5, i=1, 4)], x%d%data, 1e-14, "es: solve: result x" )
    call tc%assert_eq([(-1.0, i=1, 4)], z%d%data, 1e-14, "es: solve: result z" )

    call tc%finish()
  end subroutine

  subroutine var_init(this, idx_type)
    class(var), intent(out) :: this
    integer,    intent(in)  :: idx_type

    select case (idx_type)
      case (IDX_VERTEX)
        call this%variable_init('var_vert', 'V', g, idx_type, 0)
      case (IDX_CELL)
        call this%variable_init('var_cell', 'V', g, idx_type, 0)
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

    integer :: i

    ! init base
    call this%equation_init('eq1: y=x+1')

    ! variable selectors
    call this%x%init(x, [gtab%get_ptr()])
    call this%y%init(y, [gtab%get_ptr()])

    ! temp data
    allocate (this%y_tmp(this%y%n))

    ! provide n
    call this%provide(this%y)

    ! stencil
    call this%st%init() ! dirichlet stencil

    ! add dependency
    call this%depend(this%x)

    ! init jacobian
    this%dydx => this%init_jaco(1, 1, [this%st%get_ptr()], const = .true.)
    do i = 1, gtab%n
      call this%dydx%set(1, i, 1, 3.0)
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

    integer :: i

    ! init base
    call this%equation_init('req1: f=y-x')

    ! variable selectors
    call this%x%init(x, [gtab%get_ptr()])
    call this%y%init(y, [gtab%get_ptr()])

    call this%init_f(this%x)

    allocate (this%f_tmp(this%f%n))

    ! stencil
    call this%st%init() ! dirichlet stencil

    ! add dependencies
    call this%depend(this%x)

    ! init jacobian for x: dfdx
    this%dfdx => this%init_jaco_f(this%vdep%n, [this%st%get_ptr()], const = .true.)
    do i = 1, gtab%n
      call this%dfdx%set(1, i, 1, -1.0)
    end do

    call this%depend(this%y)

    ! init jacobian for y: dfdy
    this%dfdy => this%init_jaco_f(this%vdep%n, [this%st%get_ptr()], const = .true.)
    do i = 1, gtab%n
      call this%dfdy%set(1, i, 1, 1.0)
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

    integer :: i

    ! init base
    call this%equation_init('req1: f=y-x')

    ! variable selectors
    call this%x%init(x, [gtab%get_ptr()])
    call this%y%init(y, [gtab%get_ptr()])
    call this%z%init(z, [gtab%get_ptr()])

    ! setting main var
    call this%init_f(this%z)

    ! stencil
    call this%st%init() ! dirichlet stencil

    ! add main: z
    call this%depend(this%z)
    this%dfdz => this%init_jaco_f(this%vdep%n, [this%st%get_ptr()])

    ! add dep: x
    call this%depend(this%x)
    this%dfdx => this%init_jaco_f(this%vdep%n, [this%st%get_ptr()])

    ! add dep: y
    call this%depend(this%y)
    this%dfdy => this%init_jaco_f(this%vdep%n, [this%st%get_ptr()], const = .true.)
    do i = 1, gtab%n
      call this%dfdy%set(1, i, 1, 1.0)
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine res_equation2_eval(this)
    !! evaluate equation: f=z*x+y
    class(res_equation2), intent(inout) :: this

    integer :: i, idx1(1)
    real    :: x(1), y(1), z(1)

    do i = 1, gtab%n
      idx1 = gtab%get_idx(i)

      x = this%x%get(idx1)
      y = this%y%get(idx1)
      z = this%z%get(idx1)

      call this%f%set(idx1, z*x+y)

      call this%dfdx%set(1, i, 1, z(1))
      call this%dfdz%set(1, i, 1, x(1))
    end do
  end subroutine

end module
