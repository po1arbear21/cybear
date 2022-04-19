module test_esystem_prec_m

  use equation_m,     only: equation
  use error_m,        only: program_error
  use esystem_m,      only: esystem
  use jacobian_m,     only: jacobian
  use grid_m,         only: grid, IDX_VERTEX, IDX_EDGE, IDX_FACE, IDX_CELL
  use grid_data_m,    only: grid_data1_real
  use grid_table_m,   only: grid_table
  use grid1D_m,       only: grid1D
  use newton_m,       only: newton_opt
  use res_equation_m, only: res_equation
  use stencil_m,      only: dirichlet_stencil
  use test_case_m,    only: test_case
  use util_m,         only: int2str
  use variable_m,     only: variable_real
  use vselector_m,    only: vselector

  implicit none

  private
  public test_esystem_prec

  type, extends(variable_real) :: var
    type(grid_data1_real), pointer :: d => null()
  contains
    procedure :: init => var_init
  end type

  type, extends(equation) :: equation1
    !! [z_1]   [ 1/10 x_1          -1/5  y_2 -z_1           +1/40 y_1^2   -6]
    !! [z_2] = [          1/10 y_1                -z_2                    +7]
    !! [z_3]   [-1/5  x_1          +1/10 y_2           -z_3 -1/50 y_1 y_2 -8]

    type(vselector), pointer :: z => null()
      !! provided variable
    type(vselector), pointer :: x => null(), y => null()
      !! dependencies

    type(jacobian), pointer :: dzdx => null(), dzdy => null()
      !! jacobians
    type(jacobian), pointer :: dzdy_prec => null()
      !! preconditioner jacobian
  contains
    procedure :: init => equation1_init
    procedure :: eval => equation1_eval
  end type

  type, extends(res_equation) :: res_equation1
    !! residuum: f = 10 x_1 -2 y_2 -3 z_2 -1/10 y_1 z_1 +3== 0

    type(vselector), pointer :: x => null()
      !! main variable
    type(vselector), pointer :: y => null(), z => null()
      !! dependencies

    type(jacobian), pointer :: dfdx => null(), dfdy => null(), dfdz => null()
      !! exact jacobians
    type(jacobian), pointer :: dfdy_prec => null(), dfdz_prec => null()
      !! preconditioner jacobians
  contains
    procedure :: init => res_equation1_init
    procedure :: eval => res_equation1_eval
  end type

  type, extends(res_equation) :: res_equation2
    !! residuum: f = [      15 y_1 -   y_2 -z_1 -4 z_3 -1/20 x_1 y_2 -4]
    !!               [-2x_1        +20 y_2 -z_2        -1/30 y_2 z_3 +5]

    type(vselector), pointer :: y => null()
      !! main variable
    type(vselector), pointer :: x => null(), z => null()
      !! dependencies

    type(jacobian), pointer :: dfdx => null(), dfdy => null(), dfdz => null()
      !! exact jacobians
    type(jacobian), pointer :: dfdy_prec => null()
      !! preconditioner jacobian
  contains
    procedure :: init => res_equation2_init
    procedure :: eval => res_equation2_eval
  end type

  type(grid1D),     target :: g
  type(grid_table), target :: gtab

contains

  subroutine test_esystem_prec()
    integer         :: i
    type(test_case) :: tc

    type(var)           :: xvar, yvar(2), zvar(3)
    type(vselector)     :: xvsel, yvsel, zvsel
    type(equation1)     :: eq1
    type(esystem)       :: es
    type(res_equation1) :: req1
    type(res_equation2) :: req2
    type(newton_opt) :: nopt

    call tc%init("esystem: it solver")

    ! init grid
    call g%init("g", real([0, 1, 2, 4]))

    ! init vars
    call xvar%init("x", IDX_VERTEX)
    do i = 1, 3
      if (i < 3) call yvar(i)%init("y"//int2str(i), IDX_VERTEX)
      call zvar(i)%init("z"//int2str(i), IDX_VERTEX)
    end do

    ! init vsel
    call xvsel%init(xvar)
    call yvsel%init([(yvar(i)%get_ptr(), i = 1, 2)], 'y')
    call zvsel%init([(zvar(i)%get_ptr(), i = 1, 3)], 'z')

    ! setup grid table + ptr
    call gtab%init('all vertices', g, IDX_VERTEX, 0, initial_flags = .true.)
    call gtab%init_final()

    ! equation system
    call es%init('test eqs system', precon = .true.)

    ! eq1
    call eq1%init(zvsel, xvsel, yvsel)
    call es%add_equation(eq1)

    ! res eq1
    call req1%init(xvsel, yvsel, zvsel)
    call es%add_equation(req1)

    ! res eq2
    call req2%init(yvsel, xvsel, zvsel)
    call es%add_equation(req2)

    ! finish equation system
    call es%init_final()

    ! solve esystem
    call nopt%init(es%n, it_solver = .true.)
    call es%solve(nopt = nopt)

    ! check results (computed by matlab)
    call tc%assert_eq([( 1.917598528792069, i=1, 4)], xvar%d%data,    1e-15, "es: solve: result x" )
    call tc%assert_eq([(-2.278730533004949, i=1, 4)], yvar(1)%d%data, 1e-15, "es: solve: result y1")
    call tc%assert_eq([( 0.276521052725538, i=1, 4)], yvar(2)%d%data, 1e-15, "es: solve: result y2")

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

  subroutine equation1_init(this, z, x, y)
    class(equation1),        intent(out) :: this
    type(vselector), target, intent(in)  :: z
      !! main vselector
    type(vselector), target, intent(in)  :: x, y
      !! dependencies

    integer :: i, idx(1), iprov, idep

    ! init base
    call this%equation_init('eq1', precon = .true.)

    ! save data
    this%x => x
    this%y => y
    this%z => z

    ! provide z
    iprov = this%provide(z)

    ! vsel x: constant jacobian
    this%dzdx => this%init_jaco(iprov, this%depend(x), const = .true.)
    do i = 1, gtab%n
      idx = gtab%get_idx(i)
      call this%dzdx%set(idx, idx, reshape([0.1, 0.0, -0.2], [3, 1]))
    end do

    ! vsel y: variable jacobian, variable precon
    idep = this%depend(y)
    this%dzdy      => this%init_jaco(iprov, idep)
    this%dzdy_prec => this%init_jaco(iprov, idep, precon = .true.)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine equation1_eval(this)
    !! set z, and jacobians/preconditioners wrt y (wrt x is constant/already set/not preconditioned)
    class(equation1), intent(inout) :: this

    integer :: i, idx(1)
    real    :: x(1), y(2)

    do i = 1, gtab%n
      idx = gtab%get_idx(i)

      x = this%x%get(idx)
      y = this%y%get(idx)

      call this%z%set(idx, [ 0.1*x(1)          -0.2*y(2) +0.025*y(1)**2  -6, &
        &                             0.1*y(1)                           +7, &
        &                   -0.2*x(1)          +0.1*y(2) -0.02*y(1)*y(2) -8  ])

      call this%dzdy%set(     idx, idx, reshape([0.05*y(1), 0.1, -0.02*y(2), -0.2, 0.0, 0.1-0.02*y(1)     ], [3, 2]))
      call this%dzdy_prec%set(idx, idx, reshape([0.05*y(1), 0.1, -0.02*y(2), -0.2, 0.0,    -0.1* y(1)-y(2)], [3, 2]))
    end do
  end subroutine

  subroutine res_equation1_init(this, x, y, z)
    class(res_equation1),    intent(out) :: this
    type(vselector), target, intent(in)  :: x
      !! main vselector
    type(vselector), target, intent(in)  :: y, z
      !! dependencies

    integer :: i, idx(1), idep

    ! init base
    call this%equation_init('req1', precon = .true.)

    ! save data
    this%x => x
    this%y => y
    this%z => z

    ! init residual, set main var
    call this%init_f(this%x)

    ! vsel x: constant jacobian
    this%dfdx => this%init_jaco_f(this%depend(x), const = .true.)
    do i = 1, gtab%n
      idx = gtab%get_idx(i)
      call this%dfdx%set(idx, idx, 10.0)
    end do

    ! vsel y: variable jacobian, constant precon
    idep = this%depend(y)
    this%dfdy      => this%init_jaco_f(idep)
    this%dfdy_prec => this%init_jaco_f(idep, const = .true., precon = .true.)
    do i = 1, gtab%n
      idx = gtab%get_idx(i)
      call this%dfdy_prec%set(idx, idx, reshape([0.0, -2.0], [1, 2]))
    end do

    ! vsel z: variable jacobian, variable precon
    idep = this%depend(z)
    this%dfdz      => this%init_jaco_f(idep)
    this%dfdz_prec => this%init_jaco_f(idep, precon = .true.)

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine res_equation1_eval(this)
    !! set residual, and jacobians/preconditioners wrt y, z (wrt x is constant/already set/not preconditioned)
    class(res_equation1), intent(inout) :: this

    integer :: i, idx(1)
    real    :: x(1), y(2), z(3)

    do i = 1, gtab%n
      idx = gtab%get_idx(i)

      x = this%x%get(idx)
      y = this%y%get(idx)
      z = this%z%get(idx)

      call this%f%set(idx, [10*x(1) -2*y(2) -3*z(2) -0.1*y(1)*z(1) + 3])

      call this%dfdy%set(     idx, idx, reshape([-z(1)/10, -2.0     ], [1, 2]))
      call this%dfdz%set(     idx, idx, reshape([-y(1)/10, -3.0, 0.0], [1, 3]))
      call this%dfdz_prec%set(idx, idx, reshape([-y(1)/10, -1.0, 0.0], [1, 3]))
    end do
  end subroutine

  subroutine res_equation2_init(this, y, x, z)
    class(res_equation2),    intent(out) :: this
    type(vselector), target, intent(in)  :: y
      !! main vselector
    type(vselector), target, intent(in)  :: x, z
      !! dependencies

    integer :: idep

    ! init base
    call this%equation_init('req2', precon = .true.)

    ! save data
    this%x => x
    this%y => y
    this%z => z

    ! init residual, set main var
    call this%init_f(this%y)

    ! vsel x: variable jacobian
    this%dfdx => this%init_jaco_f(this%depend(x))

    ! vsel y: variable jacobian, variable precon
    idep = this%depend(y)
    this%dfdy      => this%init_jaco_f(idep)
    this%dfdy_prec => this%init_jaco_f(idep, precon = .true.)

    ! vsel z: variable jacobian
    this%dfdz => this%init_jaco_f(this%depend(z))

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine res_equation2_eval(this)
    !! set residuals, jacobians wrt x, y, z, and preconditioner wrt y
    class(res_equation2), intent(inout) :: this

    integer :: i, idx(1)
    real    :: x(1), y(2), z(3)

    do i = 1, gtab%n
      idx = gtab%get_idx(i)

      x = this%x%get(idx)
      y = this%y%get(idx)
      z = this%z%get(idx)

      ! residuum:      f = [        15 y_1  -   y_2  -z_1        -4 z_3  -1/20   x_1  z_2  -4]
      !                    [-2 x_1          +20 y_2        -z_2          -1/30   y_2  z_3  +5]
      call this%f%set(idx, [        15*y(1) -   y(2) -z(1)       -4*z(3) -0.05*  x(1)*z(2) -4, &
        &                   -2*x(1)         +20*y(2)       -z(2)         -1.0/30*y(2)*z(3) +5])

      call this%dfdx%set(idx, idx, reshape([-z(2)/20, -2.0], [2, 1]))

      call this%dfdy%set(     idx, idx, reshape([15.0, 0.0, -1.0, 20-z(3)/30], [2, 2]))
      call this%dfdy_prec%set(idx, idx, reshape([15.0, 0.0, -1.0, 15-z(3)   ], [2, 2]))

      call this%dfdz%set(idx, idx, reshape([-1.0, 0.0, -0.5*x(1), -1.0, -4.0, -y(2)/30], [2, 3]))
    end do
  end subroutine

end module
