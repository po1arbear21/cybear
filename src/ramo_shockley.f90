module ramo_shockley_m

  use charge_density_m,  only: charge_density
  use current_m,         only: current
  use current_density_m, only: current_density
  use device_params_m,   only: device_params
  use esystem_m,         only: esystem
  use grid_m,            only: IDX_VERTEX, IDX_EDGE, IDX_CELL
  use grid_data_m,       only: grid_data_real, allocate_grid_data1_real
  use grid0D_m,          only: get_dummy_grid
  use input_m,           only: input_section
  use jacobian_m,        only: jacobian, jacobian_ptr
  use math_m,            only: eye_real
  use normalization_m,   only: denorm
  use pardiso_m,         only: pardiso_default_params, pardiso_real
  use poisson_m,         only: poisson
  use potential_m,       only: potential
  use res_equation_m,    only: res_equation
  use semiconductor_m,   only: CR_CHARGE
  use sparse_m,          only: sparse_real
  use stencil_m,         only: dirichlet_stencil
  use voltage_m,         only: voltage
  use vselector_m,       only: vselector

  implicit none

  private
  public ramo_shockley, ramo_shockley_current

  type ramo_shockley
    !! Ramo-Shockley data object

    class(grid_data_real), allocatable :: x(:)
      !! fundamental solutions to laplace equation
    real,                  allocatable :: cap(:,:)
      !! Capacitance matrix
  contains
    procedure :: init => ramo_shockley_init
  end type

  type, extends(res_equation) :: ramo_shockley_current
    !! Ramo-Shockley current equation

    type(device_params), pointer :: par => null()

    type(vselector), allocatable :: cdens(:,:)
      !! electron/hole current density (edge direction, carrier index)
    type(vselector)              :: volt
      !! terminal voltage
    type(vselector)              :: curr
      !! terminal current

    type(dirichlet_stencil), allocatable :: st(:)
      !! coupling to cdens (edge direction)

    type(jacobian_ptr), allocatable :: jaco_cdens(:,:)
    type(jacobian),     pointer     :: jaco_volt => null()
    type(jacobian),     pointer     :: jaco_curr => null()
  contains
    procedure :: init => ramo_shockley_current_init
    procedure :: eval => ramo_shockley_current_eval
  end type

contains

  subroutine ramo_shockley_init(this, par, pot, rho, volt, poiss)
    !! initialize ramo-shockley data object (=> calculate fundamental solutions and capacitance matrix)
    class(ramo_shockley), intent(out)   :: this
    type(device_params),  intent(in)    :: par
      !! device parameters
    type(potential),      intent(inout) :: pot
      !! potential
    type(charge_density), intent(inout) :: rho
      !! charge density variable
    type(voltage),        intent(in)    :: volt(:)
      !! voltage variables
    type(poisson),        intent(in)    :: poiss
      !! poisson equation

    integer               :: i, j, k, nct, idx_dir
    integer, allocatable  :: idx1(:), idx2(:), idx(:)
    logical               :: status
    real                  :: cap, dxi, dxj
    real,    allocatable  :: x(:,:), rhs(:,:)
    type(esystem)         :: sys
    type(sparse_real)     :: df
    type(pardiso_real)    :: solver
    type(input_section)   :: params

    print "(A)", "ramo_shockley_init"

    ! init equation system
    call sys%init("ramo")

    ! add related equations to the system
    call sys%add_equation(poiss)

    ! provide variables
    call rho%reset()
    call sys%provide(rho, input = .false.)
    nct = par%nct
    do i = 1, nct
      call sys%provide(volt(i), input = .true.)
    end do

    ! finish initialization of equation system
    call sys%init_final()

    ! evaluate equation system to get jacobian
    call sys%eval()
    call sys%get_df(df)

    ! allocate memory to rhs and x
    allocate (rhs(df%nrows,sys%ninput), source = 0.0)
    allocate (  x(df%nrows,sys%ninput))

    ! set right-hand sides
    k = 0
    do i = 1, sys%input_equs%n
      do j = sys%input_i0(i), sys%input_i1(i)
        k = k + 1
        rhs(j,k) = 1.0
      end do
    end do

    ! factorize and solve
    params = pardiso_default_params()
    call solver%init(params)
    call solver%factorize(df)
    call solver%solve(rhs, x)
    call solver%destruct()

    ! save fundamental solutions
    call allocate_grid_data1_real(this%x, par%g%idx_dim, 1, nct)
    do i = 1, nct
      call this%x(i)%init(par%g, IDX_VERTEX, 0)

      ! save values in potential (correct sorting)
      call sys%set_x(x(:,i))

      ! extract values from potential
      call this%x(i)%set(pot%get())
    end do

    ! calculate capacitance matrix
    allocate (this%cap(nct,nct), source = 0.0)
    allocate (idx(par%g%idx_dim), idx1(par%g%idx_dim), idx2(par%g%idx_dim))
    do i = 1, nct
      do j = 1, i-1
        this%cap(i,j) = this%cap(j,i)
      end do
      do j = i, nct
        do idx_dir = 1, par%g%idx_dim
          do k = 1, par%poisson(IDX_EDGE,idx_dir)%n
            idx = par%poisson(IDX_EDGE,idx_dir)%get_idx(k)
            call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
            call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)
            cap = par%eps(IDX_EDGE,idx_dir)%get(idx) * par%curr_fact * par%surf(idx_dir)%get(idx) / par%g%get_len(idx, idx_dir)

            ! fundamental solution delta
            dxi = this%x(i)%get(idx2) - this%x(i)%get(idx1)
            dxj = this%x(j)%get(idx2) - this%x(j)%get(idx1)

            ! update capacitance matrix
            this%cap(i,j) = this%cap(i,j) + dxi * dxj * cap
          end do
        end do
      end do
    end do

    ! print capacitances
    do i = 1, nct
      do j = i, nct
        print "(A,ES25.16E3,A)", "C_" // par%contacts(i)%name // "_" // par%contacts(j)%name // " = ", denorm(this%cap(i,j), "F"), " F"
      end do
    end do
  end subroutine

  subroutine ramo_shockley_current_init(this, par, ramo, cdens, volt, curr)
    !! initialize Ramo-Shockley current equation
    class(ramo_shockley_current), intent(out) :: this
    type(device_params), target,  intent(in)  :: par
      !! device parameters
    type(ramo_shockley),          intent(in)  :: ramo
      !! Ramo-Shockley data object
    type(current_density),        intent(in)  :: cdens(:,:)
      !! electron/hole current density (edge direction, carrier index)
    type(voltage),                intent(in)  :: volt(:)
      !! terminal voltages
    type(current),                intent(in)  :: curr(:)
      !! terminal currents

    integer               :: i, idx_dir, ci, ict, dum(0)
    real, allocatable     :: d(:,:)
    integer, allocatable  :: idx1(:), idx2(:), idx(:), idx_bnd(:,:)
    logical               :: status

    print "(A)", "ramo_shockley_current_init"

    ! init base
    call this%equation_init("ramo_shockley")
    this%par => par

    ! create variable selectors
    allocate (this%cdens(par%g%idx_dim,2))
    do ci = par%ci0, par%ci1
      do idx_dir = 1, par%g%idx_dim
        call this%cdens(idx_dir,ci)%init(cdens(idx_dir,ci), par%transport(IDX_EDGE, idx_dir))
      end do
    end do
    call this%volt%init([(volt(ict)%get_ptr(), ict = 1, par%nct)], "voltages")
    call this%curr%init([(curr(ict)%get_ptr(), ict = 1, par%nct)], "currents")

    ! init residuals using this%curr as main variable
    call this%init_f(this%curr)

    ! init stencil
    allocate (idx_bnd(2, par%g%idx_dim))
    allocate (this%st(par%g%idx_dim))
    do idx_dir= 1, par%g%idx_dim
      call par%g%get_idx_bnd(IDX_EDGE, idx_dir, idx_bnd)
      call this%st(idx_dir)%init(get_dummy_grid(), g2 = par%g, perm = [(0, i = 1, par%g%idx_dim)], off1 = idx_bnd(1,:), off2 = idx_bnd(2,:))
    end do

    ! init jacobians
    allocate (this%jaco_cdens(par%g%idx_dim,2))
    do ci = par%ci0, par%ci1
      do idx_dir = 1, par%g%idx_dim
        this%jaco_cdens(idx_dir,ci)%p => this%init_jaco_f(this%depend(this%cdens(idx_dir,ci)), st = [this%st(idx_dir)%get_ptr()], const = .true.)
      end do
    end do
    this%jaco_volt => this%init_jaco_f(this%depend(this%volt), const = .true., dtime = .true.)
    this%jaco_curr => this%init_jaco_f(this%depend(this%curr), const = .true.)

    ! set jaco_cdens entries
    allocate (d(par%nct,1))
    allocate (idx(par%g%idx_dim), idx1(par%g%idx_dim), idx2(par%g%idx_dim))
    do idx_dir = 1, par%g%idx_dim
      do i = 1, par%transport(IDX_EDGE,idx_dir)%n
        idx = par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)
        do ict = 1, par%nct
          d(ict,1) = par%curr_fact * par%tr_surf(idx_dir)%get(idx) * (ramo%x(ict)%get(idx2) - ramo%x(ict)%get(idx1))
        end do
        do ci = par%ci0, par%ci1
          call this%jaco_cdens(idx_dir,ci)%p%set(dum, idx, CR_CHARGE(ci) * d)
        end do
      end do
    end do

    ! set jaco_volt entries
    call this%jaco_volt%set(dum, dum, - ramo%cap)

    ! set jaco_curr entries
    call this%jaco_curr%set(dum, dum, eye_real(par%nct))

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine ramo_shockley_current_eval(this)
    !! evaluate Ramo-Shockley current equation
    class(ramo_shockley_current), intent(inout) :: this

    integer           :: idx_dir, ci
    real, allocatable :: tmp(:)

    allocate (tmp(this%curr%n))

    ! calculate residuals
    call this%jaco_curr%matr%mul_vec(this%curr%get(), tmp)
    do ci = this%par%ci0, this%par%ci1
      do idx_dir = 1, this%par%g%idx_dim
        call this%jaco_cdens(idx_dir,ci)%p%matr%mul_vec(this%cdens(idx_dir,ci)%get(), tmp, fact_y = 1.0)
      end do
    end do
    call this%f%set(tmp)
  end subroutine

end module
