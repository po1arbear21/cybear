module ramo_shockley_m

  use charge_density_m,  only: charge_density
  use current_m,         only: current
  use current_density_m, only: current_density
  use device_params_m,   only: device_params, CR_CHARGE
  use esystem_m,         only: esystem
  use grid_m,            only: grid_data2_real, IDX_VERTEX, IDX_EDGE, IDX_CELL
  use grid0D_m,          only: get_dummy_grid
  use jacobian_m,        only: jacobian, jacobian_ptr
  use math_m,            only: eye_real
  use matrix_m,          only: sparse_real
  use poisson_m,         only: poisson
  use potential_m,       only: potential
  use res_equation_m,    only: res_equation
  use stencil_m,         only: dirichlet_stencil
  use voltage_m,         only: voltage
  use vselector_m,       only: vselector

  implicit none

  private
  public ramo_shockley, ramo_shockley_current

  type ramo_shockley
    !! Ramo-Shockley data object

    type(grid_data2_real), allocatable :: x(:)
      !! fundamental solutions to laplace equation
    real,                  allocatable :: cap(:,:)
      !! Capacitance matrix
  contains
    procedure :: init => ramo_shockley_init
  end type

  type, extends(res_equation) :: ramo_shockley_current
    !! Ramo-Shockley current equation

    type(device_params), pointer :: par => null()

    type(vselector) :: cdens(2,2)
      !! electron/hole current density (edge direction, carrier index)
    type(vselector) :: volt
      !! terminal voltage
    type(vselector) :: curr
      !! terminal current

    type(dirichlet_stencil) :: st(2)
      !! coupling to cdens (edge direction)

    type(jacobian_ptr)      :: jaco_cdens(2,2)
    type(jacobian), pointer :: jaco_volt
    type(jacobian), pointer :: jaco_curr
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
    type(charge_density), intent(in)    :: rho
      !! charge density variable
    type(voltage),        intent(in)    :: volt(:)
      !! voltage variables
    type(poisson),        intent(in)    :: poiss
      !! poisson equation

    integer, parameter :: eye(2,2) = reshape([1, 0, 0, 1], [2, 2])

    integer           :: i, j, k, l, nct, idx1(2), idx2(2), jdx1(2), jdx2(2), idx_dir
    real              :: cap, surf, eps, len, dxi, dxj
    real, allocatable :: x(:,:), rhs(:,:)
    type(esystem)     :: sys
    type(sparse_real) :: df

    ! init equation system
    call sys%init("ramo")

    ! add related equations to the system
    call sys%add_equation(poiss)

    ! provide variables
    call sys%provide(rho, input = .false.)
    nct = size(par%contacts)
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
    call df%factorize()
    call df%solve_mat(rhs, x)
    call df%destruct()

    ! save fundamental solution
    allocate (this%x(nct))
    do i = 1, nct
      call this%x(i)%init(par%g, IDX_VERTEX, 0)

      ! save values in potential (correct sorting)
      call sys%set_x(x(:,i))

      ! extract values from potential
      call this%x(i)%set(pot%get())
    end do

    ! calculate capacitance matrix
    allocate (this%cap(nct,nct), source = 0.0)
    do i = 1, nct
      do j = 1, i-1
        this%cap(i,j) = this%cap(j,i)
      end do
      do j = i, nct
        do k = 1, par%poisson(IDX_CELL,0)%n
          idx1 = par%poisson(IDX_CELL,0)%get_idx(k)
          do idx_dir = 1, 2
            idx2 = idx1 + eye(:,idx_dir)

            ! capacitance
            surf = 0.5 * par%g%get_surf(idx1, idx_dir)
            eps  = par%eps%get(idx1)
            len  = par%g%get_len(idx1, idx_dir)
            cap  = surf * eps / len

            ! loop over the two edges in this direction
            do l = 0, 1
              ! first and second vertex indices
              jdx1 = idx1 + l * eye(:,3-idx_dir)
              jdx2 = idx2 + l * eye(:,3-idx_dir)

              ! fundamental solution delta
              dxi = this%x(i)%get(idx2) - this%x(i)%get(idx1)
              dxj = this%x(j)%get(idx2) - this%x(i)%get(idx1)

              ! update capacitance matrix
              this%cap(i,j) = this%cap(i,j) + dxi * dxj * cap
            end do
          end do
        end do
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
    type(current_density),        intent(in)  :: cdens(2,2)
      !! electron/hole current density (edge direction, carrier index)
    type(voltage),                intent(in)  :: volt(:)
      !! terminal voltages
    type(current),                intent(in)  :: curr(:)
      !! terminal currents

    integer, parameter :: eye(2,2) = reshape([1, 0, 0, 1], [2, 2])

    integer           :: i, idx1(2), idx2(2), idx_dir, ci, ict, dum(0)
    real, allocatable :: d(:,:)

    ! init base
    call this%equation_init("ramo_shockley")
    this%par => par

    ! create variable selectors
    do ci = par%ci0, par%ci1
      do idx_dir = 1, 2
        call this%cdens(idx_dir,ci)%init(cdens(idx_dir,ci), par%transport(IDX_EDGE, idx_dir))
      end do
    end do
    call this%volt%init([(volt(ict)%get_ptr(), ict = 1, size(volt))], "voltages")
    call this%curr%init([(curr(ict)%get_ptr(), ict = 1, size(curr))], "currents")

    ! init residuals using this%curr as main variable
    call this%init_f(this%curr)

    ! init stencil
    call this%st(1)%init(get_dummy_grid(), g2 = par%g, perm = [0, 0], off1 = [1, 1], off2 = [size(par%gx%x)-1,size(par%gy%x)  ])
    call this%st(2)%init(get_dummy_grid(), g2 = par%g, perm = [0, 0], off1 = [1, 1], off2 = [size(par%gx%x),  size(par%gy%x)-1])

    ! init jacobians
    do ci = par%ci0, par%ci1
      do idx_dir = 1, 2
        this%jaco_cdens(idx_dir,ci)%p => this%init_jaco_f(this%depend(this%cdens(idx_dir,ci)), st = [this%st(idx_dir)%get_ptr()], const = .true.)
      end do
    end do
    this%jaco_volt => this%init_jaco_f(this%depend(this%volt), const = .true., dtime = .true.)
    this%jaco_curr => this%init_jaco_f(this%depend(this%curr), const = .true.)

    ! set jaco_cdens entries
    allocate (d(size(par%contacts),1))
    do idx_dir = 1, 2
      do i = 1, par%transport(IDX_EDGE,idx_dir)%n
        idx1 = par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        idx2 = idx1 + eye(:,idx_dir)
        do ict = 1, size(par%contacts)
          d(ict,1) = par%tr_surf(idx_dir)%get(idx1) * (ramo%x(ict)%get(idx2) - ramo%x(ict)%get(idx1))
        end do
        do ci = par%ci0, par%ci1
          call this%jaco_cdens(idx_dir,ci)%p%set(dum, idx1, CR_CHARGE(ci) * d)
        end do
      end do
    end do

    ! set jaco_volt entries
    call this%jaco_volt%set(dum, dum, - ramo%cap)

    ! set jaco_curr entries
    call this%jaco_curr%set(dum, dum, eye_real(size(par%contacts)))

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
      do idx_dir = 1, 2
        call this%jaco_cdens(idx_dir,ci)%p%matr%mul_vec(this%cdens(idx_dir,ci)%get(), tmp, fact_y = 1.0)
      end do
    end do
    call this%f%set(tmp)
  end subroutine

end module
