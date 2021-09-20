module example_ramo_m

  use esystem_m,                 only: esystem
  use example_charge_density_m,  only: charge_dens
  use example_contact_m,         only: contacts, grd_contacts, uncontacted
  use example_current_density_m, only: current_dens
  use example_device_m,          only: eps, grd
  use example_poisson_m,         only: poisson
  use example_potential_m,       only: pot
  use grid_m,                    only: grid_data1_real
  use grid0D_m,                  only: get_dummy_grid
  use jacobian_m,                only: jacobian, jacobian_ptr
  use math_m,                    only: eye_real
  use matrix_m,                  only: sparse_real
  use res_equation_m,            only: res_equation
  use stencil_m,                 only: dirichlet_stencil
  use vselector_m,               only: vselector

  implicit none

  private
  public ramo_cap, ramo_nu, ramo_eq
  public ramo_init

  type, extends(res_equation) :: ramo_shockley
    !! Ramo-Shockley equation I_i + S(grad(nu_i) * J) dV - C_ij dV_j/dt = 0

    type(vselector)         :: curr
      !! main variable
    type(vselector)         :: curr_dens
      !! dependency
    type(vselector)         :: volt
      !! dependency

    type(dirichlet_stencil) :: st_dir

    type(jacobian), pointer :: jaco_curr
    type(jacobian), pointer :: jaco_curr_dens
    type(jacobian), pointer :: jaco_volt
  contains
    procedure :: init => ramo_shockley_init()
    procedure :: eval => ramo_shockley_eval()
  end type

  real,                  allocatable :: ramo_cap(:,:)
  type(grid_data1_real), allocatable :: ramo_nu(:)

contains

  subroutine ramo_init()
    !! calculating the fundamental solutions nu_i of the boundary-value problem
    integer           :: i, j, k
    real              :: cap
    real, allocatable :: rhs(:,:), x(:,:)
    type(esystem)     :: sys_ramo
    type(sparse_real) :: df

    ! init equation system
    call sys_ramo%init("fundamental solutions")

    ! add related equations to the system
    call sys_ramo%add_equation(possion)

    ! provide variables
    call sys_ramo%provide(charge_dens, input = .false.)
    do i = 1, size(contacts)
      call sys_full%provide(contacts(i)%volt, input = .true.)
    end do

    ! evaluate equation system
    call sys_ramo%eva()

    ! get jacobian
    call sys_ramo%get_df(df)

    ! set right-hand sides
    k = 0
    do i = 1, sys_ramo%input_equs%n
      do j = sys_ramo%input_i0(i), sys_ramo%input_i1(i)
        k = k + 1
        rhs(j,k) = 1.0
      end do
    end do

    ! factorize and solve esystem
    call df%factorize()
    call df%solve_mat(rhs, x)

    ! allocate memory to ramo_nu and ramo_cap
    allocate(ramo_cap(size(contacts), size(contacts)), ramo_nu(size(contacts)))

    ! set the result in the esystem
    ! and set the fundamental solution as the result from the poisson equation
    do i = 1, size(contacts)
      call sys_ramo%set_x(x(i))
      call ramo_nu(i)%set(pot%x)
    end do

    ! iteration over the contacts
    do i = 1, size(ramo_mu)
      do j = i, size(ramo_mu)
        cap = 0
        ! sum over the grid
        do k = 1, size(grd%x)-1
          cap = cap + (ramo_nu(i)%get([k]) - ramo_nu(i)%get([k+1])) * (ramo_nu(j)%get([k]) - ramo_nu(j)%get([k+1]))
          cap = cap * grd%get_len([k], idx_dir = 1) * eps%get([k])
        end do
        ramo_cap(i, j) = cap
        ramo_cap(i, j) = cap
      end do
    end do
  end subroutine

  subroutine ramo_shockley_init(this)
    class(ramo_shockley), intent(out) :: this

    integer           :: i, j, k, idx1(1), idx2(1), dum(0)
    real, allocatable :: d_curr_dens(:,:), d_volt(:,:)

    ! init res_equation
    call this%equation_init("ramo_shockley")

    ! vselect the variables
    call this%curr%init([(contacts(i)%curr%get_ptr(), i = 1 , size(contacts))], "currents")
    call this%volt%init([(contacts(i)%volt%get_ptr(), i = 1 , size(contacts))], "voltages")
    call this%curr_dens%init(current_dens, [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])

    ! set main variable
    call this%init_f(this%curr)

    ! init stencils
    call this%st_dir%init(g1 = get_dummy_grid(), g2 = grd, perm = dum, off2 = [size(grd%x)-1])

    ! init jacos
    this%jaco_curr      => this%init_jaco_f(this%depend(this%curr),       [this%st_nn%get_ptr(),  (this%st_dir%get_ptr(),      i = 1, size(contacts))], const = .true.)
    this%jaco_volt      => this%init_jaco_f(this%depend(this%volt),       [this%st_em%get_ptr(),  (this%st_dir_volt%get_ptr(), i = 1, size(contacts))], const = .true.)
    this%jaco_curr_dens => this%init_jaco_f(this%depend(this%curr_dens),  [this%st_dir%get_ptr(), (this%st_em%get_ptr(),       i = 1, size(contacts))], const = .true.)

    ! allocate memory to d_volt and d_curr_dens
    allocate(d_volt(1, size(contacts)), d_curr_dens(1, size(contacts)))

    ! loop over edges
    do i = 1, size(grd%x)-1
      idx1 = [i]
      idx2 = [i+1]

      ! setting jaco_curr_dens
      do j = 1, size(contacts)
        d_curr_dens(1,i) = ramo_nu(i)%get(idx2) - ramo_nu(i)%get(idx1)
      end do
      this%jaco_curr_dens%set(dum, idx1, d_curr_dens)

      ! setting jaco_volt
      k = grd_contacts%get(idx1)
      d_curr = -eye(k:k, :)

      call this%jaco_volt%set(idx1, dum,  d_volt)

      ! setting jaco_curr

    end do
    call this%init_final()
  end subroutine

  subroutine ramo_shockley_eval(this)
    class(ramo_shockley), intent(inout) :: this

    real, allocatable :: tmp(:)

    allocate(tmp(this%curr%n))

    ! calculating the current (f)
    call this%jaco_curr%matr%mul_vec(this%curr%get(),           tmp)
    call this%jaco_curr_dens%matr%mul_vec(this%curr_dens%get(), tmp, fact_y = 1.0)
    call this%jaco_volt%matr%mul_vec(this%volt%get(),           tmp, fact_y = 1.0)
    call this%f%set(tmp)
  end subroutine

end module
