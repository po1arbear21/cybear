module example_current_m

  use esystem_m,                only: esystem
  use example_charge_density_m, only: charge_dens
  use example_contact_m,        only: contacts
  use example_device_m,         only: eps, grd
  use example_poisson_m,        only: poisson
  use example_potential_m,      only: pot
  use grid_m,                   only: grid_data1_real
  use matrix_m,                 only: sparse_real
  use variable_m,               only: variable

  implicit none

  private
  public ramo_cap, ramo_nu
  public ramo_init

  type, extends(variable) :: current
    !! electric current
    real, pointer :: x => null()
  contains
    procedure :: init => current_init
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

end module
