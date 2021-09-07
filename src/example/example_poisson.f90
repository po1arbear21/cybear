module example_poisson_m

  use example_charge_density_m, only: charge_dens
  use example_contact_m,        only: contacts, uncontacted, grd_contacts
  use example_device_m,         only: adj_v, eps, grd
  use example_potential_m,      only: pot
  use grid_m,                   only: IDX_VERTEX
  use grid0D_m,                 only: get_dummy_grid
  use jacobian_m,               only: jacobian, jacobian_ptr
  use math_m,                   only: eye_real
  use res_equation_m,           only: res_equation
  use stencil_m,                only: dirichlet_stencil, empty_stencil, near_neighb_stencil
  use vselector_m,              only: vselector

  implicit none

  private
  public pois

  type, extends(res_equation) :: poisson
    !! laplace phi = -rho

    type(vselector)           :: pot
      !! main variable
    type(vselector)           :: rho
      !! dependency
    type(vselector)           :: volt
      !! dependency

    type(dirichlet_stencil)   :: st_dir
    type(dirichlet_stencil)   :: st_dir_volt
    type(empty_stencil)       :: st_em
    type(near_neighb_stencil) :: st_nn

    type(jacobian), pointer   :: jaco_pot  => null()
    type(jacobian), pointer   :: jaco_rho  => null()
    type(jacobian), pointer   :: jaco_volt => null()
  contains
    procedure :: init => poisson_init
    procedure :: eval => poisson_eval
  end type

  type(poisson) :: pois

contains
  subroutine poisson_init(this)
    class(poisson), intent(out) :: this

    integer           :: i, idx1(1), idx2(1), j, dum(0)
    real              :: cap
    real, allocatable :: d_volt(:,:), eye(:,:)

    allocate(d_volt(1,size(contacts)), eye(size(contacts), size(contacts)))
    eye    = eye_real(size(contacts))

    ! init res_equation
    call this%equation_init("poisson")

    ! vselect the variables
    call this%volt%init([(contacts(i)%volt%get_ptr(), i = 1 , size(contacts))], "voltages")
    call this%pot%init(pot, [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])
    call this%rho%init(charge_dens)

    ! set main variable
    call this%init_f(this%pot)

    ! init stencils
    call this%st_dir%init(grd)
    call this%st_dir_volt%init(grd, g2 = get_dummy_grid(), perm = dum)
    call this%st_nn%init( grd, IDX_VERTEX, 0, IDX_VERTEX, 0)
    call this%st_em%init()

    ! init jaco
    this%jaco_pot  => this%init_jaco_f(this%depend(this%pot),  [this%st_nn%get_ptr(),  (this%st_dir%get_ptr(),      i = 1, size(contacts))], const = .true.)
    this%jaco_rho  => this%init_jaco_f(this%depend(this%rho),  [this%st_dir%get_ptr(), (this%st_em%get_ptr(),       i = 1, size(contacts))], const = .true.)
    this%jaco_volt => this%init_jaco_f(this%depend(this%volt), [this%st_em%get_ptr(),  (this%st_dir_volt%get_ptr(), i = 1, size(contacts))], const = .true.)

    ! loop over cells
    do i = 1, size(grd%x)-1
      idx1 = [i]
      idx2 = [i+1]

      ! capacity
      cap = eps%get(idx1) / grd%get_len(idx1, 1)

      ! setting jaco_pot
      if (uncontacted%flags%get(idx1)) then
        call this%jaco_pot%add(idx1, idx1,  cap)
        call this%jaco_pot%add(idx1, idx2, -cap)
      end if

      if (uncontacted%flags%get(idx2)) then
        call this%jaco_pot%add(idx2, idx1, -cap)
        call this%jaco_pot%add(idx2, idx2,  cap)
      end if
    end do

    ! loop over vertices
    do i = 1, size(grd%x)
      idx1 = [i]

      if (uncontacted%flags%get(idx1)) then
        ! setting jaco_rho
        call this%jaco_rho%set(idx1, idx1, -adj_v%get(idx1))
      else
        j = grd_contacts%get(idx1)
        d_volt = -eye(j:j, :)

        ! setting jaco_volt
        call this%jaco_volt%set(idx1, dum,  d_volt)

        ! setting jaco_pot
        call this%jaco_pot%set( idx1, idx1, 1.0)
      end if
    end do
    call this%init_final()
  end subroutine

  subroutine poisson_eval(this)
    class(poisson), intent(inout) :: this

    real, allocatable :: tmp(:)
    integer           :: i, j, idx(1)

    allocate(tmp(this%pot%n))

    ! calculating potential (f)
    call this%jaco_pot%matr%mul_vec(this%pot%get(), tmp)
    call this%jaco_rho%matr%mul_vec(this%rho%get(), tmp, fact_y = 1.0)
    call this%f%set(tmp)

    do i = 1, size(contacts)
      do j = 1, contacts(i)%conts%n
        idx = contacts(i)%conts%get_idx(j)
        call this%f%update(idx, [-contacts(i)%phi_ms-contacts(i)%volt%get()])
      end do
    end do
  end subroutine

end module
