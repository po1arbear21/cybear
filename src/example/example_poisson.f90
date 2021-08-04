module example_poisson_m

  use example_contact_m,        only: uncontacted, contacts
  use example_charge_density_m, only: e_dens
  use example_device_m,         only: grd, eps, adj_v
  use example_potential_m,      only: potential, pot
  use jacobian_m,               only: jacobian, jacobian_ptr
  use res_equation_m,           only: res_equation
  use stencil_m,                only: dirichlet_stencil, near_neighb_stencil, empty_stencil
  use vselector_m,              only: vselector
  use grid_m,                   only: IDX_VERTEX

  implicit none

  private
  public pois

  type, extends(res_equation) :: poisson
    !! laplace phi = -rho
    type(vselector)           :: pot
    type(vselector)           :: rho
    type(dirichlet_stencil)   :: st_dir
    type(near_neighb_stencil) :: st_nn
    type(empty_stencil)       :: st_em
    type(jacobian), pointer   :: jaco_pot => null()
    type(jacobian), pointer   :: jaco_rho => null()
  contains
    procedure :: init => poisson_init
    procedure :: eval => poisson_eval
  end type

  type(poisson) :: pois

contains
  subroutine poisson_init(this)
    class(poisson), intent(out) :: this

    integer :: i, idx1(1), idx2(1)
    real    :: cap

    ! init res_equation
    call this%equation_init("poisson")

    ! vselect the potential
    call this%pot%init(pot,    [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])
    call this%rho%init(e_dens, [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])

    ! set main variable
    call this%init_f(this%pot)

    ! init stencils
    call this%st_dir%init(grd)
    call this%st_nn%init(grd, IDX_VERTEX, 0, IDX_VERTEX, 0)
    call this%st_em%init()


    ! init jaco
    this%jaco_pot => this%init_jaco_f(this%depend(this%pot), [this%st_nn%get_ptr(), (this%st_dir%get_ptr(), i=1, size(contacts))], const = .true.)
    this%jaco_rho => this%init_jaco_f(this%depend(this%rho), [this%st_dir%get_ptr(), (this%st_em%get_ptr(), i=1, size(contacts))], const = .true.)

    do i = 1, size(grd%x)-1
      idx1 = [i]
      idx2 = [i+1]

      ! capacity
      cap = eps%get(idx1) / grd%get_len(idx1, 1)

      if (uncontacted%flags%get(idx1)) then
        call this%jaco_pot%add(idx1, idx1,  cap)
        call this%jaco_pot%add(idx1, idx2, -cap)
      end if

      if (uncontacted%flags%get(idx2)) then
        call this%jaco_pot%add(idx2, idx1, -cap)
        call this%jaco_pot%add(idx2, idx2,  cap)
      end if
    end do

    do i = 1, size(grd%x)
      idx1 = [i]

      if (uncontacted%flags%get(idx1)) then
        call this%jaco_rho%set(idx1, idx1, -adj_v%get(idx1))
      else
        call this%jaco_pot%set(idx1, idx1, 1.0)
      end if
    end do
    call this%init_final()
  end subroutine


  subroutine poisson_eval(this)
    class(poisson), intent(inout) :: this

    real, allocatable :: tmp(:)
    integer           :: i, j, idx(1)

    allocate(tmp(this%pot%n))
    call this%jaco_pot%matr%mul_vec(this%pot%get(), tmp)
    call this%jaco_rho%matr%mul_vec(this%rho%get(), tmp, fact_y = 1.0)
    call this%f%set(tmp)

    do i = 1, size(contacts)
      do j = 1, contacts(i)%conts%n
        idx = contacts(i)%conts%get_idx(j)
        call this%f%update(idx, [-contacts(i)%phi_ms])
      end do
    end do
  end subroutine
end module
