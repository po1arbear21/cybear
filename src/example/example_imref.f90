module example_imref_m

  use equation_m,          only: equation
  use example_contact_m,   only: contacts, uncontacted
  use example_density_m,   only: dens
  use example_device_m,    only: dop, grd, n_intrin
  use example_potential_m, only: pot
  use grid_m,              only: grid_data1_real, IDX_VERTEX
  use jacobian_m,          only: jacobian, jacobian_ptr
  use stencil_m,           only: dirichlet_stencil
  use variable_m,          only: variable

  implicit none

  private
  public calc_dens, iref

  type, extends(variable) :: imref
  !! quasi-fermi-potential
    real, pointer :: x(:) => null()
  contains
    procedure :: init => imref_init
    procedure :: calc => imref_calc
  end type

    type, extends(equation) :: calc_density
    !! n_intrin*exp(pot-iref)

    type(dirichlet_stencil) :: st

    type(jacobian), pointer :: jaco_pot   => null()
    type(jacobian), pointer :: jaco_imref => null()
  contains
    procedure :: init => calc_density_init
    procedure :: eval => calc_density_eval
  end type

  type(calc_density) :: calc_dens
  type(imref)        :: iref

contains

  subroutine imref_init(this)
    class(imref), intent(out) :: this

    type(grid_data1_real), pointer :: p

    call this%variable_init("imref", "V", g = grd, idx_type = IDX_VERTEX, idx_dir = 0)

    ! get pointer to data
    p      => this%data%get_ptr1()
    this%x => p%data
  end subroutine

  subroutine imref_calc(this)
    class(imref), intent(inout) :: this

    real, allocatable :: dens_arr(:)

    allocate(dens_arr(size(grd%x)))
    dens_arr = dens%get()
    call this%set(-log(dens_arr / n_intrin) + pot%get())
  end subroutine

  subroutine calc_density_init(this)
    class(calc_density), intent(out) :: this

    integer :: i_dep, i_prov, i

    ! init equation
    call this%equation_init("density_calc")

    ! init stencil
    call this%st%init(grd)

    ! provides density
    i_prov = this%provide(dens, [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])
    ! depends on potential
    i_dep  = this%depend(pot,   [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])
    this%jaco_pot => this%init_jaco(i_prov, i_dep, [(this%st%get_ptr(), i = 0, size(contacts))], const = .false.)
    ! depends also on imref
    i_dep  = this%depend(iref)
    this%jaco_imref => this%init_jaco(i_prov, i_dep, [(this%st%get_ptr(), i = 0, size(contacts))], const = .false.)

    call this%init_final()
  end subroutine

  subroutine calc_density_eval(this)
    class(calc_density), intent(inout) :: this

    integer :: i

    ! calculating density
    do i = 1, size(dens%x)
      call this%jaco_pot%set(  [i], [i],  n_intrin*exp(pot%get([i])-iref%get([i])))
      call this%jaco_imref%set([i], [i], -n_intrin*exp(pot%get([i])-iref%get([i])))
      call dens%set(           [i],       n_intrin*exp(pot%get([i])-iref%get([i])))
    end do
  end subroutine

end module
