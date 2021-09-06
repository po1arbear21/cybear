module example_current_density_m

  use example_contact_m,   only: contacts, uncontacted
  use example_density_m,   only: dens
  use example_device_m,    only: grd
  use example_mobility_m,  only: mobil
  use example_potential_m, only: pot
  use equation_m,          only: equation
  use grid_m,              only: grid_data1_real, IDX_EDGE, IDX_VERTEX
  use jacobian_m,          only: jacobian, jacobian_ptr
  use math_m,              only: ber, dberdx
  use stencil_m,           only: dirichlet_stencil ,near_neighb_stencil
  use variable_m,          only: variable

  implicit none

  private
  public current_dens, calc_current_dens

  type, extends(variable) :: current_density
    !! electric current_density
    real, pointer :: x(:) => null()
  contains
    procedure :: init => current_density_init
  end type

  type, extends(equation) :: calc_current_density
    !! drift diffusion

    type(dirichlet_stencil)   :: st_dir
    type(near_neighb_stencil) :: st_nn

    type(jacobian), pointer   :: jaco_dens => null()
    type(jacobian), pointer   :: jaco_mob  => null()
    type(jacobian), pointer   :: jaco_pot  => null()
  contains
    procedure :: init => calc_current_density_init
    procedure :: eval => calc_current_density_eval
  end type

  type(current_density)      :: current_dens
  type(calc_current_density) :: calc_current_dens

contains

  subroutine current_density_init(this)
    class(current_density), intent(out) :: this

    type(grid_data1_real), pointer :: p => null()

    call this%variable_init("j", "1/cm^2/s", g = grd, idx_type = IDX_EDGE, idx_dir = 1)

    ! get pointer to data
    p      => this%data%get_ptr1()
    this%x => p%data
  end subroutine

  subroutine calc_current_density_init(this)
    class(calc_current_density), intent(out) :: this

    integer :: i_dep, i_prov, i

    ! init equation
    call this%equation_init("drift_diffusion_current_dens")

    ! init stencils
    call this%st_dir%init(grd)
    call this%st_nn%init( grd, IDX_EDGE, 1, IDX_VERTEX, 0)

    ! provides current_density
    i_prov = this%provide(current_dens)

    ! depends on density
    i_dep  = this%depend(dens, [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])
    this%jaco_dens => this%init_jaco(i_prov, i_dep, [this%st_nn%get_ptr()], const = .false.)

    ! depends on potential
    i_dep  = this%depend(pot, [uncontacted%get_ptr(), (contacts(i)%conts%get_ptr() , i=1, size(contacts))])
    this%jaco_pot => this%init_jaco(i_prov, i_dep, [this%st_nn%get_ptr()], const = .false.)

    ! depends on mobility
    i_dep = this%depend(mobil)
    this%jaco_mob => this%init_jaco(i_prov, i_dep, [this%st_dir%get_ptr()], const = .false.)

    call this%init_final()
  end subroutine

  subroutine calc_current_density_eval(this)
    class(calc_current_density), intent(inout) :: this

    integer :: i, idx1(1), idx2(1)
    real    :: ber1, ber2, dber1, dber2, dens1, dens2, len, mob

    ! loop over edges
    do i = 1, size(grd%x)-1
      idx1 = [i]
      idx2 = [i+1]

      len = grd%get_len(idx1, 1)

      mob = mobil%get(idx1)

      dens1 = dens%get(idx1)
      dens2 = dens%get(idx2)

      ber1 = ber(pot%get(idx2)-pot%get(idx1))
      ber2 = ber(pot%get(idx1)-pot%get(idx2))

      dber1 = dberdx(pot%get(idx2)-pot%get(idx1))
      dber2 = dberdx(pot%get(idx1)-pot%get(idx2))

      ! curr
      call current_dens%set(   idx1,       -mob * (ber1 * dens2 - ber2 * dens1) / len)

      ! jaco_dens
      call this%jaco_dens%set( idx1, idx1,  mob * ber2 / len)
      call this%jaco_dens%set( idx1, idx2, -mob * ber1 / len)

      ! jaco_pot
      call this%jaco_pot%set(  idx1, idx1,  mob * (dber1 * dens2 + dber2 * dens1) / len)
      call this%jaco_pot%set(  idx1, idx2, -mob * (dber1 * dens2 + dber2 * dens1) / len)

      ! jaco_mob
      call this%jaco_mob%set(  idx1, idx1, -(ber1 * dens2 - ber2 * dens1) / len)
    end do
  end subroutine

end module
