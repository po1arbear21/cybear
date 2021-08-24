module example_current_density_m

  use example_density_m,   only: dens
  use example_device_m,    only: grd, mobility
  use example_potential_m, only: pot
  use equation_m,          only: equation
  use grid_m,              only: grid_data1_real, IDX_EDGE, IDX_VERTEX
  use jacobian_m,          only: jacobian, jacobian_ptr
  use math_m,              only: ber, dberdx
  use stencil_m,           only: near_neighb_stencil
  use variable_m,          only: variable

  implicit none

  private
  public curr, c_curr

  type, extends(variable) :: current_density
    !! electric current_density
    real, pointer :: x(:) => null()
  contains
    procedure :: init => current_density_init
  end type

  type, extends(equation) :: calc_curr_d
    !! drift diffusion

    type(near_neighb_stencil) :: st

    type(jacobian), pointer   :: jaco_dens => null()
    type(jacobian), pointer   :: jaco_pot => null()
  contains
    procedure :: init => calc_curr_d_init
    procedure :: eval => calc_curr_d_eval
  end type

  type(current_density) :: curr
  type(calc_curr_d)     :: c_curr

contains

  subroutine current_density_init(this)
    class(current_density), intent(out) :: this

    type(grid_data1_real), pointer :: p => null()

    call this%variable_init("j", "1/cm^2/s", g = grd, idx_type = IDX_EDGE, idx_dir = 1)

    ! get pointer to data
    p      => this%data%get_ptr1()
    this%x => p%data
  end subroutine

  subroutine calc_curr_d_init(this)
    class(calc_curr_d), intent(out) :: this

    integer :: i_dep, i_prov

    ! init equation
    call this%equation_init("drift_diffusion_curr")

    ! init stencil
    call this%st%init(grd, IDX_EDGE, 1, IDX_VERTEX, 0)

    ! provides curr
    i_prov = this%provide(curr)
    ! depends on dens
    i_dep  = this%depend(dens)
    this%jaco_dens => this%init_jaco(i_prov, i_dep, [this%st%get_ptr()], const = .false.)
    ! depends on pot
    i_dep  = this%depend(pot)
    this%jaco_pot => this%init_jaco(i_prov, i_dep, [this%st%get_ptr()], const = .false.)

    call this%init_final()
  end subroutine

  subroutine calc_curr_d_eval(this)
    class(calc_curr_d), intent(inout) :: this

    integer :: i, idx1(1), idx2(1)
    real    :: ber1, ber2, dber1, dber2

    ! loop over edges
    do i = 1, size(grd%x)-1
      idx1 = [i]
      idx2 = [i+1]

      ber1 = ber(pot%get(idx2)-pot%get(idx1))
      ber2 = ber(pot%get(idx1)-pot%get(idx2))

      dber1 = dberdx(pot%get(idx2)-pot%get(idx1))
      dber2 = dberdx(pot%get(idx1)-pot%get(idx2))

      ! curr
      call curr%set(           idx1,        mobility * (ber1 * dens%get(idx2)) - ber2 * dens%get(idx1) / grd%get_len(idx1, 1))
      ! jaco_dens
      call this%jaco_dens%set( idx1, idx1, -mobility * ber2 / grd%get_len(idx1, 1))
      call this%jaco_dens%set( idx1, idx2,  mobility * ber1 / grd%get_len(idx1, 1))
      ! jaco_pot
      call this%jaco_pot%set(  idx1, idx1, -mobility * (dber1 * dens%get(idx2) + dber2 * dens%get(idx1)) / grd%get_len(idx1, 1))
      call this%jaco_pot%add(  idx1, idx2,  mobility * (dber1 * dens%get(idx2) + dber2 * dens%get(idx1)) / grd%get_len(idx1, 1))
    end do
  end subroutine

end module
