module example_current_m

  use grid_m,              only: grid_data1_real, IDX_EDGE, IDX_VERTEX
  use example_density_m,   only: dens
  use example_device_m,    only: grd, adj_v
  use example_potential_m, only: pot
  use equation_m,          only: equation
  use jacobian_m,          only: jacobian, jacobian_ptr
  use math_m,              only: ber, dberdx
  use stencil_m,           only: near_neighb_stencil
  use variable_m,          only: variable

  implicit none

  private
  public curr, c_curr

  type, extends(variable) :: current
    !! electric current
    real, pointer :: x(:) => null()
  contains
    procedure :: init => current_init
  end type

  type, extends(equation) :: calc_curr
    !! drift diffusion
    type(near_neighb_stencil) :: st
    type(jacobian), pointer   :: jaco_dens => null()
    type(jacobian), pointer   :: jaco_adj_v => null()
    type(jacobian), pointer   :: jaco_pot => null()
  contains
    procedure :: init => calc_curr_init
    procedure :: eval => calc_curr_eval
  end type

  type(current) :: curr
  type(calc_curr) :: c_curr

contains

  subroutine current_init(this)
    class(current), intent(out) :: this

    type(grid_data1_real), pointer :: p => null()

    call this%variable_init("j", "C/s", g = grd, idx_type = IDX_EDGE, idx_dir = 1)

    ! get pointer to data
    p      => this%data%get_ptr1()
    this%x => p%data
  end subroutine

  subroutine calc_curr_init(this)
    class(calc_curr), intent(out) :: this

    integer :: i_dep, i_prov

    ! init equation
    call this%equation_init("drift_diffusion_curr")

    ! init stencil
    call this%st%init(grd, IDX_EDGE, 1, IDX_VERTEX, 0)

    ! provides curr
    i_prov = this%provide(curr)
    ! depends on dens
    i_dep  = this%depend(dens)
    this%jaco_dens => this%init_jaco(i_prov, idep, [this%st%get_ptr()], const = .false.)
    ! depends on adj_v
    i_dep  = this%depend(adj_v)
    this%jaco_adj_v => this%init_jaco(i_prov, idep, [this%st%get_ptr()], const = .false.)
    ! depends on pot
    i_dep  = this%depend(pot)
    this%jaco_pot => this%init_jaco(i_prov, idep, [this%st%get_ptr()], const = .false.)

    call this%init_final()
  end subroutine

  subroutine calc_curr_eval(this)
    class(calc_curr), intent(out) :: this

    integer :: i, idx1, idx2

    do i = 1, size(curr%x)
      idx1 = [i]
      idx2 = [i+1]

      ! jaco_dens
      call this%jaco_dens%set( idx1, idx1, -ber(pot%get(idx1)-pot%get(idx2)) / adj_v%get(idx1))
      call this%jaco_dens%add( idx2, idx1,  0.5 * ber(pot%get(idx2)-pot%get(idx1)) / adj_v%get(idx1))
      call this%jaco_dens%add( idx1, idx2,  0.5 * ber(pot%get(idx2)-pot%get(idx1)) / adj_v%get(idx1))
      ! jaco_adj_v
      call this%jaco_adj_v%set(idx1, idx1, -(ber(pot%get(idx2)-pot%get(idx1) * dens%get(idx2)) - ber(pot%get(idx1)-pot%get(idx2) * dens%get(idx1))) / adj_v%get(idx1)^2)
      ! jaco_pot
      call this%jaco_pot%set(  idx1, idx1, -(dberdx(pot%get(idx2)-pot%get(idx1) * dens%get(idx2)) + dberdx(pot%get(idx1)-pot%get(idx2) * dens%get(idx1))) / adj_v%get(idx1))
      call this%jaco_pot%add(  idx2, idx1,  0.5 * (dberdx(pot%get(idx2)-pot%get(idx1) * dens%get(idx2)) + dberdx(pot%get(idx1)-pot%get(idx2) * dens%get(idx1))) / adj_v%get(idx1))
      call this%jaco_pot%add(  idx1, idx2,  0.5 * (dberdx(pot%get(idx2)-pot%get(idx1) * dens%get(idx2)) + dberdx(pot%get(idx1)-pot%get(idx2) * dens%get(idx1))) / adj_v%get(idx1))
      ! curr
      call curr%set(           idx1,        (ber(pot%get(idx2)-pot%get(idx1) * dens%get(idx2)) - ber(pot%get(idx1)-pot%get([i+1]) * dens%get(idx1))) / adj_v%get(idx1))
    end do
  end subroutine

end module
