module example_mobility_m

  use example_device_m, only: grd, mu_0, beta_mob, v_sat
  use example_imref_m,  only: iref
  use equation_m,       only: equation
  use grid_m,           only: grid_data1_real, IDX_EDGE, IDX_VERTEX
  use jacobian_m,       only: jacobian, jacobian_ptr
  use stencil_m,        only: near_neighb_stencil
  use variable_m,       only: variable

  implicit none

  private
  public calc_mobil, mobil

  type, extends(variable) :: mobility
    !! electron mobility
    real, pointer :: x(:) => null()
  contains
    procedure :: init => mobility_init
  end type

  type, extends(equation) :: calc_mobility
    !! Caughey-Thomas mobility

    type(near_neighb_stencil) :: st

    type(jacobian), pointer   :: jaco_imref => null()
  contains
    procedure :: init => calc_mobility_init
    procedure :: eval => calc_mobility_eval
  end type

  type(mobility)      :: mobil
  type(calc_mobility) :: calc_mobil

contains

  subroutine mobility_init(this)
    class(mobility), intent(out) :: this

    type(grid_data1_real), pointer :: p

    call this%variable_init("mu", "cm^2/V/s", g = grd, idx_type = IDX_EDGE, idx_dir = 1)

    ! get pointer to data
    p      => this%data%get_ptr1()
    this%x => p%data
  end subroutine

  subroutine calc_mobility_init(this)
    class(calc_mobility), intent(out) :: this

    integer :: i_dep, i_prov

    ! init equation
    call this%equation_init("calc_mobility")

    ! init stencil
    call this%st%init(grd, IDX_EDGE, 1, IDX_VERTEX, 0)

    ! provides mobility
    i_prov = this%provide(mobil)
    ! depends on imref
    i_dep = this%depend(iref)
    ! init jaco
    this%jaco_imref => this%init_jaco(i_prov, i_dep, [this%st%get_ptr()])

    ! finalize equation
    call this%init_final()
  end subroutine


  subroutine calc_mobility_eval(this)
    class(calc_mobility), intent(inout) :: this

    integer :: abs_dif, i, idx1(1), idx2(1)
    real    :: dif, fact, val

    do i = 1, size(grd%x)-1
      idx1 = [i]
      idx2 = [i+1]

      fact = mu_0%get(idx1) / (v_sat * grd%get_len(idx1, 1))
      dif  = iref%get(idx2) - iref%get(idx1)
      ! derivative of the absolute value
      abs_dif = 1
      if (dif >= 0) abs_dif = -1

      val = -abs_dif * mu_0%get(idx1) * (fact)**beta_mob * abs(dif)**(beta_mob-1) / (1 + (fact * abs(dif))**beta_mob)**(1/beta_mob - 1)

      ! calculating the jacobi_matrix
      call this%jaco_imref%set(idx1, idx1,  val)
      call this%jaco_imref%set(idx1, idx2, -val)
      ! calculating the mobility
      call mobil%set(idx1, mu_0%get(idx1) / (1 + (fact * abs(dif))**beta_mob)**(1/beta_mob))
    end do
  end subroutine

end module
