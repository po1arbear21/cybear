module example_mobility_m

  use dual_m
  use example_contact_m, only: contacts, uncontacted
  use example_device_m,  only: grd, mu_0, beta_mob, v_sat
  use example_imref_m,   only: iref
  use equation_m,        only: equation
  use grid_m,            only: grid_data1_real, IDX_EDGE, IDX_VERTEX
  use jacobian_m,        only: jacobian, jacobian_ptr
  use stencil_m,         only: near_neighb_stencil
  use variable_m,        only: variable_real

  implicit none

  private
  public calc_mobil, mobil

  type, extends(variable_real) :: mobility
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
    this%jaco_imref => this%init_jaco(i_prov, i_dep, st = [this%st%get_ptr()])

    ! finalize equation
    call this%init_final()
  end subroutine

  subroutine calc_mobility_eval(this)
    class(calc_mobility), intent(inout) :: this

    integer      :: i, idx1(1), idx2(1)
    real         :: mob0
    type(dual_1) :: mob, delta_iref

    call delta_iref%init(1.0, 1)
    do i = 1, size(grd%x)-1
      idx1 = [i  ]
      idx2 = [i+1]

      ! low-field mobility
      mob0 = mu_0%get(idx1)

      ! imref difference
      delta_iref%x = iref%get(idx2) - iref%get(idx1)

      ! calculate mobility including derivative wrt delta_imref
      mob = mob0 / (1 + ((mob0 / (v_sat * grd%get_len(idx1, 1))) * abs(delta_iref))**beta_mob)**(1/beta_mob)

      ! calculating the mobility
      call mobil%set(idx1, mob%x)

      ! setting jaco_imref
      call this%jaco_imref%set(idx1, idx1, -mob%dx(1))
      call this%jaco_imref%set(idx1, idx2,  mob%dx(1))
    end do
  end subroutine

end module
