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
  public calc_dens, calc_iref, iref

  type, extends(variable) :: imref
  !! quasi-fermi-potential
    real, pointer :: x(:) => null()
  contains
    procedure :: init => imref_init
  end type

  type, extends(equation) :: calc_imref
  !! n_intrin * exp(pot - iref)

  type(dirichlet_stencil) :: st

  type(jacobian), pointer :: jaco_pot  => null()
  type(jacobian), pointer :: jaco_dens => null()
contains
  procedure :: init => calc_imref_init
  procedure :: eval => calc_imref_eval
end type

  type, extends(equation) :: calc_density
  !! iref = -log(dens / n_intrin) + pot

  type(dirichlet_stencil) :: st

  type(jacobian), pointer :: jaco_pot   => null()
  type(jacobian), pointer :: jaco_imref => null()
contains
  procedure :: init => calc_density_init
  procedure :: eval => calc_density_eval
end type

  type(calc_density) :: calc_dens
  type(calc_imref)   :: calc_iref
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

  subroutine calc_imref_init(this)
    class(calc_imref), intent(out) :: this

    integer :: i_dep, i_prov, i

    ! init equation
    call this%equation_init("imref_calc")

    ! init stencil
    call this%st%init(grd)

    ! provides imref
    i_prov = this%provide(iref)
    ! depends on potential
    i_dep  = this%depend(pot)
    this%jaco_pot  => this%init_jaco(i_prov, i_dep, [this%st%get_ptr()], const = .true.)
    ! depends also on density
    i_dep  = this%depend(dens)
    this%jaco_dens => this%init_jaco(i_prov, i_dep, [this%st%get_ptr()], const = .false.)

    ! set jaco_pot entries
    do i=1, size(grd%x)
      call this%jaco_pot%set([i], [i], 1.0)
    end do

    call this%init_final()
  end subroutine

  subroutine calc_imref_eval(this)
    class(calc_imref), intent(inout) :: this

    integer :: i

    do i = 1, size(grd%x)
      ! calculate imref
      call iref%set([i], -log(dens%get([i]) / n_intrin) + pot%get([i]))

      ! set jaco_dens entries
      call this%jaco_dens%set([i], [i], -1 / dens%get([i]))
    end do
  end subroutine

  subroutine calc_density_init(this)
    class(calc_density), intent(out) :: this

    integer :: i_dep, i_prov

    ! init equation
    call this%equation_init("density_calc")

    ! init stencil
    call this%st%init(grd)

    ! provides density
    i_prov = this%provide(dens)

    ! depends on potential
    i_dep = this%depend(pot)
    this%jaco_pot => this%init_jaco(i_prov, i_dep, [this%st%get_ptr()], const = .false.)

    ! depends also on imref
    i_dep = this%depend(iref)
    this%jaco_imref => this%init_jaco(i_prov, i_dep, [this%st%get_ptr()], const = .false.)

    call this%init_final()
  end subroutine

  subroutine calc_density_eval(this)
    class(calc_density), intent(inout) :: this

    integer :: i
    real    :: n

    do i = 1, size(grd%x)
      ! calculate and set density
      n = n_intrin*exp(pot%get([i])-iref%get([i]))
      call dens%set([i], n)

      ! set jacobian entries
      call this%jaco_pot%set(  [i], [i],  n)
      call this%jaco_imref%set([i], [i], -n)
    end do
  end subroutine

end module
