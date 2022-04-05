module charge_density_m

  use density_m,       only: density
  use device_params_m, only: CR_CHARGE, CR_ELEC, CR_HOLE, device_params
  use equation_m,      only: equation
  use grid_m,          only: grid_data2_real, IDX_VERTEX
  use jacobian_m,      only: jacobian
  use variable_m,      only: variable_real

  implicit none

  private
  public charge_density, calc_charge_density

  type, extends(variable_real) :: charge_density
    !! charge density
    real, pointer :: x(:,:) => null()
      !! direct pointer to data for easy access
  contains
    procedure :: init => charge_density_init
  end type

  type, extends(equation) :: calc_charge_density
    !! calculate charge density from electron/hole density and doping concentration
    type(device_params),  pointer :: par     => null()
      !! pointer to device parameters
    type(density),        pointer :: dens(:) => null()
      !! pointer to electron/hole density
    type(charge_density), pointer :: rho     => null()
      !! pointer to charge density
  contains
    procedure :: init => calc_charge_density_init
    procedure :: eval => calc_charge_density_eval
  end type

contains

  subroutine charge_density_init(this, par)
    !! initialize charge density
    class(charge_density), intent(out) :: this
    type(device_params),   intent(in)  :: par
      !! device parameters

    type(grid_data2_real), pointer :: p => null()

    ! init base
    call this%variable_init("rho", "C/cm^3", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)

    ! get pointer to data
    p      => this%data%get_ptr2()
    this%x => p%data
  end subroutine

  subroutine calc_charge_density_init(this, par, dens, rho)
    !! initialize charge density calculation equation
    class(calc_charge_density),   intent(out) :: this
    type(device_params),  target, intent(in)  :: par
      !! device parameters
    type(density),        target, intent(in)  :: dens(:)
      !! electron/hole density
    type(charge_density), target, intent(in)  :: rho
      !! charge density

    integer                 :: i, idx(2), ci, idep(2), iprov
    type(jacobian), pointer :: jaco

    ! init equation
    call this%equation_init("calc_charge_density")
    this%par  => par
    this%dens => dens
    this%rho  => rho

    ! provides charge_density
    iprov = this%provide(rho, tab = par%transport(IDX_VERTEX,0))

    ! depends on electron/hole density
    do ci = par%ci0, par%ci1
      idep(ci) = this%depend(dens(ci), tab = par%transport(IDX_VERTEX,0))

      ! init jacobian
      jaco => this%init_jaco(iprov, idep(ci), const = .true.)

      ! set jacobian entries
      do i = 1, par%transport(IDX_VERTEX,0)%n
        idx = par%transport(IDX_VERTEX,0)%get_idx(i)
        call jaco%set(idx, idx, CR_CHARGE(ci))
      end do
    end do

    call this%init_final()
  end subroutine

  subroutine calc_charge_density_eval(this)
    !! evaluate charge density calculation equation
    class(calc_charge_density), intent(inout) :: this

    integer :: i, idx(2), ci
    real    :: rho

    ! calculate charge density at each vertex in the transport region
    do i = 1, this%par%transport(IDX_VERTEX,0)%n
      idx = this%par%transport(IDX_VERTEX,0)%get_idx(i)

      ! get charge density
      rho = 0
      do ci = this%par%ci0, this%par%ci1
        rho = rho + CR_CHARGE(ci) * (this%dens(ci)%get(idx) - this%par%dop(IDX_VERTEX,0,ci)%get(idx))
      end do
      call this%rho%set(idx, rho)
    end do
  end subroutine

end module
