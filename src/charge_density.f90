module charge_density_m

  use density_m,       only: density
  use device_params_m, only: device_params
  use equation_m,      only: equation
  use error_m,         only: assert_failed, program_error
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
  use dopant_m,        only: ionization
  use jacobian_m,      only: jacobian
  use semiconductor_m, only: CR_CHARGE, DOP_DCON, DOP_ACON, DOP_CHARGE
  use variable_m,      only: variable_real

  implicit none

  private
  public charge_density, calc_charge_density

  type, extends(variable_real) :: charge_density
    !! charge density including electrons, holes and ionized dopants

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => charge_density_init
  end type

  type, extends(equation) :: calc_charge_density
    !! calculate charge density from electron/hole density and doping concentration
    type(device_params),  pointer :: par     => null()
      !! pointer to device parameters
    type(density),        pointer :: dens(:) => null()
      !! pointer to electron/hole density
    type(ionization),     pointer :: ion(:)  => null()
      !! pointer to donor/acceptor ionization ratio
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

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base
    call this%variable_init("rho", "C/cm^3", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)

    ! get pointer to data
    select case (par%g%idx_dim)
    case (1)
      p1 => this%data%get_ptr1()
      this%x1 => p1%data
    case (2)
      p2 => this%data%get_ptr2()
      this%x2 => p2%data
    case (3)
      p3 => this%data%get_ptr3()
      this%x3 => p3%data
    case default
      call program_error("Maximal 3 dimensions allowed")
    end select
  end subroutine

  subroutine calc_charge_density_init(this, par, dens, ion, rho)
    !! initialize charge density calculation equation
    class(calc_charge_density),   intent(out) :: this
    type(device_params),  target, intent(in)  :: par
      !! device parameters
    type(density),        target, intent(in)  :: dens(:)
      !! electron/hole density
    type(ionization),     target, intent(in)  :: ion(:)
      !! donor/acceptor ionization ratio
    type(charge_density), target, intent(in)  :: rho
      !! charge density

    integer                 :: i, ci, iprov, idep(2)
    type(jacobian), pointer :: jaco_dens, jaco_ion
    integer, allocatable    :: idx(:)

    ! init equation
    call this%equation_init("calc_charge_density")
    this%par  => par
    this%dens => dens
    this%ion  => ion
    this%rho  => rho

    ! provides charge_density
    iprov = this%provide(rho, tab = par%transport(IDX_VERTEX,0))

    ! depends on electron/hole density
    allocate (idx(par%g%idx_dim))
    do ci = par%ci0, par%ci1
      ! carrier density jaco
      idep(ci) = this%depend(dens(ci), tab = par%transport(IDX_VERTEX,0))
      jaco_dens => this%init_jaco(iprov, idep(ci), const = .true.)

      ! ionization jaco
      idep(ci) = this%depend(ion(ci), tab = par%transport(IDX_VERTEX,0))
      jaco_ion => this%init_jaco(iprov, idep(ci), const = .true.)

      ! set jacobian entries
      do i = 1, par%transport(IDX_VERTEX,0)%n
        idx = par%transport(IDX_VERTEX,0)%get_idx(i)
        call jaco_dens%set(idx, idx, CR_CHARGE(ci))
        call jaco_ion%set (idx, idx, DOP_CHARGE(ci)*this%par%dop(IDX_VERTEX,0,ci)%get(idx))
      end do
    end do

    call this%init_final()
  end subroutine

  subroutine calc_charge_density_eval(this)
    !! evaluate charge density calculation equation
    class(calc_charge_density), intent(inout) :: this

    integer              :: i, ci
    real                 :: rho, ion
    integer, allocatable :: idx(:)

    ! calculate charge density at each vertex in the transport region
    allocate (idx(this%par%g%idx_dim))
    do i = 1, this%par%transport(IDX_VERTEX,0)%n
      idx = this%par%transport(IDX_VERTEX,0)%get_idx(i)

      ! get charge density
      rho = 0
      do ci = this%par%ci0, this%par%ci1
        rho = rho + CR_CHARGE(ci) * this%dens(ci)%get(idx)
      end do
      do ci = DOP_DCON, DOP_ACON
        ! Ionization of dopants is only defined if the carrier type is simulated (iref dependency)
        ! 1.0 matches the old behavior
        if (this%par%ci0 <= ci .and. this%par%ci1 >= ci) then
          ion = this%ion(ci)%get(idx)
        else
          ion = 1.0
        end if

        rho = rho + DOP_CHARGE(ci) * ion * this%par%dop(IDX_VERTEX,0,ci)%get(idx)
      end do
      call this%rho%set(idx, rho)
    end do
  end subroutine

end module
