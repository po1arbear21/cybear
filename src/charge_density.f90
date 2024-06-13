module charge_density_m

  use density_m,       only: density
  use device_params_m, only: device_params
  use ionization_m,    only: ionization
  use equation_m,      only: equation
  use error_m,         only: assert_failed, program_error
  use grid_m,          only: IDX_VERTEX
  use grid_data_m,     only: grid_data1_real, grid_data2_real, grid_data3_real
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

    integer                 :: i, ci, iprov, idep
    type(jacobian), pointer :: jaco
    integer, allocatable    :: idx(:)

    print "(A)", "calc_charge_density_init"

    ! init equation
    call this%equation_init("calc_charge_density")
    this%par  => par
    this%dens => dens
    this%ion  => ion
    this%rho  => rho

    ! provides charge_density
    iprov = this%provide(rho, tab = par%transport_vct(0))

    ! depends on electron/hole density
    allocate (idx(par%g%idx_dim))
    do ci = par%ci0, par%ci1
      idep = this%depend(dens(ci), tab = par%transport_vct(0))
      jaco => this%init_jaco(iprov, idep, const = .true.)
      do i = 1, par%transport_vct(0)%n
        idx = par%transport_vct(0)%get_idx(i)
        call jaco%set(idx, idx, CR_CHARGE(ci))
      end do
    end do

    ! depends on ionization concentration (only defined if corresponding carrier density is enabled)
    if (par%smc%incomp_ion) then
      do ci = par%ci0, par%ci1
        idep = this%depend(ion(ci), tab = par%ionvert(ci))
        jaco => this%init_jaco(iprov, idep, const = .true.)
        do i = 1, par%ionvert(ci)%n
          idx = par%ionvert(ci)%get_idx(i)
          call jaco%set(idx, idx, DOP_CHARGE(ci))
        end do
      end do
    end if

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
    do i = 1, this%par%transport_vct(0)%n
      idx = this%par%transport_vct(0)%get_idx(i)

      ! electron/hole density
      rho = 0
      do ci = this%par%ci0, this%par%ci1
        rho = rho + CR_CHARGE(ci) * this%dens(ci)%get(idx)
      end do

      ! doping
      do ci = DOP_DCON, DOP_ACON
        if (this%par%dopvert(ci)%flags%get(idx)) then
          ! full ionization
          ion = this%par%dop(IDX_VERTEX,0,ci)%get(idx)

          ! incomplete ionization
          if (this%par%smc%incomp_ion .and. (ci >= this%par%ci0) .and. (ci <= this%par%ci1)) then
            if (this%par%ionvert(ci)%flags%get(idx)) then
              ion = this%ion(ci)%get(idx)
            end if
          end if

          rho = rho + DOP_CHARGE(ci) * ion
        end if
      end do

      ! save result
      call this%rho%set(idx, rho)
    end do
  end subroutine

end module
