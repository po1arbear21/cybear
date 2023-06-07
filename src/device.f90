module device_m

  use charge_density_m,  only: charge_density, calc_charge_density
  use continuity_m,      only: continuity
  use current_m,         only: current
  use current_density_m, only: current_density, calc_current_density
  use density_m,         only: density
  use device_params_m,   only: device_params
  use esystem_m,         only: esystem
  use grid_m,            only: IDX_VERTEX
  use imref_m,           only: imref, calc_imref, calc_density
  use input_m,           only: input_file
  use dopant_m,          only: ionization, calc_ionization
  use mobility_m,        only: mobility, calc_mobility
  use poisson_m,         only: poisson
  use potential_m,       only: potential
  use ramo_shockley_m,   only: ramo_shockley, ramo_shockley_current
  use semiconductor_m,   only: CR_NAME
  use voltage_m,         only: voltage
  use error_m,           only: assert_failed, program_error

  implicit none

  private
  public :: device, dev

  type device
    !! semiconductor device

    type(device_params) :: par
      !! device geometry + material parameters

    type(potential)                    :: pot
      !! electrostatic potential
    type(density)                      :: dens(2)
      !! electron/hole density (carrier index)
    type(ionization)                   :: ion(2)
      !! donor/acceptor ionization ratio (dopant index)
    type(current_density), allocatable :: cdens(:,:)
      !! electron/hole current density (direction, carrier index)
    type(imref)                        :: iref(2)
      !! electron/hole quasi-fermi potential (carrier index)
    type(mobility),        allocatable :: mob(:,:)
      !! electron/hole mobility (direction, carrier index)
    type(charge_density)               :: rho
      !! charge density
    type(voltage),         allocatable :: volt(:)
      !! terminal voltages
    type(current),         allocatable :: curr(:)
      !! terminal currents

    type(ramo_shockley) :: ramo
      !! Ramo-Shockley data object

    type(poisson)                           :: poiss
      !! poisson equation
    type(continuity)                        :: contin(2)
      !! electron/hole continuity equation (carrier index)
    type(ramo_shockley_current)             :: ramo_curr
      !! Ramo-Shockley current equation
    type(calc_imref)                        :: calc_iref(2)
      !! calculate electron/hole imref from potential and density (carrier index)
    type(calc_density)                      :: calc_dens(2)
      !! calculate electron/hole density from potential and imref (carrier index)
    type(calc_ionization)                   :: calc_ion(2)
      !! calculate stationary donor/acceptor ionization ratio from potential and imref (dopant index)
    type(calc_mobility),        allocatable :: calc_mob(:,:)
      !! calculate electron/hole mobility using Caughey-Thomas model (direction, carrier index)
    type(calc_charge_density)               :: calc_rho
      !! calculate charge density from electron/hole density and (ionized) doping concentrations
    type(calc_current_density), allocatable :: calc_cdens(:,:)
      !! calculate electron/hole current density by drift-diffusion model (direction, carrier index)

    type(esystem) :: sys_nlpe
      !! non-linear poisson equation system
    type(esystem) :: sys_dd(2)
      !! electron/hole drift-diffusion equation system
    type(esystem) :: sys_full
      !! full newton equation system
  contains
    procedure :: init     => device_init
    procedure :: destruct => device_destruct
  end type

  type(device) :: dev
    !! global device

contains

  subroutine device_init(this, filename)
    !! initialize device using device file
    class(device), intent(out) :: this
    character(*),  intent(in)  :: filename
      !! device file name

    integer          :: ci, idx_dir, ict
    type(input_file) :: file

    ! load device parameters
    call file%init(filename)
    call this%par%init(file)

    ! init variables
    allocate (this%cdens(this%par%g%idx_dim,2))
    allocate (this%mob(this%par%g%idx_dim,2))
    call this%pot%init(this%par)
    do ci = this%par%ci0, this%par%ci1
      call this%dens(ci)%init(this%par, ci)
      call this%ion(ci )%init(this%par, ci)
      call this%iref(ci)%init(this%par, ci)
      do idx_dir = 1, this%par%g%idx_dim
        call this%cdens(idx_dir,ci)%init(this%par, ci, idx_dir)
        call this%mob(idx_dir,ci)%init(this%par, ci, idx_dir)
      end do
    end do
    call this%rho%init(this%par)
    allocate (this%volt(this%par%nct))
    allocate (this%curr(this%par%nct))
    do ict = 1, this%par%nct
      call this%volt(ict)%init("V_"//this%par%contacts(ict)%name)
      call this%curr(ict)%init("I_"//this%par%contacts(ict)%name)
    end do

    ! init equations
    allocate (this%calc_cdens(this%par%g%idx_dim,2))
    allocate (this%calc_mob(this%par%g%idx_dim,2))
    call this%poiss%init(this%par, this%pot, this%rho, this%volt)
    call this%ramo%init(this%par, this%pot, this%rho, this%volt, this%poiss)
    call this%ramo_curr%init(this%par, this%ramo, this%cdens, this%volt, this%curr)
    do ci = this%par%ci0, this%par%ci1
      call this%contin(ci)%init(this%par, this%dens(ci), this%cdens(:,ci))
      call this%calc_iref(ci)%init(this%par, this%pot, this%dens(ci), this%iref(ci))
      call this%calc_dens(ci)%init(this%par, this%pot, this%dens(ci), this%iref(ci))
      call this%calc_ion( ci)%init(this%par, this%pot, this%ion(ci),  this%iref(ci))
      do idx_dir = 1, this%par%g%idx_dim
        call this%calc_mob(idx_dir,ci)%init(this%par, this%iref(ci), this%mob(idx_dir,ci))
        call this%calc_cdens(idx_dir,ci)%init(this%par, this%pot, this%dens(ci), this%cdens(idx_dir,ci), this%mob(idx_dir,ci))
      end do
    end do
    call this%calc_rho%init(this%par, this%dens, this%ion, this%rho)

    ! init non-linear poisson equation system
    call this%sys_nlpe%init("non-linear poisson")
    call this%sys_nlpe%add_equation(this%poiss)
    call this%sys_nlpe%add_equation(this%calc_rho)
    do ci = this%par%ci0, this%par%ci1
      call this%sys_nlpe%add_equation(this%calc_dens(ci))
      call this%sys_nlpe%add_equation(this%calc_ion(ci))
      call this%sys_nlpe%provide(this%iref(ci), this%par%transport(IDX_VERTEX,0))
    end do
    do ict = 1, this%par%nct
      call this%sys_nlpe%provide(this%volt(ict))
    end do
    call this%sys_nlpe%init_final()
    call this%sys_nlpe%g%output("nlpe")

    ! init drift-diffusion equation system
    do ci = this%par%ci0, this%par%ci1
      call this%sys_dd(ci)%init(CR_NAME(ci)//" drift-diffusion")
      call this%sys_dd(ci)%add_equation(this%contin(ci))
      do idx_dir = 1, this%par%g%idx_dim
        call this%sys_dd(ci)%add_equation(this%calc_cdens(idx_dir,ci))
        call this%sys_dd(ci)%add_equation(this%calc_mob(idx_dir,ci))
      end do
      call this%sys_dd(ci)%provide(this%iref(ci), this%par%transport(IDX_VERTEX,0))
      call this%sys_dd(ci)%provide(this%pot, this%par%transport(IDX_VERTEX,0))
      call this%sys_dd(ci)%init_final()
      call this%sys_dd(ci)%g%output(CR_NAME(ci)//"dd")
    end do

    ! init full-newton equation system
    call this%sys_full%init("full newton")
    call this%sys_full%add_equation(this%poiss)
    call this%sys_full%add_equation(this%calc_rho)
    do ci = this%par%ci0, this%par%ci1
      call this%sys_full%add_equation(this%contin(ci))
      do idx_dir = 1, this%par%g%idx_dim
        call this%sys_full%add_equation(this%calc_cdens(idx_dir,ci))
        call this%sys_full%add_equation(this%calc_mob(idx_dir,ci))
      end do
      call this%sys_full%add_equation(this%calc_iref(ci))
      call this%sys_full%add_equation(this%calc_ion(ci))
    end do
    call this%sys_full%add_equation(this%ramo_curr)
    do ict = 1, this%par%nct
      call this%sys_full%provide(this%volt(ict), input = .true.)
    end do

    call this%sys_full%init_final()
    call this%sys_full%g%output("full")
  end subroutine

  subroutine device_destruct(this)
    !! destruct device
    class(device), intent(inout) :: this

    call this%par%destruct()
  end subroutine

end module
