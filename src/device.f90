module device_m

  use charge_density_m,  only: charge_density, calc_charge_density
  use contact_m,         only: CT_OHMIC
  use continuity_m,      only: continuity
  use current_m,         only: current
  use current_density_m, only: current_density, calc_current_density
  use density_m,         only: density
  use device_params_m,   only: device_params
  use esystem_m,         only: esystem
  use grid_m,            only: IDX_VERTEX
  use imref_m,           only: imref, calc_imref, calc_density
  use input_m,           only: input_file
  use ionization_m,      only: ionization, calc_ionization, ion_continuity, generation_recombination, calc_generation_recombination
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
      !! donor/acceptor ionization concentration (dopant index)
    type(generation_recombination)     :: genrec(2)
      !! netto recombination rate (generation - recombination)
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
    type(continuity)                        :: contin_stat(2)
      !! electron/hole continuity equation for stationary case (carrier index)
    type(ion_continuity)                    :: ion_contin(2)
      !! ionization continuity equations (carrier index)
    type(ramo_shockley_current)             :: ramo_curr
      !! Ramo-Shockley current equation
    type(calc_imref)                        :: calc_iref(2)
      !! calculate electron/hole imref from potential and density (carrier index)
    type(calc_density)                      :: calc_dens(2)
      !! calculate electron/hole density from potential and imref (carrier index)
    type(calc_ionization)                   :: calc_ion(2)
      !! calculate stationary donor/acceptor ionization ratio from potential and imref (dopant index)
    type(calc_generation_recombination)     :: calc_genrec(2)
      !! calculate generation-recombination
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
    type(esystem) :: sys_full_stat
      !! full newton equation system for steady-state only (uses imrefs instead of densities)
    type(esystem) :: sys_full
      !! full newton equation system
  contains
    procedure :: init        => device_init
    procedure :: destruct    => device_destruct
    procedure :: equilibrium => device_equilibrium
  end type

  type(device) :: dev
    !! global device

contains

  subroutine device_init(this, projectdir, filename, T)
    !! initialize device using device file
    class(device), intent(out) :: this
    character(*),  intent(in)  :: projectdir
      !! project directory
    character(*),  intent(in)  :: filename
      !! device file name
    real,          intent(in)  :: T
      !! temperature in K

    integer          :: ci, idx_dir, ict
    type(input_file) :: file

    ! load device parameters
    call file%init(filename)
    call this%par%init(file, projectdir, T)

    ! init variables
    allocate (this%cdens(this%par%g%idx_dim,2))
    allocate (this%mob(this%par%g%idx_dim,2))
    call this%pot%init(this%par)
    do ci = this%par%ci0, this%par%ci1
      call this%dens(ci)%init(this%par, ci)
      call this%ion( ci)%init(this%par, ci)
      call this%genrec( ci)%init(this%par, ci)
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
      call this%contin(     ci)%init(this%par, .false., this%dens(ci), this%iref(ci), this%cdens(:,ci), this%genrec(ci))
      call this%contin_stat(ci)%init(this%par,  .true., this%dens(ci), this%iref(ci), this%cdens(:,ci), this%genrec(ci))
      call this%calc_iref(ci)%init(this%par, this%pot, this%dens(ci), this%iref(ci))
      call this%calc_dens(ci)%init(this%par, this%pot, this%dens(ci), this%iref(ci))
      if (this%par%smc%incomp_ion) then
        call this%calc_ion(ci)%init(this%par, this%pot, this%ion(ci),  this%iref(ci))
        call this%calc_genrec(ci)%init(this%par, this%par%smc%genrec_tau(ci), this%genrec(ci), this%pot, this%iref(ci), this%ion(ci))
        call this%ion_contin(ci)%init(this%par, this%ion(ci), this%genrec(ci))
      end if
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
      if (this%par%smc%incomp_ion) then
        call this%sys_nlpe%add_equation(this%calc_ion(ci))
      end if
      call this%sys_nlpe%provide(this%iref(ci), this%par%transport(IDX_VERTEX,0))
    end do
    do ict = 1, this%par%nct
      call this%sys_nlpe%provide(this%volt(ict))
    end do
    call this%sys_nlpe%init_final()
    ! call this%sys_nlpe%g%output("nlpe")

    ! init drift-diffusion equation system
    do ci = this%par%ci0, this%par%ci1
      call this%sys_dd(ci)%init(CR_NAME(ci)//" drift-diffusion")
      call this%sys_dd(ci)%add_equation(this%contin(ci))
      do idx_dir = 1, this%par%g%idx_dim
        call this%sys_dd(ci)%add_equation(this%calc_cdens(idx_dir,ci))
        call this%sys_dd(ci)%add_equation(this%calc_mob(idx_dir,ci))
      end do
      if (this%par%smc%incomp_ion) then
        call this%sys_dd(ci)%add_equation(this%ion_contin(ci))
        call this%sys_dd(ci)%add_equation(this%calc_genrec(ci))
      end if
      call this%sys_dd(ci)%provide(this%iref(ci), this%par%transport(IDX_VERTEX,0))
      call this%sys_dd(ci)%provide(this%pot, this%par%transport(IDX_VERTEX,0))
      call this%sys_dd(ci)%init_final()
      call this%sys_dd(ci)%g%output(CR_NAME(ci)//"dd")
    end do

    ! init stationary full-newton equation system
    call this%sys_full_stat%init("full newton stat")
    call this%sys_full_stat%add_equation(this%poiss)
    call this%sys_full_stat%add_equation(this%calc_rho)
    do ci = this%par%ci0, this%par%ci1
      call this%sys_full_stat%add_equation(this%contin_stat(ci))
      call this%sys_full_stat%add_equation(this%calc_dens(ci))
      do idx_dir = 1, this%par%g%idx_dim
        call this%sys_full_stat%add_equation(this%calc_cdens(idx_dir,ci))
        call this%sys_full_stat%add_equation(this%calc_mob(idx_dir,ci))
      end do
      if (this%par%smc%incomp_ion) then
        call this%sys_full_stat%add_equation(this%ion_contin(ci))
        call this%sys_full_stat%add_equation(this%calc_genrec(ci))
      end if
    end do
    call this%sys_full_stat%add_equation(this%ramo_curr)
    do ict = 1, this%par%nct
      call this%sys_full_stat%provide(this%volt(ict), input = .true.)
    end do
    call this%sys_full_stat%init_final()
    call this%sys_full_stat%g%output("full_stat")

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
      if (this%par%smc%incomp_ion) then
        call this%sys_full%add_equation(this%ion_contin(ci))
        call this%sys_full%add_equation(this%calc_genrec(ci))
      end if
      call this%sys_full%add_equation(this%calc_iref(ci))
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

  function device_equilibrium(this) result(equi)
    !! check whether equilibrium voltage is applied
    class(device), intent(in) :: this
    logical                   :: equi
      !! return true if ohmic contacts are in equilibrium, otherwise false

    integer :: ict
    real    :: V0

    V0 = huge(1.0)
    equi = .true.
    do ict = 1, this%par%nct
      if (this%par%contacts(ict)%type /= CT_OHMIC) cycle
      if (V0 == huge(1.0)) then
        V0 = this%volt(ict)%x
        cycle
      end if
      if (this%volt(ict)%x /= V0) then
        equi = .false.
        exit
      end if
    end do
  end function

end module
