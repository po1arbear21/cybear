module device_m

  use charge_density_m,  only: charge_density, calc_charge_density
  use contact_m,         only: CT_SCHOTTKY
  use continuity_m,      only: continuity
  use current_m,         only: current
  use current_density_m, only: current_density, calc_current_density
  use density_m,         only: density
  use device_params_m,   only: device_params
  use electric_field_m,  only: electric_field, electric_field_abs, calc_efield, calc_efield_abs
  use esystem_m,         only: esystem
  use grid_m,            only: IDX_VERTEX
  use imref_m,           only: chemical_pot, imref, calc_chemical_pot_dens, calc_chemical_pot_imref, calc_imref, calc_density
  use input_m,           only: input_file
  use ionization_m,      only: ionization, calc_ionization, ion_continuity, generation_recombination, calc_generation_recombination
  use mobility_m,        only: mobility, calc_mobility
  use poisson_m,         only: poisson
  use potential_m,       only: potential
  use ramo_shockley_m,   only: ramo_shockley, ramo_shockley_current
  use schottky_m,        only: schottky_bc, calc_schottky_bc
  use semiconductor_m,   only: CR_NAME
  use voltage_m,         only: voltage

  implicit none

  private
  public :: device, dev

  type device
    !! semiconductor device

    type(device_params) :: par
      !! device geometry + material parameters

    type(potential)                             :: pot
      !! electrostatic potential
    type(electric_field_abs)                    :: efield_abs
      !! absolute electric field strength
    type(electric_field),           allocatable :: efield(:)
      !! directional electric field strength (direction)
    type(density)                               :: dens(2)
      !! electron/hole density (carrier index)
    type(ionization)                            :: ion(2)
      !! donor/acceptor ionization concentration (dopant index)
    type(generation_recombination)              :: genrec(2)
      !! netto recombination rate (generation - recombination)
    type(current_density),          allocatable :: cdens(:,:)
      !! electron/hole current density (direction, carrier index)
    type(chemical_pot)                          :: eta(2)
      !! electron/hole chemical potential (carrier index)
    type(imref)                                 :: iref(2)
      !! electron/hole quasi-fermi potential (carrier index)
    type(mobility),                 allocatable :: mob(:,:)
      !! electron/hole mobility (direction, carrier index)
    type(charge_density)                        :: rho
      !! charge density
    type(voltage),                  allocatable :: volt(:)
      !! terminal voltages
    type(current),                  allocatable :: curr(:)
      !! terminal currents
    type(schottky_bc),     allocatable :: sbc(:,:)
      !! Schottky boundary current variable per (contact, carrier);
      !! initialized only for CT_SCHOTTKY contacts

    type(ramo_shockley) :: ramo
      !! Ramo-Shockley data object

    type(poisson)                                    :: poiss
      !! poisson equation
    type(continuity)                                 :: contin(2)
      !! electron/hole continuity equation (carrier index)
    type(ion_continuity)                             :: ion_contin(2)
      !! ionization continuity equations (dopant index)
    type(ramo_shockley_current)                      :: ramo_curr
      !! Ramo-Shockley current equation
    type(calc_chemical_pot_dens)                     :: calc_eta_dens(2)
      !! calculate electron/hole chemical potential from density (carrier index)
    type(calc_chemical_pot_imref)                    :: calc_eta_iref(2)
      !! calculate electron/hole chemical potential from imref (carrier index)
    type(calc_imref)                                 :: calc_iref(2)
      !! calculate electron/hole imref from potential and chemical potential (carrier index)
    type(calc_density)                               :: calc_dens(2)
      !! calculate electron/hole density from potential and imref (carrier index)
    type(calc_ionization)                            :: calc_ion(2)
      !! calculate stationary donor/acceptor ionization ratio from potential and imref (dopant index)
    type(calc_generation_recombination)              :: calc_genrec(2)
      !! calculate generation-recombination (dopant index)
    type(calc_mobility),                 allocatable :: calc_mob(:,:)
      !! calculate electron/hole mobility using Caughey-Thomas model (direction, carrier index)
    type(calc_charge_density)                        :: calc_rho
      !! calculate charge density from electron/hole density and (ionized) doping concentrations
    type(calc_current_density),          allocatable :: calc_cdens(:,:)
      !! calculate electron/hole current density by drift-diffusion model (direction, carrier index)
    type(calc_efield),          allocatable :: calc_efield(:)
      !! calculate electric field from potential gradient (direction)
    type(calc_efield_abs)                   :: calc_efield_abs
      !! calculate absolute electric field strength from components
    type(calc_schottky_bc),     allocatable :: calc_sbc(:,:)
      !! equation that fills sbc(ict,ci); allocated for all contacts but
      !! only init+registered for CT_SCHOTTKY

    type(esystem)              :: sys_nlpe
      !! non-linear poisson equation system
    type(esystem), allocatable :: sys_dd(:)
      !! electron/hole drift-diffusion equation system
    type(esystem)              :: sys_full
      !! full newton equation system
  contains
    procedure :: init        => device_init
    procedure :: destruct    => device_destruct
  end type

  type(device) :: dev
    !! global device

contains

  subroutine device_init(this, filename, T)
    !! initialize device using device file
    class(device), intent(out) :: this
    character(*),  intent(in)  :: filename
      !! device file name
    real,          intent(in)  :: T
      !! temperature in K

    integer                   :: ci, dir, idx_dir, ict, ivar, i0, i1, tmp
    character(:), allocatable :: vname
    type(input_file)          :: file

    ! load device parameters
    call file%init(filename)
    call this%par%init(file, T)

    ! init variables
    allocate (this%cdens(this%par%g%idx_dim,2))
    allocate (this%mob(this%par%g%idx_dim,2))
    call this%pot%init(this%par)
    if (this%par%smc%ii_pf .or. this%par%smc%ii_tun) then
      call this%efield_abs%init(this%par)
    end if
    do ci = this%par%ci0, this%par%ci1
      call this%dens(ci)%init(this%par, ci)
      call this%ion(ci)%init(this%par, ci)
      call this%genrec(ci)%init(this%par, ci)
      call this%eta(ci)%init(this%par, ci)
      call this%iref(ci)%init(this%par, ci)
      do idx_dir = 1, this%par%g%idx_dim
        call this%cdens(idx_dir,ci)%init(this%par, ci, idx_dir)
        call this%mob(idx_dir,ci)%init(this%par, ci, idx_dir)
      end do
    end do
    call this%rho%init(this%par)
    allocate (this%efield(this%par%g%dim))
    do idx_dir = 1, this%par%g%dim
      call this%efield(idx_dir)%init(this%par, idx_dir)
    end do
    allocate (this%volt(this%par%nct))
    allocate (this%curr(this%par%nct))
    do ict = 1, this%par%nct
      call this%volt(ict)%init("V_"//this%par%contacts(ict)%name)
      call this%curr(ict)%init("I_"//this%par%contacts(ict)%name)
    end do

    ! init Schottky boundary current variables (per contact, per carrier)
    ! must be initialized BEFORE continuity_init (passed as input)
    allocate (this%sbc(this%par%nct, 2))
    do ci = this%par%ci0, this%par%ci1
      do ict = 1, this%par%nct
        if (this%par%contacts(ict)%type == CT_SCHOTTKY) then
          call this%sbc(ict, ci)%init(this%par, ict, ci)
        end if
      end do
    end do

    ! init equations
    allocate (this%calc_cdens(this%par%g%idx_dim,2))
    if (this%par%smc%mob) allocate (this%calc_mob(this%par%g%idx_dim,2))
    call this%poiss%init(this%par, this%pot, this%rho, this%volt)
    call this%ramo%init(this%par, this%pot, this%rho, this%volt, this%poiss)
    call this%ramo_curr%init(this%par, this%ramo, this%cdens, this%volt, this%curr)
    do ci = this%par%ci0, this%par%ci1
      call this%contin(ci)%init(this%par, this%dens(ci), this%cdens(:,ci), this%genrec(ci), this%sbc(:,ci))
      call this%calc_eta_dens(ci)%init(this%par, this%dens(ci), this%eta(ci))
      call this%calc_iref(ci)%init(this%par, this%eta(ci), this%pot, this%iref(ci))
      call this%calc_eta_iref(ci)%init(this%par, this%pot, this%iref(ci), this%eta(ci))
      call this%calc_dens(ci)%init(this%par, this%eta(ci), this%dens(ci))
      if (this%par%smc%incomp_ion) then
        call this%calc_ion(ci)%init(this%par, this%eta(ci), this%ion(ci))
        call this%calc_genrec(ci)%init(this%par, this%genrec(ci), this%eta(ci), this%ion(ci), this%efield_abs)
        call this%ion_contin(ci)%init(this%par, this%ion(ci), this%genrec(ci))
      end if
      do idx_dir = 1, this%par%g%idx_dim
        if (this%par%smc%mob) call this%calc_mob(idx_dir,ci)%init(this%par, this%iref(ci), this%mob(idx_dir,ci))
        call this%calc_cdens(idx_dir,ci)%init(this%par, this%pot, this%dens(ci), this%cdens(idx_dir,ci), this%mob(idx_dir,ci))
      end do
    end do
    call this%calc_rho%init(this%par, this%dens, this%ion, this%rho)
    allocate (this%calc_efield(this%par%g%dim))
    do idx_dir = 1, this%par%g%dim
      call this%calc_efield(idx_dir)%init(this%par, this%efield(idx_dir), this%pot)
    end do

    ! Initialize calc_efield_abs only when field-assisted ionization is enabled
    if (this%par%smc%ii_pf .or. this%par%smc%ii_tun) then
      call this%calc_efield_abs%init(this%par, this%efield, this%efield_abs)
    end if

    ! init calc_schottky_bc equations (per Schottky contact, per carrier)
    allocate (this%calc_sbc(this%par%nct, 2))
    do ci = this%par%ci0, this%par%ci1
      do ict = 1, this%par%nct
        if (this%par%contacts(ict)%type == CT_SCHOTTKY) then
          call this%calc_sbc(ict, ci)%init(this%par, this%dens(ci), &
            & this%efield, this%sbc(ict, ci))
        end if
      end do
    end do

    ! init non-linear poisson equation system
    call this%sys_nlpe%init("non-linear poisson")
    call this%sys_nlpe%add_equation(this%poiss)
    call this%sys_nlpe%add_equation(this%calc_rho)
    do ci = this%par%ci0, this%par%ci1
      call this%sys_nlpe%add_equation(this%calc_eta_iref(ci))
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
    ! print general info about equation system
    print "(2A,I0,A)", this%sys_nlpe%name, ": ", this%sys_nlpe%n, " variables"
    tmp = 0
    do ivar = 1, this%sys_nlpe%g%imvar%n
      tmp = max(tmp, len(this%sys_nlpe%g%nodes%d(this%sys_nlpe%g%imvar%d(ivar))%p%v%name))
    end do
    allocate (character(len=tmp) :: vname)
    do ivar = 1, this%sys_nlpe%g%imvar%n
      i0 = this%sys_nlpe%i0(this%sys_nlpe%res2block(ivar)%d(1))
      i1 = this%sys_nlpe%i1(this%sys_nlpe%res2block(ivar)%d(size(this%sys_nlpe%res2block(ivar)%d)))
      vname(:) = this%sys_nlpe%g%nodes%d(this%sys_nlpe%g%imvar%d(ivar))%p%v%name
      print "(A,I0,A,I0)", this%sys_nlpe%name // ": " // vname // " goes from ", i0, " to ", i1
    end do
    deallocate (vname)
    ! output dependency graph to PDF file in run folder
    call this%sys_nlpe%g%output("nlpe")

    ! init drift-diffusion equation systems
    allocate (this%sys_dd(this%par%ci0:this%par%ci1))
    do ci = this%par%ci0, this%par%ci1
      call this%sys_dd(ci)%init(CR_NAME(ci)//" drift-diffusion")
      call this%sys_dd(ci)%add_equation(this%contin(ci))
      if (this%par%smc%ii_pf .or. this%par%smc%ii_tun) then
        call this%sys_dd(ci)%add_equation(this%calc_efield_abs)
      end if
      do idx_dir = 1, this%par%g%idx_dim
        call this%sys_dd(ci)%add_equation(this%calc_cdens(idx_dir,ci))
        if (this%par%smc%mob) call this%sys_dd(ci)%add_equation(this%calc_mob(idx_dir,ci))
      end do
      if (this%par%smc%incomp_ion) then
        call this%sys_dd(ci)%add_equation(this%calc_eta_iref(ci))
        call this%sys_dd(ci)%add_equation(this%ion_contin(ci))
        call this%sys_dd(ci)%add_equation(this%calc_genrec(ci))
      end if
      ! add E-field equations (consumed by calc_schottky_bc; lagged dependency)
      do dir = 1, this%par%g%dim
        call this%sys_dd(ci)%add_equation(this%calc_efield(dir))
      end do
      ! add per-Schottky-contact boundary equations
      do ict = 1, this%par%nct
        if (this%par%contacts(ict)%type == CT_SCHOTTKY) then
          call this%sys_dd(ci)%add_equation(this%calc_sbc(ict, ci))
        end if
      end do
      call this%sys_dd(ci)%provide(this%iref(ci), this%par%transport(IDX_VERTEX,0))
      call this%sys_dd(ci)%provide(this%pot, this%par%poisson(IDX_VERTEX,0))
      call this%sys_dd(ci)%init_final()
      ! print general info about equation system
      print "(2A,I0,A)", this%sys_dd(ci)%name, ": ", this%sys_dd(ci)%n, " variables"
      tmp = 0
      do ivar = 1, this%sys_dd(ci)%g%imvar%n
        tmp = max(tmp, len(this%sys_dd(ci)%g%nodes%d(this%sys_dd(ci)%g%imvar%d(ivar))%p%v%name))
      end do
      allocate (character(len=tmp) :: vname)
      do ivar = 1, this%sys_dd(ci)%g%imvar%n
        i0 = this%sys_dd(ci)%i0(this%sys_dd(ci)%res2block(ivar)%d(1))
        i1 = this%sys_dd(ci)%i1(this%sys_dd(ci)%res2block(ivar)%d(size(this%sys_dd(ci)%res2block(ivar)%d)))
        vname(:) = this%sys_dd(ci)%g%nodes%d(this%sys_dd(ci)%g%imvar%d(ivar))%p%v%name
        print "(A,I0,A,I0)", this%sys_dd(ci)%name // ": " // vname // " goes from ", i0, " to ", i1
      end do
      deallocate (vname)
      ! output dependency graph to PDF file in run folder
      call this%sys_dd(ci)%g%output(CR_NAME(ci)//"dd")
    end do

    ! init full-newton equation system
    call this%sys_full%init("full newton")
    call this%sys_full%add_equation(this%poiss)
    call this%sys_full%add_equation(this%calc_rho)
    if (this%par%smc%ii_pf .or. this%par%smc%ii_tun) then
      call this%sys_full%add_equation(this%calc_efield_abs)
    end if
    do ci = this%par%ci0, this%par%ci1
      call this%sys_full%add_equation(this%contin(ci))
      do idx_dir = 1, this%par%g%idx_dim
        call this%sys_full%add_equation(this%calc_cdens(idx_dir,ci))
        if (this%par%smc%mob) call this%sys_full%add_equation(this%calc_mob(idx_dir,ci))
      end do
      if (this%par%smc%incomp_ion) then
        call this%sys_full%add_equation(this%ion_contin(ci))
        call this%sys_full%add_equation(this%calc_genrec(ci))
      end if
      call this%sys_full%add_equation(this%calc_eta_dens(ci))
      call this%sys_full%add_equation(this%calc_iref(ci))
    end do
    call this%sys_full%add_equation(this%ramo_curr)
    ! add electric field calculation equations
    do dir = 1, this%par%g%dim
      call this%sys_full%add_equation(this%calc_efield(dir))
    end do
    ! add per-Schottky-contact boundary equations
    do ci = this%par%ci0, this%par%ci1
      do ict = 1, this%par%nct
        if (this%par%contacts(ict)%type == CT_SCHOTTKY) then
          call this%sys_full%add_equation(this%calc_sbc(ict, ci))
        end if
      end do
    end do
    do ict = 1, this%par%nct
      call this%sys_full%provide(this%volt(ict), input = .true.)
    end do
    call this%sys_full%init_final()
    ! print general info about equation system
    print "(2A,I0,A)", this%sys_full%name, ": ", this%sys_full%n, " variables"
    tmp = 0
    do ivar = 1, this%sys_full%g%imvar%n
      tmp = max(tmp, len(this%sys_full%g%nodes%d(this%sys_full%g%imvar%d(ivar))%p%v%name))
    end do
    allocate (character(len=tmp) :: vname)
    do ivar = 1, this%sys_full%g%imvar%n
      i0 = this%sys_full%i0(this%sys_full%res2block(ivar)%d(1))
      i1 = this%sys_full%i1(this%sys_full%res2block(ivar)%d(size(this%sys_full%res2block(ivar)%d)))
      vname(:) = this%sys_full%g%nodes%d(this%sys_full%g%imvar%d(ivar))%p%v%name
      print "(A,I0,A,I0)", this%sys_full%name // ": " // vname // " goes from ", i0, " to ", i1
    end do
    deallocate (vname)
    ! output dependency graph to PDF file in run folder
    call this%sys_full%g%output("full")
  end subroutine

  subroutine device_destruct(this)
    !! destruct device
    class(device), intent(inout) :: this

    call this%par%destruct()
  end subroutine

end module
