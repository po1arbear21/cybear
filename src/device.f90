module device_m

  use charge_density_m,  only: charge_density, calc_charge_density
  use continuity_m,      only: continuity
  use current_m,         only: current
  use current_density_m, only: current_density, calc_current_density
  use density_m,         only: density
  use device_params_m,   only: device_params
  use electric_field_m,  only: electric_field, calc_efield
  use esystem_m,         only: esystem
  use grid_m,            only: IDX_VERTEX
  use imref_m,           only: imref, calc_imref, calc_density
  use input_m,           only: input_file
  use ionization_m,      only: ionization, calc_ionization, ion_continuity, generation_recombination, calc_generation_recombination
  use beam_generation_m,  only: beam_generation, beam_position, calc_beam_generation
  use srh_recombination_m, only: srh_recombination, calc_srh_recombination, &
    &                             surface_srh_recombination, calc_surface_srh
  use mobility_m,        only: mobility, calc_mobility
  use poisson_m,         only: poisson
  use potential_m,       only: potential
  use ramo_shockley_m,   only: ramo_shockley, ramo_shockley_current
  use semiconductor_m,   only: CR_NAME, CR_ELEC, CR_HOLE
  use voltage_m,         only: voltage

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
    type(beam_generation)              :: bgen(2)
      !! external beam generation rate (STEM-EBIC)
    type(beam_position)                :: beam_pos
      !! beam position for STEM-EBIC sweep (Y_BEAM)
    type(current_density), allocatable :: cdens(:,:)
      !! electron/hole current density (direction, carrier index)
    type(imref)                        :: iref(2)
      !! electron/hole quasi-fermi potential (carrier index)
    type(mobility),        allocatable :: mob(:,:)
      !! electron/hole mobility (direction, carrier index)
    type(charge_density)               :: rho
      !! charge density
    type(electric_field),  allocatable :: efield(:)
      !! electric field components (direction)
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
    type(calc_beam_generation)               :: calc_bgen(2)
      !! calculate external beam generation rate (STEM-EBIC)
    type(srh_recombination)                  :: srh
      !! SRH recombination rate (single variable for both carriers)
    type(calc_srh_recombination)             :: calc_srh
      !! calculate SRH recombination rate
    type(surface_srh_recombination)          :: surf_srh
      !! surface SRH recombination rate (single variable for both carriers)
    type(calc_surface_srh)                   :: calc_surf_srh
      !! calculate surface SRH recombination rate
    type(calc_mobility),        allocatable :: calc_mob(:,:)
      !! calculate electron/hole mobility using Caughey-Thomas model (direction, carrier index)
    type(calc_charge_density)               :: calc_rho
      !! calculate charge density from electron/hole density and (ionized) doping concentrations
    type(calc_current_density), allocatable :: calc_cdens(:,:)
      !! calculate electron/hole current density by drift-diffusion model (direction, carrier index)
    type(calc_efield),          allocatable :: calc_efield(:)
      !! calculate electric field from potential gradient (direction)

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

    integer          :: ci, idx_dir, ict, dir
    type(input_file) :: file

    ! load device parameters
    call file%init(filename)
    call this%par%init(file, T)

    ! init variables
    allocate (this%cdens(this%par%g%idx_dim,2))
    allocate (this%mob(this%par%g%idx_dim,2))
    call this%pot%init(this%par)
    do ci = this%par%ci0, this%par%ci1
      call this%dens(ci)%init(this%par, ci)
      call this%ion( ci)%init(this%par, ci)
      call this%genrec( ci)%init(this%par, ci)
      if (this%par%has_beam_gen) call this%bgen(ci)%init(this%par, ci)
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
    ! init beam position for STEM-EBIC sweep
    if (this%par%has_beam_gen) then
      call this%beam_pos%init("Y_BEAM")
    end if

    ! init equations
    allocate (this%calc_cdens(this%par%g%idx_dim,2))
    if (this%par%smc%mob) allocate (this%calc_mob(this%par%g%idx_dim,2))
    call this%poiss%init(this%par, this%pot, this%rho, this%volt)
    call this%ramo%init(this%par, this%pot, this%rho, this%volt, this%poiss)
    call this%ramo_curr%init(this%par, this%ramo, this%cdens, this%volt, this%curr)
    ! Initialize SRH recombination if enabled (single instance for both carriers)
    if (this%par%smc%srh) then
      call this%srh%init(this%par)
      call this%calc_srh%init(this%par, this%dens(CR_ELEC), this%dens(CR_HOLE), this%srh)
    end if

    ! Initialize surface SRH recombination if enabled
    if (this%par%smc%surf_recom) then
      call this%surf_srh%init(this%par)
      call this%calc_surf_srh%init(this%par, this%dens(CR_ELEC), this%dens(CR_HOLE), this%surf_srh)
    end if

    do ci = this%par%ci0, this%par%ci1
      ! Initialize continuity with all optional physics (beam, SRH, surface SRH)
      if (this%par%has_beam_gen .and. this%par%smc%srh .and. this%par%smc%surf_recom) then
        call this%contin(ci)%init(this%par, this%dens(ci), this%cdens(:,ci), this%genrec(ci), &
          & bgen=this%bgen(ci), srh=this%srh, &
          & surf_srh=this%surf_srh, surf_idx=this%calc_surf_srh%surf_idx, A_surf=this%calc_surf_srh%A_surf)
        call this%calc_bgen(ci)%init(this%par, this%bgen(ci), this%beam_pos)
      elseif (this%par%has_beam_gen .and. this%par%smc%srh) then
        call this%contin(ci)%init(this%par, this%dens(ci), this%cdens(:,ci), this%genrec(ci), bgen=this%bgen(ci), srh=this%srh)
        call this%calc_bgen(ci)%init(this%par, this%bgen(ci), this%beam_pos)
      elseif (this%par%has_beam_gen .and. this%par%smc%surf_recom) then
        call this%contin(ci)%init(this%par, this%dens(ci), this%cdens(:,ci), this%genrec(ci), &
          & bgen=this%bgen(ci), &
          & surf_srh=this%surf_srh, surf_idx=this%calc_surf_srh%surf_idx, A_surf=this%calc_surf_srh%A_surf)
        call this%calc_bgen(ci)%init(this%par, this%bgen(ci), this%beam_pos)
      elseif (this%par%has_beam_gen) then
        call this%contin(ci)%init(this%par, this%dens(ci), this%cdens(:,ci), this%genrec(ci), bgen=this%bgen(ci))
        call this%calc_bgen(ci)%init(this%par, this%bgen(ci), this%beam_pos)
      elseif (this%par%smc%srh .and. this%par%smc%surf_recom) then
        call this%contin(ci)%init(this%par, this%dens(ci), this%cdens(:,ci), this%genrec(ci), &
          & srh=this%srh, &
          & surf_srh=this%surf_srh, surf_idx=this%calc_surf_srh%surf_idx, A_surf=this%calc_surf_srh%A_surf)
      elseif (this%par%smc%srh) then
        call this%contin(ci)%init(this%par, this%dens(ci), this%cdens(:,ci), this%genrec(ci), srh=this%srh)
      elseif (this%par%smc%surf_recom) then
        call this%contin(ci)%init(this%par, this%dens(ci), this%cdens(:,ci), this%genrec(ci), &
          & surf_srh=this%surf_srh, surf_idx=this%calc_surf_srh%surf_idx, A_surf=this%calc_surf_srh%A_surf)
      else
        call this%contin(ci)%init(this%par, this%dens(ci), this%cdens(:,ci), this%genrec(ci))
      end if
      call this%calc_iref(ci)%init(this%par, this%pot, this%dens(ci), this%iref(ci))
      call this%calc_dens(ci)%init(this%par, this%pot, this%dens(ci), this%iref(ci))
      if (this%par%smc%incomp_ion) then
        call this%calc_ion(ci)%init(this%par, this%pot, this%ion(ci),  this%iref(ci))
        call this%calc_genrec(ci)%init(this%par, this%par%smc%ii_tau(ci), this%genrec(ci), this%pot, this%iref(ci), this%ion(ci))
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
      call this%calc_efield(idx_dir)%init(this%par, this%pot, this%efield(idx_dir))
    end do

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
    print "(2A,I0,A)", this%sys_nlpe%name, ": ", this%sys_nlpe%n, " variables"
    call this%sys_nlpe%g%output("nlpe")

    ! init drift-diffusion equation systems
    allocate (this%sys_dd(this%par%ci0:this%par%ci1))
    do ci = this%par%ci0, this%par%ci1
      call this%sys_dd(ci)%init(CR_NAME(ci)//" drift-diffusion")
      call this%sys_dd(ci)%add_equation(this%contin(ci))
      do idx_dir = 1, this%par%g%idx_dim
        call this%sys_dd(ci)%add_equation(this%calc_cdens(idx_dir,ci))
        if (this%par%smc%mob) call this%sys_dd(ci)%add_equation(this%calc_mob(idx_dir,ci))
      end do
      if (this%par%smc%incomp_ion) then
        call this%sys_dd(ci)%add_equation(this%ion_contin(ci))
        call this%sys_dd(ci)%add_equation(this%calc_genrec(ci))
      end if
      if (this%par%has_beam_gen) then
        call this%sys_dd(ci)%add_equation(this%calc_bgen(ci))
      end if
      if (this%par%smc%srh) then
        call this%sys_dd(ci)%add_equation(this%calc_srh)
      end if
      if (this%par%smc%surf_recom) then
        call this%sys_dd(ci)%add_equation(this%calc_surf_srh)
      end if
      ! Provide the OTHER carrier's density as external input (for SRH or surface SRH)
      if (this%par%smc%srh .or. this%par%smc%surf_recom) then
        if (ci == CR_ELEC) then
          call this%sys_dd(ci)%provide(this%dens(CR_HOLE), this%par%transport(IDX_VERTEX,0))
        else
          call this%sys_dd(ci)%provide(this%dens(CR_ELEC), this%par%transport(IDX_VERTEX,0))
        end if
      end if
      call this%sys_dd(ci)%provide(this%iref(ci), this%par%transport(IDX_VERTEX,0))
      call this%sys_dd(ci)%provide(this%pot, this%par%transport(IDX_VERTEX,0))
      call this%sys_dd(ci)%init_final()
      print "(2A,I0,A)", this%sys_dd(ci)%name, ": ", this%sys_dd(ci)%n, " variables"
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
        if (this%par%smc%mob) call this%sys_full%add_equation(this%calc_mob(idx_dir,ci))
      end do
      if (this%par%smc%incomp_ion) then
        call this%sys_full%add_equation(this%ion_contin(ci))
        call this%sys_full%add_equation(this%calc_genrec(ci))
      end if
      if (this%par%has_beam_gen) then
        call this%sys_full%add_equation(this%calc_bgen(ci))
      end if
      call this%sys_full%add_equation(this%calc_iref(ci))
    end do
    ! SRH recombination (single equation for both carriers - add OUTSIDE carrier loop)
    if (this%par%smc%srh) then
      call this%sys_full%add_equation(this%calc_srh)
    end if
    ! Surface SRH recombination (single equation for both carriers - add OUTSIDE carrier loop)
    if (this%par%smc%surf_recom) then
      call this%sys_full%add_equation(this%calc_surf_srh)
    end if
    call this%sys_full%add_equation(this%ramo_curr)
    ! add electric field calculation equations
    do dir = 1, this%par%g%dim
      call this%sys_full%add_equation(this%calc_efield(dir))
    end do
    do ict = 1, this%par%nct
      call this%sys_full%provide(this%volt(ict), input = .true.)
    end do
    ! register beam position as input for STEM-EBIC sweep
    if (this%par%has_beam_gen) then
      call this%sys_full%provide(this%beam_pos, input = .true.)
    end if
    call this%sys_full%init_final()
    print "(2A,I0,A)", this%sys_full%name, ": ", this%sys_full%n, " variables"
    call this%sys_full%g%output("full")
  end subroutine

  subroutine device_destruct(this)
    !! destruct device
    class(device), intent(inout) :: this

    call this%par%destruct()
  end subroutine

end module
