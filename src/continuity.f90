module continuity_m

  use current_density_m, only: current_density
  use density_m,         only: density
  use device_params_m,   only: device_params
  use error_m,           only: assert_failed, program_error
  use grid_m,            only: IDX_VERTEX, IDX_EDGE
  use imref_m,           only: imref
  use jacobian_m,        only: jacobian, jacobian_ptr
  use ionization_m,      only: generation_recombination
  use potential_m,       only: potential
  use res_equation_m,    only: res_equation
  use semiconductor_m,   only: CR_CHARGE, CR_NAME, DOS_PARABOLIC, DIST_MAXWELL, CR_ELEC, CR_HOLE
  use contact_m,         only: CT_OHMIC, CT_GATE, CT_SCHOTTKY
  use normalization_m,   only: norm, denorm
  use stencil_m,         only: dirichlet_stencil, empty_stencil, near_neighb_stencil, stencil_ptr
  use vselector_m,       only: vselector
  use schottky_m,        only: schottky_injection_mb, schottky_velocity, schottky_injection_mb_bias

  implicit none

  private
  public continuity

  type, extends(res_equation) :: continuity
    !! dn/dt + div j = 0

    integer :: ci
      !! carrier index

    type(device_params), pointer :: par => null()
      !! pointer to device parameters
    type(vselector)              :: dens
      !! main variable: density
    type(vselector)              :: iref
      !! main variable: quasi-fermi potential
    type(vselector), allocatable :: cdens(:)
      !! dependencies: current densities in 2 directions
    type(vselector)              :: genrec
      !! generation - recombination
    type(vselector)              :: pot
      !! dependency: potential (for bias-dependent Schottky BC)

    real, allocatable :: b(:)
      !! right-hand side

    type(dirichlet_stencil)                :: st_dir
    type(near_neighb_stencil), allocatable :: st_nn(:)
    type(empty_stencil)                    :: st_em

    type(jacobian),     pointer     :: jaco_dens     => null()
    type(jacobian),     pointer     :: jaco_dens_t   => null()
    type(jacobian_ptr), allocatable :: jaco_cdens(:)
    type(jacobian),     pointer     :: jaco_genrec   => null()
    type(jacobian),     pointer     :: jaco_pot      => null()
      !! jacobian for potential dependency (Schottky BC)
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

contains

  subroutine continuity_init(this, par, stat, dens, iref, cdens, genrec, pot)
    !! initialize continuity equation
    class(continuity),              intent(out) :: this
    type(device_params), target,    intent(in)  :: par
      !! device parameters
    logical,                        intent(in)  :: stat
      !! stationary? if true, use imref as main variable
    type(density),                  intent(in)  :: dens
      !! electron/hole density
    type(imref),                    intent(in)  :: iref
      !! electron/hole quasi-fermi potential
    type(current_density),          intent(in)  :: cdens(:)
      !! electron/hole current density
    type(generation_recombination), intent(in)  :: genrec
      !! generation-recombination rate
    type(potential), optional,      intent(in)  :: pot
      !! potential (optional, for bias-dependent Schottky BC)

    integer              :: ci, i, ict, idx_dir, idens, idx_dim, igenrec, j, ipot
    integer, allocatable :: idx(:), idx1(:), idx2(:), icdens(:)
    logical              :: status
    real                 :: surf, F, dF
    real                 :: A_ct, v_surf, n0b  ! Variables for Schottky BC
    type(stencil_ptr), allocatable :: st_dens_ct(:), st_dens_t_ct(:), st_cdens_ct(:)

    print "(A)", "continuity_init"

    ci = dens%ci

    ! init base
    if (stat) then
      call this%equation_init(CR_NAME(dens%ci)//"continuity_stat")
    else
      call this%equation_init(CR_NAME(dens%ci)//"continuity")
    end if
    this%par => par
    this%ci  = ci

    idx_dim = par%g%idx_dim
    allocate (idx(idx_dim), idx1(idx_dim), idx2(idx_dim), icdens(idx_dim))
    allocate (this%cdens(idx_dim), this%st_nn(idx_dim), this%jaco_cdens(idx_dim))

    ! init variable selectors
    call this%dens%init(dens, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
    do idx_dir = 1, idx_dim
      call this%cdens(idx_dir)%init(cdens(idx_dir), par%transport(IDX_EDGE,idx_dir))
    end do
    if (par%smc%incomp_ion) call this%genrec%init(genrec, par%ionvert(ci))

    ! init potential selector if provided (for bias-dependent Schottky BC)
    if (present(pot)) then
      call this%pot%init(pot, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
    end if

    ! init residuals using this%dens or this%iref as main variable
    if (stat) then
      call this%iref%init(iref, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
      call this%init_f(this%iref)
    else
      call this%init_f(this%dens)
    end if

    ! init stencils
    call this%st_dir%init(par%g)
    do idx_dir = 1, idx_dim
      call this%st_nn(idx_dir)%init(par%g, IDX_VERTEX, 0, IDX_EDGE, idx_dir)
    end do
    call this%st_em%init()

    ! setup contact-specific stencils
    allocate(st_cdens_ct(par%nct), st_dens_ct(par%nct), st_dens_t_ct(par%nct))

    ! For time-dependent jaco_dens_t: differentiate by contact type
    do ict = 1, par%nct
      st_dens_ct(ict) = this%st_dir%get_ptr()
      if (par%contacts(ict)%type == CT_SCHOTTKY) then
        st_dens_t_ct(ict) = this%st_dir%get_ptr()  ! Time and genrec terms for Schottky
      else
        st_dens_t_ct(ict) = this%st_em%get_ptr()   ! No time/genrec for Ohmic/Gate
      end if
    end do

    ! dependencies
    idens = this%depend(this%dens)
    do idx_dir = 1, idx_dim
      icdens(idx_dir) = this%depend(this%cdens(idx_dir))
    end do
    if (par%smc%incomp_ion) igenrec = this%depend(this%genrec)
    if (present(pot)) ipot = this%depend(this%pot)

    ! init jacobians
    this%jaco_dens   => this%init_jaco_f(idens, &
      & st = [this%st_em%get_ptr(), (st_dens_ct(ict), ict = 1, par%nct)], &
      & const = .true., dtime = .false.)
    if (.not. stat) then
      this%jaco_dens_t => this%init_jaco_f(idens, &
        & st = [this%st_dir%get_ptr(), (st_dens_t_ct(ict), ict = 1, par%nct)], &
        & const = .true., dtime = .true. )
    end if
    do idx_dir = 1, idx_dim
      ! setup contact-specific stencils for this direction
      do ict = 1, par%nct
        if (par%contacts(ict)%type == CT_SCHOTTKY) then
          st_cdens_ct(ict) = this%st_nn(idx_dir)%get_ptr()  ! Edge contributions for Schottky
        else
          st_cdens_ct(ict) = this%st_em%get_ptr()           ! No edge contributions for Ohmic/Gate
        end if
      end do

      this%jaco_cdens(idx_dir)%p => this%init_jaco_f(icdens(idx_dir), &
        & st = [this%st_nn(idx_dir)%get_ptr(), (st_cdens_ct(ict), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
    end do
    if (par%smc%incomp_ion) then
      this%jaco_genrec => this%init_jaco_f(igenrec, &
        & st = [this%st_dir%get_ptr(), (st_dens_t_ct(ict), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
    end if

    ! init potential jacobian for bias-dependent Schottky BC
    if (present(pot)) then
      ! Check if any contacts are Schottky type
      if (any(par%contacts(1:par%nct)%type == CT_SCHOTTKY)) then
        ! For now, use empty stencil for interior and dir stencil for Schottky contacts
        ! This will be refined when we implement the actual bias dependency
        this%jaco_pot => this%init_jaco_f(ipot, &
          & st = [this%st_em%get_ptr(), (merge(this%st_dir%get_ptr(), this%st_em%get_ptr(), &
          & par%contacts(ict)%type == CT_SCHOTTKY), ict = 1, par%nct)], &
          & const = .false., dtime = .false.)
      end if
    end if

    ! set current density jacobian entries
    do idx_dir = 1, idx_dim
      do i = 1, par%transport(IDX_EDGE,idx_dir)%n
        ! get indices of edge, first and second vertex
        idx = par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        ! get edge surface
        surf = par%tr_surf(idx_dir)%get(idx)

        ! Debug edge surfaces near contacts
        if (par%ict%get(idx1) > 0 .or. par%ict%get(idx2) > 0) then
          if (par%ict%get(idx1) == 1 .or. par%ict%get(idx2) == 1) then  ! Near Schottky
            print *, "    Edge ", idx, " surface = ", surf, " (normalized)"
          end if
        end if

        ! set values if uncontacted or Schottky
        ! Include edge contribution if vertex is interior (0) or Schottky contact
        if (par%ict%get(idx1) == 0) then
          call this%jaco_cdens(idx_dir)%p%set(idx1, idx,  surf)
        else if (par%contacts(par%ict%get(idx1))%type == CT_SCHOTTKY) then
          call this%jaco_cdens(idx_dir)%p%set(idx1, idx,  surf)
        end if

        if (par%ict%get(idx2) == 0) then
          call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -surf)
        else if (par%contacts(par%ict%get(idx2))%type == CT_SCHOTTKY) then
          call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -surf)
        end if
      end do
    end do

    ! density time derivative factor
    if (.not. stat) then
      do i = 1, par%transport_vct(0)%n
        idx1 = par%transport_vct(0)%get_idx(i)
        call this%jaco_dens_t%set(idx1, idx1, par%tr_vol%get(idx1))
      end do
    end if

    ! generation-recombination
    if (par%smc%incomp_ion) then
      do i = 1, par%ionvert(ci)%n
        idx1 = par%ionvert(ci)%get_idx(i)
        call this%jaco_genrec%set(idx1, idx1, - par%tr_vol%get(idx1))
      end do
    end if

    ! boundary conditions: Dirichlet for Ohmic/Gate, Robin for Schottky
    allocate (this%b(this%f%n), source = 0.0)
    j = par%transport_vct(0)%n
    do ict = 1, par%nct
      do i = 1, par%transport_vct(ict)%n
        j = j + 1
        idx1 = par%transport_vct(ict)%get_idx(i)

        if (par%contacts(ict)%type == CT_SCHOTTKY) then
          ! Robin BC for Schottky: Add boundary flux terms
          ! Get Schottky BC parameters
          call schottky_injection_mb(par, ci, ict, n0b)
          v_surf = schottky_velocity(par, ci, ict)
          A_ct = par%get_ct_surf(ict, idx1)

          ! Debug output for first vertex of each contact
          if (i == 1) then
            print *, "DEBUG Robin BC at contact ", ict, " for carrier ", ci
            print *, "  n0b (normalized) = ", n0b
            print *, "  n0b (physical) = ", denorm(n0b, "1/cm^3"), " 1/cm^3"
            print *, "  v_surf (normalized) = ", v_surf
            print *, "  v_surf (physical) = ", denorm(v_surf, "cm/s"), " cm/s"
            print *, "  A_ct = ", A_ct
            print *, "  Jacobian diagonal += ", A_ct * v_surf
            print *, "  RHS += ", A_ct * v_surf * n0b
          end if

          ! Add Robin BC terms (use SET not ADD for contact vertices)
          call this%jaco_dens%set(idx1, idx1, A_ct * v_surf)
          this%b(j) = A_ct * v_surf * n0b

          ! Add potential derivative for bias-dependent BC
          if (associated(this%jaco_pot)) then
            ! Placeholder: zero derivative for now (will be updated in eval)
            call this%jaco_pot%set(idx1, idx1, 0.0)
          end if
        else
          ! Dirichlet BC for Ohmic/Gate contacts
          call this%jaco_dens%set(idx1, idx1, 1.0)
          if ((par%smc%dos == DOS_PARABOLIC) .and. (par%smc%dist == DIST_MAXWELL)) then
            this%b(j) = sqrt(par%smc%edos(1) * par%smc%edos(2)) * exp(- CR_CHARGE(ci) * par%contacts(ict)%phims - 0.5 * par%smc%band_gap)
          else
            call par%smc%get_dist(- CR_CHARGE(ci) * (par%contacts(ict)%phims - par%smc%band_edge(ci)), 0, F, dF)
            this%b(j) = par%smc%edos(ci) * F
          end if
        end if
      end do
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine continuity_eval(this)
    !! evaluate continuity equation
    class(continuity), intent(inout) :: this

    integer              :: idx_dir, i, j, ict
    integer, allocatable :: idx1(:)
    real, allocatable    :: tmp(:)
    real                 :: A_ct, v_surf, n0b, dn0b_dE, phi0
    real                 :: E_field  ! Electric field for barrier lowering
    real                 :: phi0_arr(1)  ! Temporary array for vselector get

    allocate (tmp(this%f%n))
    allocate (idx1(this%par%g%idx_dim))

    ! Update boundary conditions for Schottky contacts (bias-dependent)
    if (associated(this%jaco_pot)) then
      j = this%par%transport_vct(0)%n
      do ict = 1, this%par%nct
        if (this%par%contacts(ict)%type == CT_SCHOTTKY) then
          ! Get surface velocity once per contact
          v_surf = schottky_velocity(this%par, this%ci, ict)
        end if

        do i = 1, this%par%transport_vct(ict)%n
          j = j + 1
          idx1 = this%par%transport_vct(ict)%get_idx(i)

          if (this%par%contacts(ict)%type == CT_SCHOTTKY) then
            ! Get contact area for this vertex
            A_ct = this%par%get_ct_surf(ict, idx1)

            ! TODO: Calculate electric field at contact and dE/dphi
            ! For now, use placeholder zero field
            E_field = 0.0

            ! Get bias-dependent injection density and derivative w.r.t. E field
            call schottky_injection_mb_bias(this%par, this%ci, ict, E_field, n0b, dn0b_dE)

            ! Get potential at contact vertex for linearization
            phi0_arr = this%pot%get(idx1)
            phi0 = phi0_arr(1)

            ! Update RHS with linearization: b = Av(n0b - dn0b/dphi * phi0)
            ! Note: dn0b/dphi = dn0b/dE * dE/dphi, but dE/dphi = 0 for now (placeholder)
            this%b(j) = A_ct * v_surf * (n0b - dn0b_dE * 0.0)  ! No phi coupling yet

            ! Update potential Jacobian entry
            ! Will be -A*v*dn0b/dphi when dE/dphi is implemented
            call this%jaco_pot%set(idx1, idx1, 0.0)  ! Placeholder until dE/dphi available
          end if
        end do
      end do

      ! Materialize the Jacobian matrix after all entries are set
      call this%jaco_pot%set_matr(const = .false., nonconst = .true.)
    end if

    call this%jaco_dens%matr%mul_vec(this%dens%get(), tmp)
    do idx_dir = 1, this%par%g%idx_dim
      call this%jaco_cdens(idx_dir)%p%matr%mul_vec(this%cdens(idx_dir)%get(), tmp, fact_y = 1.0)
    end do
    if (this%par%smc%incomp_ion) call this%jaco_genrec%matr%mul_vec(this%genrec%get(), tmp, fact_y = 1.0)

    ! Add potential contribution for bias-dependent Schottky BC
    if (associated(this%jaco_pot)) then
      call this%jaco_pot%matr%mul_vec(this%pot%get(), tmp, fact_y = 1.0)
    end if

    call this%f%set(tmp - this%b)
  end subroutine

end module
