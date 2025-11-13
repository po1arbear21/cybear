module continuity_m

  use current_density_m, only: current_density
  use density_m,         only: density
  use device_params_m,   only: device_params
  use electric_field_m,  only: electric_field
  use error_m,           only: assert_failed, program_error
  use grid_m,            only: IDX_VERTEX, IDX_EDGE
  use imref_m,           only: imref
  use jacobian_m,        only: jacobian, jacobian_ptr
  use ionization_m,      only: generation_recombination
  use res_equation_m,    only: res_equation
  use semiconductor_m,   only: CR_CHARGE, CR_NAME, DOS_PARABOLIC, DIST_MAXWELL, CR_ELEC, CR_HOLE
  use contact_m,         only: CT_OHMIC, CT_GATE, CT_SCHOTTKY
  use normalization_m,   only: norm, denorm
  use stencil_m,         only: dirichlet_stencil, empty_stencil, near_neighb_stencil, stencil_ptr
  use vselector_m,       only: vselector
  use schottky_m,        only: schottky_injection_mb, schottky_velocity, schottky_injection_mb_bias, schottky_injection, schottky_tunnel_current

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
    type(electric_field), pointer :: efield(:) => null()
      !! electric field components (for Schottky barrier lowering)
    type(vselector)              :: n0b_inj
      !! dependency: Schottky injection density (field-dependent)
    type(vselector)              :: jtn_inj
      !! dependency: Schottky tunneling current density (field-dependent)

    real, allocatable :: b(:)
      !! right-hand side

    type(dirichlet_stencil)                :: st_dir
    type(near_neighb_stencil), allocatable :: st_nn(:)
    type(empty_stencil)                    :: st_em

    type(jacobian),     pointer     :: jaco_dens     => null()
    type(jacobian),     pointer     :: jaco_dens_t   => null()
    type(jacobian_ptr), allocatable :: jaco_cdens(:)
    type(jacobian),     pointer     :: jaco_genrec   => null()
    type(jacobian),     pointer     :: jaco_n0b      => null()
      !! jacobian for n0B dependency (field-dependent Schottky BC)
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

contains

  subroutine continuity_init(this, par, stat, dens, iref, cdens, genrec, efield, n0b, jtn_current)
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
    type(electric_field), optional, target, intent(in) :: efield(:)
      !! electric field components (optional, for Schottky barrier lowering)
    type(schottky_injection), optional, intent(in) :: n0b
      !! Schottky injection density (optional, for field-dependent BC)
    type(schottky_tunnel_current), optional, intent(in) :: jtn_current
      !! Tunneling current density (optional, for Schottky contacts with tunneling)

    integer              :: ci, i, ict, idx_dir, idens, idx_dim, igenrec, j, in0b, ijtn
    integer, allocatable :: idx(:), idx1(:), idx2(:), icdens(:)
    logical              :: status
    real                 :: surf, F, dF
    real                 :: A_ct, v_surf, n0b_val  ! Variables for Schottky BC
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

    ! store electric field reference if provided (for Schottky barrier lowering)
    if (present(efield)) then
      this%efield => efield
    end if

    ! init n0b selector if provided (for field-dependent Schottky BC)
    if (present(n0b)) then
      ! Only initialize for Schottky contact vertices
      ! We need all transport vertices, but only Schottky ones will use n0b
      call this%n0b_inj%init(n0b, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
    end if

    ! init jtn selector if provided (for tunneling current at Schottky contacts)
    if (present(jtn_current)) then
      ! Only initialize for Schottky contact vertices with tunneling enabled
      call this%jtn_inj%init(jtn_current, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
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

    ! n0b dependency
    if (present(n0b)) then
      in0b = this%depend(this%n0b_inj)
    end if

    ! jtn dependency
    if (present(jtn_current)) then
      ijtn = this%depend(this%jtn_inj)
    end if

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

    ! init n0b jacobian for field-dependent Schottky BC
    if (present(n0b)) then
      ! Check if any contacts are Schottky type
      if (any(par%contacts(1:par%nct)%type == CT_SCHOTTKY)) then
        ! n0B dependency only at Schottky contact vertices
        this%jaco_n0b => this%init_jaco_f(in0b, &
          & st = [this%st_em%get_ptr(), (merge(this%st_dir%get_ptr(), this%st_em%get_ptr(), &
          & par%contacts(ict)%type == CT_SCHOTTKY), ict = 1, par%nct)], &
          & const = .true., dtime = .false.)
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
          v_surf = schottky_velocity(par, ci, ict)
          A_ct = par%get_ct_surf(ict, idx1)

          ! Add Robin BC terms (use SET not ADD for contact vertices)
          call this%jaco_dens%set(idx1, idx1, A_ct * v_surf)

          if (present(n0b)) then
            ! Field-dependent n0B: set Jacobian entry
            call this%jaco_n0b%set(idx1, idx1, -A_ct * v_surf)
            ! RHS is set to 0 here, tunneling current added in eval()
            this%b(j) = 0.0
          else
            ! Constant n0B: set RHS
            call schottky_injection_mb(par, ci, ict, n0b_val)
            this%b(j) = A_ct * v_surf * n0b_val

            ! Debug output for Schottky BC (only print for first vertex at each contact)
            if (i == 1) then
              print "(A,A,A,I2,A)", "DEBUG_BC_INIT: Schottky BC at contact ", trim(par%contacts(ict)%name), &
                                    " (", ict, "):"
              print "(A,ES14.6,A)", "  A_ct           = ", A_ct, " (normalized)"
              print "(A,ES14.6,A)", "  A_ct (denorm)  = ", denorm(A_ct, "nm"), " nm"
              print "(A,ES14.6,A)", "  v_surf         = ", v_surf, " (normalized)"
              print "(A,ES14.6,A)", "  v_surf (denorm)= ", denorm(v_surf, "cm/s"), " cm/s"
              print "(A,ES14.6,A)", "  n0B (denorm)   = ", denorm(n0b_val, "cm^-3"), " cm^-3"
              print "(A,ES14.6)",   "  A*v*n0B (norm) = ", A_ct * v_surf * n0b_val
              print "(A,ES14.6,A)", "  RHS term       = ", denorm(this%b(j), "1/s"), " 1/s"
            end if
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

  subroutine add_tunneling_contribution(this, tmp)
    !! Add tunneling current contribution to residual for Schottky contacts
    class(continuity), intent(in)    :: this
    real,              intent(inout) :: tmp(:)

    integer :: ict, i, j
    integer :: idx(this%par%g%idx_dim)
    real :: A_ct, J_tn_contrib
    real, allocatable :: jtn_vec(:)

    ! Get the full tunneling current vector
    jtn_vec = this%jtn_inj%get()

    ! Loop over all contacts
    j = this%par%transport_vct(0)%n
    do ict = 1, this%par%nct
      ! Only process Schottky contacts with tunneling enabled
      if (this%par%contacts(ict)%type == CT_SCHOTTKY .and. &
          this%par%contacts(ict)%tunneling) then

        do i = 1, this%par%transport_vct(ict)%n
          j = j + 1
          idx = this%par%transport_vct(ict)%get_idx(i)

          ! Get contact surface area
          A_ct = this%par%get_ct_surf(ict, idx)

          ! Add tunneling current as boundary flux
          ! J_tn is current leaving semiconductor (positive for electrons leaving)
          ! Robin BC convention: J = A_ct * v_surf * (n - n0b) + A_ct * J_tn
          J_tn_contrib = A_ct * jtn_vec(j)
          tmp(j) = tmp(j) + J_tn_contrib
        end do
      else
        ! Skip non-Schottky or non-tunneling contacts
        j = j + this%par%transport_vct(ict)%n
      end if
    end do

    deallocate(jtn_vec)
  end subroutine

  subroutine continuity_eval(this)
    !! evaluate continuity equation
    class(continuity), intent(inout) :: this

    integer              :: idx_dir
    integer              :: idx1(this%par%g%idx_dim)
    real, allocatable    :: tmp(:)

    allocate (tmp(this%f%n))

    call this%jaco_dens%matr%mul_vec(this%dens%get(), tmp)
    do idx_dir = 1, this%par%g%idx_dim
      call this%jaco_cdens(idx_dir)%p%matr%mul_vec(this%cdens(idx_dir)%get(), tmp, fact_y = 1.0)
    end do
    if (this%par%smc%incomp_ion) call this%jaco_genrec%matr%mul_vec(this%genrec%get(), tmp, fact_y = 1.0)

    ! Add n0b contribution for field-dependent Schottky BC
    if (associated(this%jaco_n0b)) then
      call this%jaco_n0b%matr%mul_vec(this%n0b_inj%get(), tmp, fact_y = 1.0)
      if (this%par%nct > 0 .and. this%par%contacts(1)%type == CT_SCHOTTKY) then
        if (this%par%transport_vct(1)%n > 0) then
          idx1 = this%par%transport_vct(1)%get_idx(1)
          print "(A,2ES14.6)", "DEBUG_SCHOTTKY: n, n0B(E) = ", &
            denorm(this%dens%get(idx1), "cm^-3"), denorm(this%n0b_inj%get(idx1), "cm^-3")
        end if
      end if
    end if

    ! Add tunneling current contribution for Schottky contacts with tunneling enabled
    if (this%jtn_inj%n > 0) then
      call add_tunneling_contribution(this, tmp)
    end if

    call this%f%set(tmp - this%b)
  end subroutine

end module
