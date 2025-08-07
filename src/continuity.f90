module continuity_m

  use current_density_m, only: current_density
  use density_m,         only: density
  use device_params_m,   only: device_params
  use error_m,           only: assert_failed, program_error
  use grid_m,            only: IDX_VERTEX, IDX_EDGE
  use imref_m,           only: imref
  use jacobian_m,        only: jacobian, jacobian_ptr
  use ionization_m,      only: generation_recombination
  use normalization_m,   only: denorm
  use potential_m,       only: potential
  use res_equation_m,    only: res_equation
  use semiconductor_m,   only: CR_CHARGE, CR_NAME, DOS_PARABOLIC, DIST_MAXWELL, CR_ELEC, CR_HOLE
  use contact_m,         only: CT_SCHOTTKY, CT_OHMIC, CT_GATE
  use stencil_m,         only: dirichlet_stencil, empty_stencil, near_neighb_stencil, stencil_ptr
  use vselector_m,       only: vselector

  implicit none

  private
  public continuity

  type, extends(res_equation) :: continuity
    !! dn/dt + div j = 0

    integer :: ci
      !! carrier index

    type(device_params), pointer :: par => null()
      !! pointer to device parameters
    type(potential),     pointer :: pot => null()
      !! pointer to electrostatic potential (for bias-dependent BC)
    type(vselector)              :: dens
      !! main variable: density
    type(vselector)              :: iref
      !! main variable: quasi-fermi potential
    type(vselector), allocatable :: cdens(:)
      !! dependencies: current densities in 2 directions
    type(vselector)              :: genrec
      !! generation - recombination

    real, allocatable :: b(:)
      !! right-hand side
    real, allocatable :: b_edge(:,:)
      !! edge contributions to RHS for Robin BC (edge_index, direction)

    type(dirichlet_stencil)                :: st_dir
    type(near_neighb_stencil), allocatable :: st_nn(:)
    type(empty_stencil)                    :: st_em

    type(jacobian),     pointer     :: jaco_dens     => null()
    type(jacobian),     pointer     :: jaco_dens_t   => null()
    type(jacobian_ptr), allocatable :: jaco_cdens(:)
    type(jacobian),     pointer     :: jaco_genrec   => null()
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

contains

  subroutine continuity_init(this, par, stat, pot, dens, iref, cdens, genrec)
    !! initialize continuity equation
    class(continuity),              intent(out) :: this
    type(device_params), target,    intent(in)  :: par
      !! device parameters
    logical,                        intent(in)  :: stat
      !! stationary? if true, use imref as main variable
    type(potential), target,        intent(in)  :: pot
      !! electrostatic potential (for bias-dependent BC)
    type(density),                  intent(in)  :: dens
      !! electron/hole density
    type(imref),                    intent(in)  :: iref
      !! electron/hole quasi-fermi potential
    type(current_density),          intent(in)  :: cdens(:)
      !! electron/hole current density
    type(generation_recombination), intent(in)  :: genrec
      !! generation-recombination rate

    integer              :: ci, i, ict, idx_dir, idens, idx_dim, igenrec, j
    integer, allocatable :: idx(:), idx1(:), idx2(:), icdens(:)
    integer              :: ict1, ict2
    logical              :: status
    real                 :: surf, F, dF, dx, S, D, n0
    type(stencil_ptr), allocatable :: st_cdens(:)

    print "(A)", "continuity_init"

    ci = dens%ci

    ! init base
    if (stat) then
      call this%equation_init(CR_NAME(dens%ci)//"continuity_stat")
    else
      call this%equation_init(CR_NAME(dens%ci)//"continuity")
    end if
    this%par => par
    this%pot => pot
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

    ! dependencies
    idens = this%depend(this%dens)
    do idx_dir = 1, idx_dim
      icdens(idx_dir) = this%depend(this%cdens(idx_dir))
    end do
    if (par%smc%incomp_ion) igenrec = this%depend(this%genrec)

    ! init jacobians
    this%jaco_dens   => this%init_jaco_f(idens, &
      & st = [this%st_em%get_ptr(), (this%st_dir%get_ptr(), ict = 1, par%nct)], &
      & const = .true., dtime = .false.)
    if (.not. stat) then
      this%jaco_dens_t => this%init_jaco_f(idens, &
        & st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .true. )
    end if
    ! Build stencil arrays based on contact types
    do idx_dir = 1, idx_dim
      ! Create stencil array: index 0 for uncontacted, 1:nct for contacts
      allocate(st_cdens(0:par%nct))

      ! Uncontacted vertices always use near_neighb stencil
      st_cdens(0) = this%st_nn(idx_dir)%get_ptr()

      ! For each contact, choose stencil based on contact type
      do ict = 1, par%nct
        select case (par%contacts(ict)%type)
        case (CT_SCHOTTKY)
          ! Schottky contacts need edge contributions for Robin BC
          st_cdens(ict) = this%st_nn(idx_dir)%get_ptr()
        case (CT_OHMIC, CT_GATE)
          ! Ohmic/Gate contacts use empty stencil (no edge contributions)
          st_cdens(ict) = this%st_em%get_ptr()
        case default
          ! Default to empty for unknown contact types
          st_cdens(ict) = this%st_em%get_ptr()
        end select
      end do

      this%jaco_cdens(idx_dir)%p => this%init_jaco_f(icdens(idx_dir), &
        & st = st_cdens, const = .true., dtime = .false.)

      deallocate(st_cdens)
    end do
    if (par%smc%incomp_ion) then
      this%jaco_genrec => this%init_jaco_f(igenrec, &
        & st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
    end if

    ! allocate edge RHS storage for Robin BC
    allocate (this%b_edge(maxval([(par%transport(IDX_EDGE,idx_dir)%n, idx_dir = 1, idx_dim)]), idx_dim), &
              source = 0.0)

    ! set current density jacobian entries
    do idx_dir = 1, idx_dim
      do i = 1, par%transport(IDX_EDGE,idx_dir)%n
        ! get indices of edge, first and second vertex
        idx = par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        ! get edge surface and length
        surf = par%tr_surf(idx_dir)%get(idx)
        dx = par%g%get_len(idx, idx_dir)

        ! check contact types
        ict1 = par%ict%get(idx1)
        ict2 = par%ict%get(idx2)

        ! handle different vertex combinations
        if (ict1 == 0 .and. ict2 == 0) then
          ! both uncontacted: standard flux continuity
          call this%jaco_cdens(idx_dir)%p%set(idx1, idx,  surf)
          call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -surf)

        else if (ict1 > 0 .and. ict2 == 0) then
          if (par%contacts(ict1)%type == CT_SCHOTTKY) then
            ! Schottky contact (vertex 1) to interior (vertex 2): Robin BC
            n0 = calculate_equilibrium_density(par, ci, ict1)
            S = calculate_thermionic_velocity(par, ci, ict1)
            D = par%mob0(IDX_EDGE, idx_dir, ci)%get(idx)  ! D = mobility in normalized units

            print *, "DEBUG Robin BC EDGE ", i, " (contact->interior):"
            print *, "  Contact vertex idx1 =", idx1, "(ict=", ict1, ")"
            print *, "  Interior vertex idx2 =", idx2
            print *, "  Contact type =", par%contacts(ict1)%type
            print *, "  Contact name =", trim(par%contacts(ict1)%name)
            print *, "  Carrier index ci =", ci
            print *, "  n0 (equilibrium) =", n0
            print *, "  S (thermionic) =", S, "denorm =", denorm(S, "cm/s"), "cm/s"
            print *, "  D (diffusion) =", D, "denorm =", denorm(D, "cm^2/s"), "cm²/s"
            print *, "  dx (edge length) =", dx, "denorm =", denorm(dx, "cm"), "cm"
            print *, "  surf (edge area) =", surf, "denorm =", denorm(surf, "cm^2"), "cm²"
            print *, "  D/dx =", D/dx
            print *, "  D/dx + S =", D/dx + S
            print *, "  Matrix[idx1,edge] = -(D/dx + S)*surf =", -(D/dx + S) * surf
            print *, "  Matrix[idx2,edge] = (D/dx)*surf =", (D/dx) * surf
            print *, "  RHS b_edge = -S*n0*surf =", -S * n0 * surf
            print *, ""

            ! Robin BC: -D(n2 - n1)/dx = S(n1 - n0)
            ! Rearrange: D/dx * n2 - (D/dx + S) * n1 = -S * n0
            call this%jaco_cdens(idx_dir)%p%set(idx1, idx, -(D/dx + S) * surf)  ! contact vertex
            call this%jaco_cdens(idx_dir)%p%set(idx2, idx,  (D/dx) * surf)      ! interior vertex

            ! store RHS contribution
            this%b_edge(i, idx_dir) = -S * n0 * surf
          end if

        else if (ict2 > 0 .and. ict1 == 0) then
          if (par%contacts(ict2)%type == CT_SCHOTTKY) then
            ! Interior (vertex 1) to Schottky contact (vertex 2): Robin BC
            n0 = calculate_equilibrium_density(par, ci, ict2)
            S = calculate_thermionic_velocity(par, ci, ict2)
            D = par%mob0(IDX_EDGE, idx_dir, ci)%get(idx)  ! D = mobility in normalized units

            print *, "DEBUG Robin BC EDGE ", i, " (interior->contact):"
            print *, "  Interior vertex idx1 =", idx1
            print *, "  Contact vertex idx2 =", idx2, "(ict=", ict2, ")"
            print *, "  Contact type =", par%contacts(ict2)%type
            print *, "  Contact name =", trim(par%contacts(ict2)%name)
            print *, "  Carrier index ci =", ci
            print *, "  n0 (equilibrium) =", n0
            print *, "  S (thermionic) =", S, "denorm =", denorm(S, "cm/s"), "cm/s"
            print *, "  D (diffusion) =", D, "denorm =", denorm(D, "cm^2/s"), "cm²/s"
            print *, "  dx (edge length) =", dx, "denorm =", denorm(dx, "cm"), "cm"
            print *, "  surf (edge area) =", surf, "denorm =", denorm(surf, "cm^2"), "cm²"
            print *, "  Matrix[idx1,edge] = (D/dx)*surf =", (D/dx) * surf
            print *, "  Matrix[idx2,edge] = -(D/dx + S)*surf =", -(D/dx + S) * surf
            print *, "  RHS b_edge = +S*n0*surf =", S * n0 * surf
            print *, ""

            ! Robin BC with opposite orientation
            call this%jaco_cdens(idx_dir)%p%set(idx1, idx, (D/dx) * surf)       ! interior vertex
            call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -(D/dx + S) * surf)  ! contact vertex

            ! store RHS contribution (opposite sign due to orientation)
            this%b_edge(i, idx_dir) = S * n0 * surf
          end if

        else
          ! other cases: keep standard behavior for uncontacted vertices
          if (ict1 == 0) call this%jaco_cdens(idx_dir)%p%set(idx1, idx,  surf)
          if (ict2 == 0) call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -surf)
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

    ! dirichlet conditions
    allocate (this%b(this%f%n), source = 0.0)
    j = par%transport_vct(0)%n
    do ict = 1, par%nct
      select case (par%contacts(ict)%type)
      case (CT_SCHOTTKY)
        ! CRITICAL FIX: Do NOT apply Dirichlet BC for Schottky contacts!
        ! The Robin BC from edge assembly fully determines the contact density
        ! Applying both Dirichlet and Robin creates an over-constrained system
        print *, "DEBUG: Schottky contact ", ict, " - NO Dirichlet BC"
        print *, "  Contact name: ", trim(par%contacts(ict)%name)
        print *, "  Robin BC from edges will determine contact density"
        do i = 1, par%transport_vct(ict)%n
          j = j + 1
          idx1 = par%transport_vct(ict)%get_idx(i)
          ! CRITICAL: Do NOT set matrix row to identity!
          ! Do NOT fix the contact density!
          ! The Robin BC edge assembly handles everything
          print *, "    Contact vertex ", i, ": idx1=", idx1, ", equation j=", j
          print *, "    (no Dirichlet constraint applied)"
        end do
      case default  ! CT_OHMIC, CT_GATE
        ! Apply Dirichlet BC for Ohmic/Gate contacts
        do i = 1, par%transport_vct(ict)%n
          j = j + 1
          idx1 = par%transport_vct(ict)%get_idx(i)
          call this%jaco_dens%set(idx1, idx1, 1.0)
          ! Ohmic contact: infinite injection (equilibrium density)
          if ((par%smc%dos == DOS_PARABOLIC) .and. (par%smc%dist == DIST_MAXWELL)) then
            this%b(j) = sqrt(par%smc%edos(1) * par%smc%edos(2)) * exp(- CR_CHARGE(ci) * par%contacts(ict)%phims - 0.5 * par%smc%band_gap)
          else
            call par%smc%get_dist(- CR_CHARGE(ci) * (par%contacts(ict)%phims - par%smc%band_edge(ci)), 0, F, dF)
            this%b(j) = par%smc%edos(ci) * F
          end if
        end do
      end select
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine update_robin_bc(this)
    !! Update Robin BC edge contributions based on current potential
    !! This makes the BC bias-dependent for Schottky contacts
    class(continuity), intent(inout) :: this

    integer :: idx_dir, i, ci, ict1, ict2
    integer, allocatable :: idx(:), idx1(:), idx2(:)
    logical :: status
    real :: n0, S, D, dx, surf, V_contact

    ci = this%ci
    allocate(idx(this%par%g%idx_dim), idx1(this%par%g%idx_dim), idx2(this%par%g%idx_dim))

    ! Loop over all edges to update Robin BC contributions
    do idx_dir = 1, this%par%g%idx_dim
      do i = 1, this%par%transport(IDX_EDGE,idx_dir)%n
        idx = this%par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        ict1 = this%par%ict%get(idx1)
        ict2 = this%par%ict%get(idx2)

        ! Update Robin BC for Schottky contacts
        ! Only process edges connecting Schottky contact to interior vertex
        if (ict1 > 0 .and. ict2 == 0) then
          if (this%par%contacts(ict1)%type == CT_SCHOTTKY) then
            surf = this%par%tr_surf(idx_dir)%get(idx)
            dx = this%par%g%get_len(idx, idx_dir)
            D = this%par%mob0(IDX_EDGE, idx_dir, ci)%get(idx)

            ! Get current voltage at contact vertex 1
            V_contact = this%pot%get(idx1)
            ! Calculate bias-dependent equilibrium density
            n0 = calculate_equilibrium_density_bias(this%par, ci, ict1, V_contact)
            ! Calculate bias-dependent thermionic velocity
            S = calculate_thermionic_velocity_bias(this%par, ci, ict1, V_contact)
            ! Update RHS contribution
            this%b_edge(i, idx_dir) = -S * n0 * surf
          end if
        else if (ict2 > 0 .and. ict1 == 0) then
          if (this%par%contacts(ict2)%type == CT_SCHOTTKY) then
            surf = this%par%tr_surf(idx_dir)%get(idx)
            dx = this%par%g%get_len(idx, idx_dir)
            D = this%par%mob0(IDX_EDGE, idx_dir, ci)%get(idx)

            ! Get current voltage at contact vertex 2
            V_contact = this%pot%get(idx2)
            ! Calculate bias-dependent equilibrium density
            n0 = calculate_equilibrium_density_bias(this%par, ci, ict2, V_contact)
            ! Calculate bias-dependent thermionic velocity
            S = calculate_thermionic_velocity_bias(this%par, ci, ict2, V_contact)
            ! Update RHS contribution (opposite sign)
            this%b_edge(i, idx_dir) = S * n0 * surf
          end if
        end if
      end do
    end do
  end subroutine

  subroutine continuity_eval(this)
    !! evaluate continuity equation
    class(continuity), intent(inout) :: this

    integer           :: idx_dir, i, j, ict1, ict2, j1, j2
    integer, allocatable :: idx(:), idx1(:), idx2(:)
    logical           :: status
    real, allocatable :: tmp(:)

    ! CRITICAL: Update Robin BC with current potential before evaluation
    ! This makes Schottky contacts bias-dependent
    call update_robin_bc(this)

    allocate (tmp(this%f%n))
    allocate (idx(this%par%g%idx_dim), idx1(this%par%g%idx_dim), idx2(this%par%g%idx_dim))

    call this%jaco_dens%matr%mul_vec(this%dens%get(), tmp)
    do idx_dir = 1, this%par%g%idx_dim
      call this%jaco_cdens(idx_dir)%p%matr%mul_vec(this%cdens(idx_dir)%get(), tmp, fact_y = 1.0)
    end do
    if (this%par%smc%incomp_ion) call this%jaco_genrec%matr%mul_vec(this%genrec%get(), tmp, fact_y = 1.0)

    ! Add edge contributions from Robin BC for Schottky contacts
    print *, "DEBUG: Adding Robin BC edge contributions to residual"
    do idx_dir = 1, this%par%g%idx_dim
      do i = 1, this%par%transport(IDX_EDGE,idx_dir)%n
        if (abs(this%b_edge(i, idx_dir)) > 1.0e-15) then  ! Only process if there's a contribution
          idx = this%par%transport(IDX_EDGE,idx_dir)%get_idx(i)
          call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
          call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

          ict1 = this%par%ict%get(idx1)
          ict2 = this%par%ict%get(idx2)

          print *, "  Edge ", i, ": b_edge =", this%b_edge(i, idx_dir)
          print *, "    idx1=", idx1, ", ict1=", ict1
          print *, "    idx2=", idx2, ", ict2=", ict2

          ! CRITICAL FIX: Only add contribution to the INTERIOR (uncontacted) vertex
          ! The Robin BC contribution should only affect the equation of the interior vertex
          ! Not both vertices (which was causing double-counting)

          if (ict1 == 0) then
            ! Vertex 1 is interior, vertex 2 is contact
            j1 = this%dens%itab%get(idx1)
            if (j1 > 0 .and. j1 <= this%par%transport_vct(0)%n) then
              print *, "    Adding to interior vertex 1: j1=", j1
              print *, "    tmp(j1) before:", tmp(j1)
              print *, "    b_edge contribution:", -this%b_edge(i, idx_dir)
              tmp(j1) = tmp(j1) - this%b_edge(i, idx_dir)
              print *, "    tmp(j1) after:", tmp(j1)
            end if
          else if (ict2 == 0) then
            ! Vertex 2 is interior, vertex 1 is contact
            j2 = this%dens%itab%get(idx2)
            if (j2 > 0 .and. j2 <= this%par%transport_vct(0)%n) then
              print *, "    Adding to interior vertex 2: j2=", j2
              print *, "    tmp(j2) before:", tmp(j2)
              print *, "    b_edge contribution:", -this%b_edge(i, idx_dir)
              tmp(j2) = tmp(j2) - this%b_edge(i, idx_dir)
              print *, "    tmp(j2) after:", tmp(j2)
            end if
          end if
          ! If both vertices are contacted or both uncontacted, no Robin BC contribution
        end if
      end do
    end do

    call this%f%set(tmp - this%b)
  end subroutine

  function calculate_equilibrium_density(par, ci, ict) result(n0)
    !! Calculate equilibrium carrier density at Schottky contact
    !! n₀ = N_c * exp(-q*Φ_B/kT)
    type(device_params), intent(in) :: par
    integer,             intent(in) :: ci   ! carrier index (CR_ELEC or CR_HOLE)
    integer,             intent(in) :: ict  ! contact index
    real                            :: n0

    real :: phi_B

    ! barrier_height is already in normalized units from input file parsing
    phi_B = par%contacts(ict)%barrier_height

    ! Schottky barrier reduces carrier density for both electrons and holes
    ! n₀ = N_c * exp(-Φ_B/kT) where Φ_B > 0 is the barrier height
    if ((par%smc%dos == DOS_PARABOLIC) .and. (par%smc%dist == DIST_MAXWELL)) then
      ! The barrier always reduces density, regardless of carrier type
      n0 = par%smc%edos(ci) * exp(-phi_B)
    else
      ! For general distribution, barrier reduces occupancy
      call par%smc%get_dist(-phi_B, 0, n0, phi_B)  ! use phi_B as dummy for dF
      n0 = par%smc%edos(ci) * n0
    end if
  end function

  function calculate_thermionic_velocity(par, ci, ict) result(S)
    !! Calculate thermionic emission velocity
    !! S = (A*T²)/(q*N_c) - NO barrier factor here!
    !! The barrier is included in n₀, not in S
    use normalization_m, only: norm, denorm

    type(device_params), intent(in) :: par
    integer,             intent(in) :: ci   ! carrier index
    integer,             intent(in) :: ict  ! contact index
    real                            :: S

    real :: A_star_phys, phi_B_norm, T, N_c, exp_factor
    real :: N_c_physical, S_physical, A_star_norm, phi_B_phys
    real, parameter :: ELEMENTARY_CHARGE = 1.602176634e-19  ! C (exact NIST value)
    logical, save :: first_call = .true.

    ! Get parameters (already normalized from input file parsing!)
    A_star_norm = par%contacts(ict)%richardson_const  ! Already normalized
    phi_B_norm = par%contacts(ict)%barrier_height     ! Already normalized
    T = par%T                                         ! K (physical)
    N_c = par%smc%edos(ci)                           ! normalized density

    ! CRITICAL FIX: Values from input file are ALREADY NORMALIZED
    ! Denormalize them to get physical values for display
    A_star_phys = denorm(A_star_norm, "A/cm^2/K^2")
    phi_B_phys = denorm(phi_B_norm, "eV")

    ! IMPORTANT: The thermionic emission velocity S does NOT include the barrier!
    ! The barrier effect is included in n₀ (equilibrium density), not in S
    ! S = (A*T²)/(q*N_c) is just the emission velocity coefficient
    exp_factor = 1.0  ! No barrier factor in velocity calculation!

    ! Denormalize N_c to physical units (cm⁻³)
    N_c_physical = denorm(N_c, "1/cm^3")

    if (first_call) then
      print *, "DEBUG calculate_thermionic_velocity:"
      print *, "  A_star (physical) =", A_star_phys, "A/cm²/K²"
      print *, "  A_star (normalized) =", A_star_norm
      print *, "  phi_B (physical) =", phi_B_phys, "eV"
      print *, "  phi_B (normalized) =", phi_B_norm
      print *, "  T =", T, "K"
      print *, "  N_c (normalized) =", N_c
      print *, "  N_c_physical =", N_c_physical, "cm⁻³"
      print *, "  CR_CHARGE(ci) =", CR_CHARGE(ci)
      print *, "  exp_factor =", exp_factor
      print *, "  ELEMENTARY_CHARGE =", ELEMENTARY_CHARGE, "C"
    end if

    ! Calculate thermionic velocity in physical units
    ! S = (A*T²)/(q*N_c) - WITHOUT barrier factor
    ! A* in A/cm²/K², T in K, exp_factor=1.0, EC in C, N_c in cm⁻³
    ! Units: [A/cm²/K²][K²]/([C][cm⁻³]) = [A/cm²]/[C/cm³] = [A·cm]/[C]
    ! Since A = C/s, we get [C·cm/s]/[C] = cm/s ✓
    S_physical = (A_star_phys * T**2 * exp_factor) / (ELEMENTARY_CHARGE * N_c_physical)  ! cm/s

    ! Normalize to simulator units
    S = norm(S_physical, "cm/s")

    if (first_call) then
      print *, "  S_physical =", S_physical, "cm/s"
      print *, "  S (normalized) =", S
      first_call = .false.
    end if

    ! Note: This properly handles pre-normalized input values
  end function

  function calculate_equilibrium_density_bias(par, ci, ict, V_contact) result(n0)
    !! Calculate equilibrium carrier density at Schottky contact with bias
    !! n₀(V) = N_c * exp(-(Φ_B - qV)/kT)
    type(device_params), intent(in) :: par
    integer,             intent(in) :: ci   ! carrier index (CR_ELEC or CR_HOLE)
    integer,             intent(in) :: ict  ! contact index
    real,                intent(in) :: V_contact ! electrostatic potential at contact
    real                            :: n0

    real :: phi_B, V_bias

    ! barrier_height is already in normalized units from input file parsing
    phi_B = par%contacts(ict)%barrier_height

    ! V_contact is the electrostatic potential (normalized)
    ! For bias-dependent emission, we need the voltage relative to equilibrium
    ! This is already included in the potential
    V_bias = V_contact

    ! Schottky barrier with bias: n₀ = N_c * exp(-(Φ_B - qV)/kT)
    ! In normalized units: n₀ = N_c * exp(-(Φ_B - V))
    if ((par%smc%dos == DOS_PARABOLIC) .and. (par%smc%dist == DIST_MAXWELL)) then
      ! Include bias effect: barrier is reduced by applied voltage
      n0 = par%smc%edos(ci) * exp(-(phi_B - CR_CHARGE(ci) * V_bias))
    else
      ! For general distribution with bias
      call par%smc%get_dist(-(phi_B - CR_CHARGE(ci) * V_bias), 0, n0, phi_B)
      n0 = par%smc%edos(ci) * n0
    end if

    ! Debug output for first call
    if (abs(V_bias) > 1.0e-6) then
      print *, "DEBUG: Bias-dependent n0:"
      print *, "  V_contact =", denorm(V_contact, "V")
      print *, "  V_bias =", denorm(V_bias, "V")
      print *, "  phi_B =", denorm(phi_B, "V")
      print *, "  n0(V) =", denorm(n0, "1/cm^3")
      print *, "  n0(V)/n0(0) =", exp(CR_CHARGE(ci) * V_bias)
    end if
  end function

  function calculate_thermionic_velocity_bias(par, ci, ict, V_contact) result(S)
    !! Calculate thermionic emission velocity with bias
    !! For now, S is voltage-independent (could add image force lowering later)
    use normalization_m, only: norm, denorm

    type(device_params), intent(in) :: par
    integer,             intent(in) :: ci   ! carrier index
    integer,             intent(in) :: ict  ! contact index
    real,                intent(in) :: V_contact ! electrostatic potential at contact
    real                            :: S

    ! For now, use voltage-independent thermionic velocity
    ! (Could add image force barrier lowering here if needed)
    S = calculate_thermionic_velocity(par, ci, ict)

    ! Future enhancement: include voltage-dependent barrier lowering
    ! due to image force effect (Schottky effect)
  end function

end module
