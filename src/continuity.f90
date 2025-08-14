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

    type(dirichlet_stencil)                :: st_dir
    type(near_neighb_stencil), allocatable :: st_nn(:)
    type(near_neighb_stencil)              :: st_vv     ! Single vertex-to-vertex stencil for Robin BC
    type(empty_stencil)                    :: st_em
    type(stencil_ptr), allocatable         :: st_cdens(:,:)  ! (0:nct, idx_dir)
    type(stencil_ptr), allocatable         :: st_dens(:)    ! Stencil array for jaco_dens

    type(jacobian),     pointer     :: jaco_dens     => null()
    type(jacobian),     pointer     :: jaco_dens_t   => null()
    type(jacobian_ptr), allocatable :: jaco_cdens(:)
    type(jacobian),     pointer     :: jaco_genrec   => null()
  contains
    procedure :: init => continuity_init
    procedure :: eval => continuity_eval
  end type

contains

  subroutine continuity_init(this, par, stat, dens, iref, cdens, genrec)
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

    integer              :: ci, i, ict, idx_dir, idens, idx_dim, igenrec, j, k, j_interior
    integer, allocatable :: idx(:), idx1(:), idx2(:), idx3(:), icdens(:)
    integer              :: ict1, ict2
    logical              :: status
    real                 :: surf, F, dF, dx, S, D, n0

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
    allocate (idx(idx_dim), idx1(idx_dim), idx2(idx_dim), idx3(idx_dim), icdens(idx_dim))
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
      call this%st_nn(idx_dir)%init(par%g, IDX_VERTEX, 0, IDX_EDGE, idx_dir)  ! for jaco_cdens (divergence)
    end do
    call this%st_em%init()

    ! Initialize vertex-to-vertex near-neighbor stencil for Robin BC
    ! This allows connections from vertices to their neighboring vertices
    ! Using IDX_VERTEX for both source and target with direction 0 means vertex-to-vertex connections
    call this%st_vv%init(par%g, IDX_VERTEX, 0, IDX_VERTEX, 0)

    ! DEBUG: Check if stencil actually allows vertex-to-vertex connections
    print *, "DEBUG: st_vv nmax =", this%st_vv%nmax
    print *, "DEBUG: VERTEX->VERTEX max neighb =", par%g%get_max_neighb(IDX_VERTEX,0,IDX_VERTEX,0)

    ! Test the stencil enumeration for the first vertex
    if (par%nct > 0 .and. par%transport_vct(1)%n > 0) then
      idx1 = par%transport_vct(1)%get_idx(1)
      print *, "DEBUG: Testing st_vv for vertex", idx1
      do k = 1, this%st_vv%nmax
        call this%st_vv%get(idx1, k, idx2, status)
        if (status) then
          print *, "  Neighbor", k, ":", idx2, "valid=", status
        end if
      end do
    end if

    ! dependencies
    idens = this%depend(this%dens)
    do idx_dir = 1, idx_dim
      icdens(idx_dir) = this%depend(this%cdens(idx_dir))
    end do
    if (par%smc%incomp_ion) igenrec = this%depend(this%genrec)

    ! init jacobians
    ! Build per-contact stencil array for jaco_dens
    allocate(this%st_dens(0:par%nct))
    this%st_dens(0) = this%st_em%get_ptr()  ! Interior vertices use empty stencil
    do ict = 1, par%nct
      select case (par%contacts(ict)%type)
      case (CT_SCHOTTKY)
        ! Schottky contacts need vertex-to-vertex stencil for Robin BC
        ! This allows connections from contact vertices to neighboring interior vertices
        this%st_dens(ict) = this%st_vv%get_ptr()
      case default  ! CT_OHMIC, CT_GATE
        ! Ohmic/Gate contacts use dirichlet stencil (diagonal only)
        this%st_dens(ict) = this%st_dir%get_ptr()
      end select
    end do

    this%jaco_dens   => this%init_jaco_f(idens, &
      & st = this%st_dens, const = .true., dtime = .false.)  ! const=.true. since st_vv has nmax=3 (sufficient)
    if (.not. stat) then
      this%jaco_dens_t => this%init_jaco_f(idens, &
        & st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .true. )
    end if
    ! Build stencil arrays based on contact types
    ! Allocate stencil array as component to ensure proper lifetime
    allocate(this%st_cdens(0:par%nct, idx_dim))

    do idx_dir = 1, idx_dim
      ! Uncontacted vertices always use near_neighb stencil
      this%st_cdens(0, idx_dir) = this%st_nn(idx_dir)%get_ptr()

      ! For each contact, choose stencil based on contact type
      do ict = 1, par%nct
        select case (par%contacts(ict)%type)
        case (CT_SCHOTTKY)
          ! Schottky contacts need edge contributions for Robin BC
          this%st_cdens(ict, idx_dir) = this%st_nn(idx_dir)%get_ptr()
        case (CT_OHMIC, CT_GATE)
          ! Ohmic/Gate contacts use empty stencil (no edge contributions)
          this%st_cdens(ict, idx_dir) = this%st_em%get_ptr()
        case default
          ! Default to empty for unknown contact types
          this%st_cdens(ict, idx_dir) = this%st_em%get_ptr()
        end select
      end do

      this%jaco_cdens(idx_dir)%p => this%init_jaco_f(icdens(idx_dir), &
        & st = this%st_cdens(:, idx_dir), const = .true., dtime = .false.)
    end do
    if (par%smc%incomp_ion) then
      this%jaco_genrec => this%init_jaco_f(igenrec, &
        & st = [this%st_dir%get_ptr(), (this%st_em%get_ptr(), ict = 1, par%nct)], &
        & const = .true., dtime = .false.)
    end if

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
          ! Contact at vertex 1, interior at vertex 2
          ! Standard one-sided flux (only interior vertex gets contribution)
          ! Contact vertex boundary condition is handled separately
          call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -surf)  ! interior vertex only

        else if (ict2 > 0 .and. ict1 == 0) then
          ! Interior at vertex 1, contact at vertex 2
          ! Standard one-sided flux (only interior vertex gets contribution)
          ! Contact vertex boundary condition is handled separately
          call this%jaco_cdens(idx_dir)%p%set(idx1, idx, surf)  ! interior vertex only

        else
          ! other cases: keep standard behavior for uncontacted vertices
          if (ict1 == 0) call this%jaco_cdens(idx_dir)%p%set(idx1, idx,  surf)
          if (ict2 == 0) call this%jaco_cdens(idx_dir)%p%set(idx2, idx, -surf)
        end if
      end do
    end do

    ! Allocate RHS before any assembly
    allocate (this%b(this%f%n), source = 0.0)

    ! Apply Robin BC for Schottky contacts
    ! Robin BC: (D/dx)(n_i - n_c) - S(n_c - n_0) = 0
    ! Rearranged: -(D/dx + S)*n_c + (D/dx)*n_i + S*n_0 = 0
    j = par%transport_vct(0)%n
    do ict = 1, par%nct
      if (par%contacts(ict)%type == CT_SCHOTTKY) then
        ! Calculate Robin BC parameters for this contact
        n0 = calculate_equilibrium_density(par, ci, ict)
        S = calculate_thermionic_velocity(par, ci, ict)

        do i = 1, par%transport_vct(ict)%n
          j = j + 1
          idx1 = par%transport_vct(ict)%get_idx(i)  ! Contact vertex

          print *, "DEBUG: Processing Schottky contact vertex:", idx1

          ! Find connected edges and add Robin BC contributions
          do idx_dir = 1, idx_dim
            do k = 1, par%transport(IDX_EDGE, idx_dir)%n
              idx = par%transport(IDX_EDGE, idx_dir)%get_idx(k)

              ! Get edge vertices
              call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx2, status)
              if (.not. status) cycle
              call par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx3, status)
              if (.not. status) cycle

              ! print *, "  Edge", k, "connects vertices:", idx2, "and", idx3
              ! print *, "  Contact vertex idx1:", idx1
              ! print *, "  idx1==idx2?", all(idx1 == idx2), "idx1==idx3?", all(idx1 == idx3)
              ! print *, "  ict(idx2)=", par%ict%get(idx2), "ict(idx3)=", par%ict%get(idx3)

              ! Check if this edge connects our contact to an interior vertex
              if (all(idx1 == idx2) .and. par%ict%get(idx3) == 0) then
                ! Our contact (idx1) is at vertex 1 of edge, interior at vertex 2
                print *, "  Found edge connecting Schottky to interior (case 1)"
                dx = par%g%get_len(idx, idx_dir)
                D = par%mob0(IDX_EDGE, idx_dir, ci)%get(idx)  ! Diffusivity = mobility in normalized units
                surf = par%tr_surf(idx_dir)%get(idx)  ! Get surface factor for flux calculation

                print *, "    D=", D, "dx=", dx, "D/dx=", D/dx, "S=", S, "surf=", surf

                ! Add Robin BC contributions to matrix with surface factor
                ! Use add method to accumulate contributions (multiple edges may connect to same contact)
                call this%jaco_dens%add(idx1, idx1, -(D/dx + S) * surf)
                call this%jaco_dens%add(idx1, idx3, (D/dx) * surf)

                ! Add RHS contribution for this edge
                j_interior = this%dens%itab%get(idx1)  ! Get equation index for contact vertex
                if (j_interior > 0) then
                  this%b(j_interior) = this%b(j_interior) - S * n0 * surf  ! MINUS sign: RHS = -S*n0*surf
                end if

              else if (all(idx1 == idx3) .and. par%ict%get(idx2) == 0) then
                ! Our contact (idx1) is at vertex 2 of edge, interior at vertex 1
                print *, "  Found edge connecting Schottky to interior (case 2)"
                dx = par%g%get_len(idx, idx_dir)
                D = par%mob0(IDX_EDGE, idx_dir, ci)%get(idx)  ! Diffusivity = mobility in normalized units
                surf = par%tr_surf(idx_dir)%get(idx)  ! Get surface factor for flux calculation

                print *, "    D=", D, "dx=", dx, "D/dx=", D/dx, "S=", S, "surf=", surf

                ! Add Robin BC contributions to matrix with surface factor
                ! Use add method to accumulate contributions (multiple edges may connect to same contact)
                call this%jaco_dens%add(idx1, idx1, -(D/dx + S) * surf)
                call this%jaco_dens%add(idx1, idx2, (D/dx) * surf)

                ! Add RHS contribution for this edge
                j_interior = this%dens%itab%get(idx1)  ! Get equation index for contact vertex
                if (j_interior > 0) then
                  this%b(j_interior) = this%b(j_interior) - S * n0 * surf  ! MINUS sign: RHS = -S*n0*surf
                end if
              end if
            end do
          end do

          ! RHS contributions are now added per-edge in the loop above
        end do
      else
        ! Non-Schottky contacts: just increment j
        j = j + par%transport_vct(ict)%n
      end if
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
    j = par%transport_vct(0)%n
    do ict = 1, par%nct
      select case (par%contacts(ict)%type)
      case (CT_SCHOTTKY)
        ! Schottky contacts: Add regularization to prevent singular matrix
        ! The Robin BC contributions added above will dominate this regularization
        do i = 1, par%transport_vct(ict)%n
          j = j + 1
          idx1 = par%transport_vct(ict)%get_idx(i)

          ! Add diagonal regularization to ensure matrix is non-singular
          ! This will be dominated by the Robin BC terms added earlier
          ! Increased from 1e-10 to 1e-6 for better numerical stability
          call this%jaco_dens%add(idx1, idx1, 1.0e-6)

          ! Set initial guess for contact density
          ! Use equilibrium density as starting point
          n0 = calculate_equilibrium_density(par, ci, ict)
          ! Clamp n0 to prevent numerical issues with extremely small values
          n0 = max(n0, 1.0e-20)
          ! Get the actual equation index for this contact vertex
          j_interior = this%dens%itab%get(idx1)
          if (j_interior > 0) then
            this%b(j_interior) = this%b(j_interior) - 1.0e-6 * n0  ! Use same sign as Robin BC
          end if
        end do
      case default  ! CT_OHMIC, CT_GATE
        ! Apply Dirichlet BC for Ohmic/Gate contacts
        do i = 1, par%transport_vct(ict)%n
          j = j + 1
          idx1 = par%transport_vct(ict)%get_idx(i)
          call this%jaco_dens%set(idx1, idx1, 1.0)
          ! Get the actual equation index for this contact vertex
          j_interior = this%dens%itab%get(idx1)
          if (j_interior > 0) then
            ! Ohmic contact: infinite injection (equilibrium density)
            if ((par%smc%dos == DOS_PARABOLIC) .and. (par%smc%dist == DIST_MAXWELL)) then
              this%b(j_interior) = sqrt(par%smc%edos(1) * par%smc%edos(2)) * exp(- CR_CHARGE(ci) * par%contacts(ict)%phims - 0.5 * par%smc%band_gap)
            else
              call par%smc%get_dist(- CR_CHARGE(ci) * (par%contacts(ict)%phims - par%smc%band_edge(ci)), 0, F, dF)
              this%b(j_interior) = par%smc%edos(ci) * F
            end if
          end if
        end do
      end select
    end do

    ! finish initialization
    call this%init_final()
  end subroutine


  subroutine continuity_eval(this)
    !! evaluate continuity equation
    class(continuity), intent(inout) :: this

    integer           :: idx_dir, i, j, ict1, ict2, j1, j2
    integer, allocatable :: idx(:), idx1(:), idx2(:)
    logical           :: status
    real, allocatable :: tmp(:)

    allocate (tmp(this%f%n))
    allocate (idx(this%par%g%idx_dim), idx1(this%par%g%idx_dim), idx2(this%par%g%idx_dim))

    call this%jaco_dens%matr%mul_vec(this%dens%get(), tmp)
    do idx_dir = 1, this%par%g%idx_dim
      call this%jaco_cdens(idx_dir)%p%matr%mul_vec(this%cdens(idx_dir)%get(), tmp, fact_y = 1.0)
    end do
    if (this%par%smc%incomp_ion) call this%jaco_genrec%matr%mul_vec(this%genrec%get(), tmp, fact_y = 1.0)

    ! Robin BC contributions are now handled in the density block during init
    ! No need for separate edge contributions here

    call this%f%set(tmp - this%b)

    ! Check Robin BC residual for debugging (only for Schottky contacts)
    call check_robin_bc_residual(this)
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

  subroutine check_robin_bc_residual(this)
    !! Check if Robin BC is satisfied at Schottky contacts (for debugging)
    class(continuity), intent(in) :: this

    integer              :: idx_dir, i, ict1, ict2, ci, j_contact, j_interior
    integer, allocatable :: idx(:), idx1(:), idx2(:)
    real                 :: n_contact, n_interior, n0, S, D, dx, surf
    real                 :: Rc, D_over_dx, D_plus_S_over_dx
    real, allocatable   :: dens_array(:)
    logical              :: status, printed_header

    ci = this%ci
    printed_header = .false.
    allocate(idx(this%par%g%idx_dim), idx1(this%par%g%idx_dim), idx2(this%par%g%idx_dim))

    ! Get density data array
    allocate(dens_array(this%dens%n))
    dens_array = this%dens%get()

    ! Loop through edges to find Schottky contacts
    do idx_dir = 1, this%par%g%idx_dim
      do i = 1, this%par%transport(IDX_EDGE,idx_dir)%n
        ! Get edge and vertex indices
        idx = this%par%transport(IDX_EDGE,idx_dir)%get_idx(i)
        call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 1, idx1, status)
        call this%par%g%get_neighb(IDX_EDGE, idx_dir, IDX_VERTEX, 0, idx, 2, idx2, status)

        ict1 = this%par%ict%get(idx1)
        ict2 = this%par%ict%get(idx2)

        ! Check for Schottky contact edges
        if (ict1 > 0 .and. ict2 == 0) then
          if (this%par%contacts(ict1)%type == CT_SCHOTTKY) then
            ! Get edge parameters
            surf = this%par%tr_surf(idx_dir)%get(idx)
            dx = this%par%g%get_len(idx, idx_dir)
            D = this%par%mob0(IDX_EDGE, idx_dir, ci)%get(idx)
            ! Contact at vertex 1, interior at vertex 2
            n0 = calculate_equilibrium_density(this%par, ci, ict1)
            S = calculate_thermionic_velocity(this%par, ci, ict1)
            ! Get density values at vertices using the equation indices
            j_contact = this%dens%itab%get(idx1)
            j_interior = this%dens%itab%get(idx2)
            if (j_contact > 0) n_contact = dens_array(j_contact)
            if (j_interior > 0) n_interior = dens_array(j_interior)

            ! Calculate residual: Rc = (D/dx + S)*n_contact - (D/dx)*n_interior - S*n0
            D_over_dx = D/dx
            D_plus_S_over_dx = D_over_dx + S
            Rc = D_plus_S_over_dx * n_contact - D_over_dx * n_interior - S * n0

            ! ! Print header once
            ! if (.not. printed_header) then
            !   print '(A)', "=== ROBIN BC RESIDUAL CHECK ==="
            !   print '(A)', "Contact: " // trim(this%par%contacts(ict1)%name)
            !   print '(A,ES12.3)', "  Fixed n0_B = ", denorm(n0,"cm^-3")
            !   print '(A,ES12.3,A)', "  S = ", S, " (normalized)"
            !   print '(A,ES12.3)', "  D/dx = ", D_over_dx
            !   printed_header = .true.
            ! end if

            ! Print residual if significant
            if (abs(Rc) > 1.0e-10) then
              print '(A,I0,A)', "Edge ", i, " (contact->interior):"
              print '(A,ES12.3)', "  n_contact = ", n_contact
              print '(A,ES12.3)', "  n_interior = ", n_interior
              print '(A,ES12.3)', "  n_contact/n0_B = ", n_contact/n0
              print '(A,ES12.3)', "  Residual Rc = ", Rc
            end if
          end if
        else if (ict2 > 0 .and. ict1 == 0) then
          if (this%par%contacts(ict2)%type == CT_SCHOTTKY) then
            ! Get edge parameters
            surf = this%par%tr_surf(idx_dir)%get(idx)
            dx = this%par%g%get_len(idx, idx_dir)
            D = this%par%mob0(IDX_EDGE, idx_dir, ci)%get(idx)
            ! Contact at vertex 2, interior at vertex 1
            n0 = calculate_equilibrium_density(this%par, ci, ict2)
            S = calculate_thermionic_velocity(this%par, ci, ict2)
            ! Get density values at vertices using the equation indices
            j_contact = this%dens%itab%get(idx2)
            j_interior = this%dens%itab%get(idx1)
            if (j_contact > 0) n_contact = dens_array(j_contact)
            if (j_interior > 0) n_interior = dens_array(j_interior)

            ! Calculate residual with opposite orientation
            D_over_dx = D/dx
            D_plus_S_over_dx = D_over_dx + S
            Rc = D_plus_S_over_dx * n_contact - D_over_dx * n_interior - S * n0

            ! Print header once
            ! if (.not. printed_header) then
            !   print '(A)', "=== ROBIN BC RESIDUAL CHECK ==="
            !   print '(A)', "Contact: " // trim(this%par%contacts(ict2)%name)
            !   print '(A,ES12.3)', "  Fixed n0_B = ", denorm(n0,"cm^-3")
            !   print '(A,ES12.3,A)', "  S = ", S, " (normalized)"
            !   print '(A,ES12.3)', "  D/dx = ", D_over_dx
            !   printed_header = .true.
            ! end if

            ! Print residual if significant
            if (abs(Rc) > 1.0e-10) then
              print '(A,I0,A)', "Edge ", i, " (interior->contact):"
              print '(A,ES12.3)', "  n_contact = ", n_contact
              print '(A,ES12.3)', "  n_interior = ", n_interior
              print '(A,ES12.3)', "  n_contact/n0_B = ", n_contact/n0
              print '(A,ES12.3)', "  Residual Rc = ", Rc
            end if
          end if
        end if
      end do
    end do

    if (printed_header) then
      print '(A)', ""
    end if

    deallocate(idx, idx1, idx2, dens_array)
  end subroutine


end module
