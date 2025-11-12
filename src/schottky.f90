module schottky_m

  use device_params_m,  only: device_params
  use normalization_m,  only: norm, denorm
  use semiconductor_m,  only: CR_ELEC, CR_HOLE, CR_NAME
  use grid_m,           only: IDX_VERTEX, IDX_EDGE, IDX_CELL
  use variable_m,       only: variable_real
  use equation_m,       only: equation
  use jacobian_m,       only: jacobian, jacobian_ptr
  use electric_field_m, only: electric_field
  use grid_data_m,      only: grid_data1_real, grid_data2_real, grid_data3_real
  use vselector_m,      only: vselector
  use stencil_m,        only: dirichlet_stencil, empty_stencil
  use error_m,          only: program_error
  use contact_m,        only: CT_SCHOTTKY
  use math_m,           only: PI
  use quad_m,           only: quad

  implicit none

  private
  public :: schottky_injection_mb, schottky_velocity
  public :: schottky_injection_mb_bias, schottky_barrier_lowering
  public :: get_schottky_contact_normal_dir
  public :: schottky_injection, calc_schottky_injection
  public :: barrier_lowering

  type, extends(variable_real) :: barrier_lowering
    !! Barrier lowering delta_phi_b scalar field at vertices
  contains
    procedure :: init => barrier_lowering_init
  end type

  type, extends(variable_real) :: schottky_injection
    !! Bias-dependent injection density n0B at Schottky contacts
    !! Only exists at Schottky contact vertices

    integer :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => schottky_injection_init
  end type

  type, extends(equation) :: calc_schottky_injection
    !! Calculate bias-dependent injection density n0B from electric field
    !! n0B = n0B_base * exp(delta_phi_b(E)) for electrons
    !! n0B = n0B_base * exp(-delta_phi_b(E)) for holes

    type(device_params), pointer :: par => null()
      !! pointer to device parameters
    integer :: ci
      !! carrier index

    type(electric_field), pointer :: efield(:) => null()
      !! electric field components (array of size dim)
    type(vselector) :: n0b
      !! injection density variable

    integer, allocatable :: contact_normal(:)
      !! normal direction for each contact (0 if not Schottky)

    type(dirichlet_stencil) :: st_dir
      !! stencil for contact vertices
    type(empty_stencil) :: st_em
      !! empty stencil

    type(jacobian_ptr), allocatable :: jaco_efield(:)
      !! jacobian for E-field dependencies (one per direction)

    type(barrier_lowering), pointer :: delta_phi_b_var => null()
      !! barrier lowering variable (for output)
  contains
    procedure :: init => calc_schottky_injection_init
    procedure :: eval => calc_schottky_injection_eval
  end type

contains

  function get_schottky_contact_normal_dir(par, ict) result(normal_dir)
    !! Determine the normal direction of a Schottky contact for barrier lowering
    !! Analyzes vertex coordinates to find the constant dimension
    !! Returns: 1 for x-normal, 2 for y-normal, 3 for z-normal
    !! Returns: 0 if contact is not axis-aligned

    type(device_params), intent(in) :: par
    integer,             intent(in) :: ict     ! Contact index
    integer                         :: normal_dir

    integer :: i, dir, idx(par%g%idx_dim)
    integer :: n_vertices
    real    :: coord_min(3), coord_max(3), coord_range(3)
    real    :: vertex_pos(3)
    real    :: tolerance, min_range, char_length

    ! Initialize
    normal_dir = 0
    coord_min = 1e30
    coord_max = -1e30
    tolerance = 1e-6  ! Relative tolerance for planarity

    ! Get number of transport vertices for this contact
    n_vertices = par%transport_vct(ict)%n

    if (n_vertices == 0) then
      print *, "WARNING: Schottky contact ", ict, " has no transport vertices"
      return
    end if

    ! Analyze coordinate ranges
    do i = 1, n_vertices
      idx = par%transport_vct(ict)%get_idx(i)
      call par%g%get_vertex(idx, vertex_pos(1:par%g%dim))

      do dir = 1, par%g%dim
        coord_min(dir) = min(coord_min(dir), vertex_pos(dir))
        coord_max(dir) = max(coord_max(dir), vertex_pos(dir))
      end do
    end do

    ! Calculate range in each direction
    do dir = 1, par%g%dim
      coord_range(dir) = coord_max(dir) - coord_min(dir)
    end do

    ! Find direction with minimum range (normal to contact)
    min_range = huge(1.0)
    do dir = 1, par%g%dim
      if (coord_range(dir) < min_range) then
        min_range = coord_range(dir)
        normal_dir = dir
      end if
    end do

    ! Validate planarity
    if (par%g%dim > 1) then
      char_length = maxval(coord_range(1:par%g%dim))
      if (min_range > tolerance * char_length .and. char_length > 0.0) then
        print *, "WARNING: Schottky contact ", ict, " not perfectly axis-aligned for barrier lowering"
        print *, "  Range in normal direction: ", denorm(min_range, "nm"), " nm"
      end if
    end if

  end function

  subroutine schottky_injection_mb(par, ci, ict, ninj)
    !! Calculate equilibrium density n0 for Schottky contact Robin BC
    !! This is the thermionic emission density at zero bias
    !! Electrons: n0 = Nc * exp(-phi_Bn)
    !! Holes:     p0 = Nv * exp(-phi_Bp) where phi_Bp = E_g - phi_Bn

    type(device_params), intent(in)  :: par
    integer,             intent(in)  :: ci      ! Carrier index (CR_ELEC or CR_HOLE)
    integer,             intent(in)  :: ict     ! Contact index
    real,                intent(out) :: ninj    ! Equilibrium density n0 (normalized)

    real :: phi_Bn, phi_Bp
    logical, save :: first_call = .true.
    integer, save :: call_count = 0

    ! Get normalized barrier height (already converted in device_params)
    phi_Bn = par%contacts(ict)%phi_b

    if (ci == CR_ELEC) then
      ! Electrons: n0 = Nc * exp(-phi_Bn)
      ninj = par%smc%edos(CR_ELEC) * exp(-phi_Bn)

      ! Debug output - only print first few calls and then periodically
      call_count = call_count + 1
      if (first_call .or. mod(call_count, 100) == 0) then
        print "(A,A,A,I2,A)", "DEBUG_SCHOTTKY: Contact ", trim(par%contacts(ict)%name), &
                              " (", ict, ") n0B calculation:"
        print "(A,ES14.6,A)", "  phi_b (denorm) = ", denorm(phi_Bn, "eV"), " eV"
        print "(A,ES14.6,A)", "  phi_b (norm)   = ", phi_Bn, " (normalized)"
        print "(A,ES14.6,A)", "  Nc (denorm)    = ", denorm(par%smc%edos(CR_ELEC), "cm^-3"), " cm^-3"
        print "(A,ES14.6,A)", "  Nc (norm)      = ", par%smc%edos(CR_ELEC), " (normalized)"
        print "(A,ES14.6)",   "  exp(-phi_Bn)   = ", exp(-phi_Bn)
        print "(A,ES14.6,A)", "  n0B (denorm)   = ", denorm(ninj, "cm^-3"), " cm^-3"
        print "(A,ES14.6,A)", "  n0B (norm)     = ", ninj, " (normalized)"
        first_call = .false.
      end if
    else  ! CR_HOLE
      ! Holes: barrier from valence band
      phi_Bp = par%smc%band_gap - phi_Bn
      ninj = par%smc%edos(CR_HOLE) * exp(-phi_Bp)
    end if
  end subroutine

  function schottky_velocity(par, ci, ict) result(s)
    !! Calculate thermionic emission velocity at Schottky contact
    !! Using Richardson constant: v_surf = A*T^2/(q*Nc)

    type(device_params), intent(in) :: par
    integer,             intent(in) :: ci   ! Carrier index
    integer,             intent(in) :: ict  ! Contact index
    real                            :: s
    logical, save :: first_call_v = .true.

    ! Check if Richardson constant is provided and > 0
    if (par%contacts(ict)%A_richardson > 0.0) then
      ! Calculate normalized surface velocity
      ! v_surf = A*T^2/(q*Nc) where q is handled by normalization
      ! T must be normalized, Nc is already normalized
      s = par%contacts(ict)%A_richardson * norm(par%T, "K") * norm(par%T, "K") / par%smc%edos(ci)

      if (first_call_v) then
        print "(A,A,A,I2,A)", "DEBUG_SCHOTTKY: Contact ", trim(par%contacts(ict)%name), &
                              " (", ict, ") surface velocity:"
        print "(A,ES14.6,A)", "  A_richardson  = ", par%contacts(ict)%A_richardson, " A/cm^2/K^2"
        print "(A,ES14.6,A)", "  Temperature    = ", par%T, " K"
        print "(A,ES14.6,A)", "  v_surf (denorm)= ", denorm(s,"cm/s"), " cm/s"
        print "(A,ES14.6,A)", "  v_surf (norm)  = ", s, " (normalized)"
        first_call_v = .false.
      end if
    else
      ! Default thermal velocity estimate (v_th/4)
      s = 0.25  ! v_th/4 in normalized units
      if (first_call_v) then
        print "(A,A)", "DEBUG_SCHOTTKY: Using default v_surf = ", "0.25 (v_th/4)"
        first_call_v = .false.
      end if
    end if
  end function

  subroutine schottky_injection_mb_bias(par, ci, ict, E_field, eps_r, ninj, dninj_dE)
    !! Calculate bias-dependent equilibrium density and derivative for Schottky contact
    !! Includes image force barrier lowering effect
    !! NOTE: This is a legacy function, not currently used in main flow.
    !! The field-dependent barrier lowering is handled in calc_schottky_injection_eval

    type(device_params), intent(in)  :: par
    integer,             intent(in)  :: ci      ! Carrier index (CR_ELEC or CR_HOLE)
    integer,             intent(in)  :: ict     ! Contact index
    real,                intent(in)  :: E_field ! Electric field at interface (normalized)
    real,                intent(in)  :: eps_r   ! Relative permittivity
    real,                intent(out) :: ninj    ! Injection density (normalized)
    real,                intent(out) :: dninj_dE ! Derivative d(ninj)/dE (normalized)

    real :: ninj_base, delta_phi_b, d_delta_phi_dE

    ! Get base injection without field
    call schottky_injection_mb(par, ci, ict, ninj_base)

    ! Calculate barrier lowering and its derivative
    call schottky_barrier_lowering(par, ict, E_field, eps_r, delta_phi_b, d_delta_phi_dE)

    ! Apply barrier lowering to injection density
    if (ci == CR_ELEC) then
      ! Electrons: lower barrier increases injection
      ninj = ninj_base * exp(delta_phi_b)
      dninj_dE = ninj * d_delta_phi_dE
    else  ! CR_HOLE
      ! Holes: opposite effect
      ninj = ninj_base * exp(-delta_phi_b)
      dninj_dE = -ninj * d_delta_phi_dE
    end if
  end subroutine

  subroutine schottky_barrier_lowering(par, ict, E_field, eps_r, delta_phi_b, d_delta_phi_dE)
    !! Calculate image force barrier lowering (Schottky effect)
    !! Implements Δφ_b = sqrt(q*|E|/(4π*ε))

    type(device_params), intent(in)  :: par
    integer,             intent(in)  :: ict           ! Contact index
    real,                intent(in)  :: E_field       ! Electric field (normalized)
    real,                intent(in)  :: eps_r         ! Relative permittivity
    real,                intent(out) :: delta_phi_b   ! Barrier lowering (normalized)
    real,                intent(out) :: d_delta_phi_dE ! Derivative d(Δφ_b)/dE

    real :: E_smooth, eps_smooth
    real :: gamma, q_norm
    integer, save :: debug_count = 0

    ! Check if IFBL is disabled for this contact
    if (.not. par%contacts(ict)%ifbl) then
      delta_phi_b = 0.0
      d_delta_phi_dE = 0.0
      return
    else
    end if

    ! Smoothing parameter to avoid kink at E=0
    eps_smooth = 1e-10

    ! Smooth |E| to avoid numerical issues: |E| ≈ sqrt(E² + ε²)
    E_smooth = sqrt(E_field**2 + eps_smooth**2)

    ! Normalized charge (already in normalized units)
    q_norm = 1.0

    ! Image force coefficient: γ = sqrt(q/(4π*ε₀*εᵣ))
    ! In normalized units
    gamma = sqrt(q_norm / (4.0 * PI * eps_r))

    ! Barrier lowering: Δφ_b = γ*sqrt(|E|)
    delta_phi_b = gamma * sqrt(E_smooth)

    ! Derivative: d(Δφ_b)/dE = γ * 0.5 * E / (|E| * sqrt(|E|))
    ! Using smoothed version to avoid division by zero
    if (E_smooth > eps_smooth) then
      d_delta_phi_dE = gamma * 0.5 * E_field / (E_smooth * sqrt(E_smooth))
    else
      d_delta_phi_dE = 0.0
    end if

  end subroutine

  pure subroutine calc_wkb_transmission(E, phi_b, F, m_tn, T_wkb, dT_dE, dT_dphi, dT_dF)
    !! Calculate WKB transmission probability through triangular Schottky barrier
    !! Implements: T_wkb = exp[-(4/3) * sqrt(2*m_tn) * (phi_b - E)^(3/2) / |F|]
    !!
    !! Following ATLAS Tsu-Esaki model (Eq. 3-174) with triangular barrier approximation
    !! Coefficient (4/3) assumes h-bar normalization (not h directly)
    !!
    !! All inputs/outputs in normalized units (kT/q normalization)

    real, intent(in)  :: E         !! Longitudinal kinetic energy (normalized)
    real, intent(in)  :: phi_b     !! Schottky barrier height (normalized)
    real, intent(in)  :: F         !! Electric field magnitude (normalized)
    real, intent(in)  :: m_tn      !! Tunneling effective mass ratio (m*/m0)
    real, intent(out) :: T_wkb     !! Transmission probability [0,1]
    real, intent(out) :: dT_dE     !! Derivative dT/dE
    real, intent(out) :: dT_dphi   !! Derivative dT/d(phi_b)
    real, intent(out) :: dT_dF     !! Derivative dT/dF

    real :: gamma               ! Tunneling exponent
    real :: F_smooth           ! Smoothed field magnitude
    real :: eps_smooth         ! Smoothing parameter
    real :: delta_E            ! Barrier height relative to energy: (phi_b - E)
    real :: coeff              ! WKB coefficient: (4/3) * sqrt(2*m_tn*PI)
    real :: dgamma_dE, dgamma_dphi, dgamma_dF  ! Exponent derivatives

    ! Smoothing parameter to avoid F=0 singularity
    eps_smooth = 1e-10

    ! Case 1: Above barrier - classical transmission
    if (E >= phi_b) then
      T_wkb = 1.0
      dT_dE = 0.0
      dT_dphi = 0.0
      dT_dF = 0.0
      return
    end if

    ! Case 2: Zero or very small field - no tunneling possible
    if (abs(F) < eps_smooth) then
      T_wkb = 0.0
      dT_dE = 0.0
      dT_dphi = 0.0
      dT_dF = 0.0
      return
    end if

    ! Case 3: Under-barrier tunneling
    delta_E = phi_b - E
    F_smooth = sqrt(F**2 + eps_smooth**2)

    ! WKB coefficient with h-bar normalization
    coeff = (4.0 / 3.0) * sqrt(2.0 * m_tn)

    ! Tunneling exponent: gamma = coeff * (phi_b - E)^(3/2) / |F|
    gamma = coeff * delta_E**1.5 / F_smooth

    ! Transmission probability
    T_wkb = exp(-gamma)

    ! Derivatives of gamma
    ! dgamma/dE = -coeff * (3/2) * (phi_b - E)^(1/2) / |F|
    dgamma_dE = -coeff * 1.5 * sqrt(delta_E) / F_smooth

    ! dgamma/d(phi_b) = coeff * (3/2) * (phi_b - E)^(1/2) / |F|
    dgamma_dphi = coeff * 1.5 * sqrt(delta_E) / F_smooth

    ! dgamma/dF = -coeff * (phi_b - E)^(3/2) * F / (|F|^3)
    ! Using F_smooth for numerical stability
    dgamma_dF = -coeff * delta_E**1.5 * F / (F_smooth**3)

    ! Transmission derivatives using chain rule: dT/dx = T * (-dgamma/dx)
    dT_dE = T_wkb * (-dgamma_dE)
    dT_dphi = T_wkb * (-dgamma_dphi)
    dT_dF = T_wkb * (-dgamma_dF)

  end subroutine

  subroutine barrier_lowering_init(this, par)
    !! Initialize barrier_lowering variable
    class(barrier_lowering), intent(out) :: this
    type(device_params),     intent(in)  :: par
    call this%variable_init("delta_phi_b", "eV", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
  end subroutine

  subroutine schottky_injection_init(this, par, ci)
    !! Initialize schottky_injection variable for storing n0B
    class(schottky_injection), intent(out) :: this
    type(device_params),       intent(in)  :: par
      !! device parameters
    integer,                   intent(in)  :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base - only exists at vertices
    call this%variable_init(CR_NAME(ci)//"n0b", "cm^-3", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
    this%ci = ci

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

  subroutine calc_schottky_injection_init(this, par, ci, efield, n0b, delta_phi_b)
    !! Initialize equation for calculating n0B from electric field
    class(calc_schottky_injection), intent(out) :: this
    type(device_params), target,    intent(in)  :: par
      !! device parameters
    integer,                        intent(in)  :: ci
      !! carrier index
    type(electric_field), target,   intent(in)  :: efield(:)
      !! electric field components
    type(schottky_injection), target, intent(in) :: n0b
      !! injection density variable
    type(barrier_lowering), target, intent(in) :: delta_phi_b
      !! barrier lowering variable for output

    integer :: ict, i, iprov, iprov_delta, idx(par%g%idx_dim), dir, idep
    integer :: normal_dir
    logical :: has_schottky
    integer, allocatable :: iprov_array(:), idep_array(:,:)

    print "(A)", "calc_schottky_injection_init for "//CR_NAME(ci)

    ! Check if we have any Schottky contacts
    has_schottky = .false.
    do ict = 1, par%nct
      if (par%contacts(ict)%type == CT_SCHOTTKY) then
        has_schottky = .true.
        exit
      end if
    end do

    if (.not. has_schottky) then
      print *, "  No Schottky contacts found, skipping initialization"
      return
    end if

    ! init base equation
    call this%equation_init("calc_"//CR_NAME(ci)//"n0b")
    this%par => par
    this%ci = ci

    ! init stencils
    call this%st_dir%init(par%g)
    call this%st_em%init()

    ! Allocate and determine normal directions for each contact
    allocate(this%contact_normal(par%nct))
    this%contact_normal = 0  ! Initialize all to 0 (not Schottky or invalid)

    do ict = 1, par%nct
      if (par%contacts(ict)%type == CT_SCHOTTKY) then
        this%contact_normal(ict) = get_schottky_contact_normal_dir(par, ict)
        if (this%contact_normal(ict) > 0 .and. this%contact_normal(ict) <= par%g%dim) then
          print *, "  Schottky contact ", ict, " normal direction: ", this%contact_normal(ict)
        else
          print *, "  WARNING: Could not determine normal for Schottky contact ", ict
        end if
      end if
    end do

    ! Store reference to electric field components and delta_phi_b variable
    this%efield => efield
    this%delta_phi_b_var => delta_phi_b

    ! Create n0b selector for all Schottky contact vertices
    ! We'll include all contacts but only Schottky ones will have non-zero values
    call this%n0b%init(n0b, [(par%transport_vct(ict)%get_ptr(), ict = 1, par%nct)])

    ! Provide n0B at Schottky contact vertices (provide for all, but only Schottky will be set)
    iprov = this%provide(n0b, [(par%transport_vct(ict)%get_ptr(), ict = 1, par%nct)])

    ! Provide delta_phi_b at all transport vertices for output
    iprov_delta = this%provide(delta_phi_b, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])

    ! Set up dependencies and Jacobians for each E-field direction that's needed
    allocate(this%jaco_efield(par%g%dim))

    ! Create dependencies for each E-field component that Schottky contacts need
    do dir = 1, par%g%dim
      if (any(this%contact_normal == dir)) then
        ! This direction is used by at least one Schottky contact
        print *, "  Setting up E-field dependency for direction ", dir

        ! Create dependency for this E-field component at all contact vertices
        ! Only Schottky contacts with this normal will actually use it
        idep = this%depend(efield(dir), [(par%transport_vct(ict)%get_ptr(), &
                                         ict = 1, par%nct)])

        ! Create Jacobian for this E-field component
        this%jaco_efield(dir)%p => this%init_jaco(iprov, idep, &
          st = [(this%st_dir%get_ptr(), ict = 1, par%nct)], &
          const = .false.)
      end if
    end do

    ! finish initialization
    call this%init_final()
  end subroutine

  subroutine calc_schottky_injection_eval(this)
    !! Evaluate n0B = n0B_base * exp(±delta_phi_b(E))
    class(calc_schottky_injection), intent(inout) :: this

    integer :: ict, i, j, idx(this%par%g%idx_dim), normal_dir
    integer :: edge_idx(this%par%g%idx_dim)
    logical :: edge_status
    real :: E_field, n0b_base, n0b_val, delta_phi_b, d_delta_phi_dE
    real :: eps_r
    real, allocatable :: tmp(:)

    ! Check if properly initialized
    if (.not. allocated(this%contact_normal)) return
    if (.not. allocated(this%jaco_efield)) return
    if (.not. associated(this%efield)) return

    allocate(tmp(this%n0b%n))
    tmp = 0.0

    ! Clear all Jacobians first
    do normal_dir = 1, this%par%g%dim
      if (associated(this%jaco_efield(normal_dir)%p)) then
        call this%jaco_efield(normal_dir)%p%reset()
      end if
    end do

    ! Counter for position in n0b array
    j = 0

    ! Loop over all contacts (n0b includes all, but we only set Schottky ones)
    do ict = 1, this%par%nct

      ! Process each vertex of this contact
      do i = 1, this%par%transport_vct(ict)%n
        idx = this%par%transport_vct(ict)%get_idx(i)
        j = j + 1  ! Increment position in n0b array

        ! Check if this is a Schottky contact with valid normal
        if (this%par%contacts(ict)%type /= CT_SCHOTTKY) then
          tmp(j) = 0.0  ! Non-Schottky contacts have zero injection
          cycle
        end if

        normal_dir = this%contact_normal(ict)
        if (normal_dir <= 0 .or. normal_dir > this%par%g%dim) then
          tmp(j) = 0.0  ! Invalid normal direction
          cycle
        end if

        ! Get permittivity from the edge connected to this vertex in the normal direction
        ! First find the edge connected to this vertex in the normal direction
        call this%par%g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, normal_dir, idx, 1, edge_idx, edge_status)

        if (edge_status) then
          ! Successfully found the edge, get permittivity from it
          eps_r = this%par%eps(IDX_EDGE, normal_dir)%get(edge_idx)
        else
          ! Fallback: use material default (should ideally get from semiconductor parameters)
          ! For now, use Silicon default
          eps_r = 11.7
          print "(A,I2,A)", "WARNING: Could not find edge for vertex at contact ", ict, &
                           ", using default eps_r = 11.7"
        end if

        ! Get the E-field component in the normal direction for this contact
        E_field = this%efield(normal_dir)%get(idx)


        ! Get base injection density (without field)
        call schottky_injection_mb(this%par, this%ci, ict, n0b_base)

        ! Calculate barrier lowering
        call schottky_barrier_lowering(this%par, ict, E_field, eps_r, delta_phi_b, d_delta_phi_dE)

        ! Store delta_phi_b in output variable
        if (associated(this%delta_phi_b_var)) then
          call this%delta_phi_b_var%set(idx, delta_phi_b)
        end if

        ! Apply barrier lowering to get n0B
        if (this%ci == CR_ELEC) then
          n0b_val = n0b_base * exp(delta_phi_b)
          ! Set Jacobian entry for the correct E-field component
          if (associated(this%jaco_efield(normal_dir)%p)) then
            call this%jaco_efield(normal_dir)%p%set(idx, idx, n0b_val * d_delta_phi_dE)
          end if
        else  ! CR_HOLE
          n0b_val = n0b_base * exp(-delta_phi_b)
          ! Set Jacobian entry for the correct E-field component
          if (associated(this%jaco_efield(normal_dir)%p)) then
            call this%jaco_efield(normal_dir)%p%set(idx, idx, -n0b_val * d_delta_phi_dE)
          end if
        end if

        ! Store in temporary array at position j
        tmp(j) = n0b_val
      end do
    end do

    ! Materialize all Jacobians
    do normal_dir = 1, this%par%g%dim
      if (associated(this%jaco_efield(normal_dir)%p)) then
        call this%jaco_efield(normal_dir)%p%set_matr(const = .false., nonconst = .true.)
      end if
    end do

    ! Set the n0B values
    call this%n0b%set(tmp)

    deallocate(tmp)
  end subroutine

end module
