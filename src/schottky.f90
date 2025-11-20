module schottky_m

  use device_params_m,  only: device_params
  use potential_m,      only: potential
  use density_m,        only: density
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
  public :: schottky_injection, schottky_tunnel_current, calc_schottky_injection
  public :: barrier_lowering
  public :: test_tsuesaki_vs_te_comparison

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

  type, extends(variable_real) :: schottky_tunnel_current
    !! Tunneling current density J_tn at Schottky contacts
    !! Only exists at Schottky contact vertices with tunneling enabled

    integer :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    real, pointer :: x1(:)     => null()
      !! direct pointer to data for easy access (only used if idx_dim == 1)
    real, pointer :: x2(:,:)   => null()
      !! direct pointer to data for easy access (only used if idx_dim == 2)
    real, pointer :: x3(:,:,:) => null()
      !! direct pointer to data for easy access (only used if idx_dim == 3)
  contains
    procedure :: init => schottky_tunnel_current_init
  end type

  type, extends(equation) :: calc_schottky_injection
    !! Calculate bias-dependent injection density n0B from electric field
    !! n0B = n0B_base * exp(delta_phi_b(E)) for electrons
    !! n0B = n0B_base * exp(-delta_phi_b(E)) for holes

    type(device_params), pointer :: par => null()
      !! pointer to device parameters
    integer :: ci
      !! carrier index

    type(potential), pointer :: pot => null()
      !! pointer to electrostatic potential variable
    type(electric_field), pointer :: efield(:) => null()
      !! electric field components (array of size dim)
    type(vselector) :: n0b
      !! injection density variable
    type(vselector) :: jtn
      !! tunneling current density variable (J_tn at contact vertices)
    type(density), pointer :: dens => null()
      !! carrier density variable (for diagnostics/validation)

    integer, allocatable :: contact_normal(:)
      !! normal direction for each contact (0 if not Schottky)

    type(dirichlet_stencil) :: st_dir
      !! stencil for contact vertices
    type(empty_stencil) :: st_em
      !! empty stencil

    type(jacobian_ptr), allocatable :: jaco_efield_n0b(:)
      !! jacobian for n0B E-field dependencies: ∂n0B/∂E (IFBL only, one per direction)
    type(jacobian_ptr), allocatable :: jaco_efield_jtn(:)
      !! jacobian for J_tn E-field dependencies: ∂J_tn/∂E (tunneling, one per direction)
    type(jacobian_ptr) :: jaco_pot_jtn
      !! jacobian for J_tn potential coupling: ∂J_tn/∂phi_k (tunneling)

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
    real :: A_rich

    ! Check if Richardson constant is provided and > 0

    if (ci == CR_ELEC) then
      A_rich = par%contacts(ict)%A_richardson_n
    else
      A_rich = par%contacts(ict)%A_richardson_p
    end if

    if (A_rich > 0.0) then
      ! Calculate normalized surface velocity
      ! v_surf = A*T^2/(q*Nc) where q is handled by normalization
      ! T must be normalized, Nc is already normalized
      s = A_rich * norm(par%T, "K") * norm(par%T, "K") / par%smc%edos(ci)

      if (first_call_v) then
        print "(A,A,A,I2,A,A,A)", "DEBUG_SCHOTTKY: Contact ", trim(par%contacts(ict)%name), &
                              " (", ict, ") ", CR_NAME(ci), " surface velocity:"
        print "(A,ES14.6,A)", "  A_richardson  = ", A_rich, " A/cm^2/K^2"
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

  pure function log1p_exp(x) result(y)
    !! Numerically stable computation of log(1 + exp(x))
    !!
    !! For large positive x:  log(1 + exp(x)) ≈ x
    !! For large negative x:  log(1 + exp(x)) ≈ exp(x)
    !! For moderate x:        Use direct formula
    !!
    !! This avoids overflow when exp(x) → ∞ and maintains accuracy

    real, intent(in) :: x
    real :: y

    if (x > 30.0) then
      ! For x > 30, exp(x) >> 1, so log(1 + exp(x)) ≈ log(exp(x)) = x
      y = x
    else if (x < -30.0) then
      ! For x < -30, exp(x) << 1, so log(1 + exp(x)) ≈ exp(x)
      y = exp(x)
    else
      ! For moderate x, use direct computation
      y = log(1.0 + exp(x))
    end if

  end function log1p_exp

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

    ! WKB coefficient for triangular barrier
    !
    ! Standard formula: γ = (4/3) * (√(2m*)/ℏ) * (φ - E)^(3/2) / |F|
    !
    ! In normalized units (kT/q for energy, device-specific for field):
    !   coeff = (4/3) * √(2 m*/m₀)
    !
    ! where m_tn = m*/m₀ is the tunneling effective mass ratio.
    !
    ! Note: Full normalization includes additional factors from the normalization
    ! scheme (e.g., √(m₀ kT / ℏ²) * ε₀εᵣ / q) which should be absorbed into the
    ! electric field normalization. Verify against experimental data or TCAD tools.
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

  subroutine tsu_esaki_integrand(E, p, f, dfdE, dfdp)
    !! Integrand for Tsu-Esaki tunneling current integration
    !! Implements: f(E) = T_wkb(E) * N(E)
    !! where N(E) = ln[(1+exp((E_F-E)/kT)) / (1+exp((E_F-E-qV)/kT))]
    !!
    !! This function matches the interface required by quad_m
    !! All outputs are REQUIRED (cannot be omitted)
    !!
    !! Parameter array structure:
    !!   p(1) = phi_eff: Effective barrier height (phi_bn for e-, E_g-phi_bn for h+)
    !!   p(2) = efield : Electric field magnitude
    !!   p(3) = m_tn   : Tunneling effective mass ratio
    !!   p(4) = phi_k  : Electrostatic potential at contact interface (from Poisson)
    !!   p(5) = ci     : Carrier index (1=electron, 2=hole)
    !!   p(6) = E_g    : Band gap (normalized)

    real, intent(in)  :: E       !! Energy (integration variable, normalized)
    real, intent(in)  :: p(:)    !! Parameter array [phi_eff, efield, m_tn, phi_k, ci, E_g]
    real, intent(out) :: f       !! Integrand value
    real, intent(out) :: dfdE    !! Derivative wrt energy
    real, intent(out) :: dfdp(:) !! Derivatives wrt parameters

    real :: phi_eff, efield, m_tn, phi_k, E_g
    integer :: ci
    real :: T_wkb, dT_dE, dT_dphi, dT_dF
    real :: N_E, dN_dE, dN_dphi  ! Occupancy difference (logarithmic form) and derivatives
    real :: f_s, f_m             ! Fermi-Dirac functions (for derivative calculation)
    real :: phi_semi             ! Semiconductor potential for occupancy

    ! Extract parameters from array
    phi_eff = p(1)
    efield  = p(2)
    m_tn    = p(3)
    phi_k   = p(4)
    ci      = int(p(5))
    E_g     = p(6)


    call calc_wkb_transmission(E, phi_eff, efield, m_tn, T_wkb, dT_dE, dT_dphi, dT_dF)

    ! Calculate occupancy difference using logarithmic form
    ! Energy reference: Metal Fermi level (E_F = 0)
    !
    ! Electrons (CB): phi_semi = phi_k (CB edge relative to metal EF)
    ! Holes (VB):     phi_semi = -E_g + phi_k (VB edge relative to metal EF)
    !
    ! N(E) for electrons: log(1 + exp(-E)) - log(1 + exp(-E - phi_semi))
    ! N(E) for holes:     log(1 + exp(E)) - log(1 + exp(E + phi_semi))

    if (ci == CR_ELEC) then
      ! ELECTRONS
      ! Electron occupancy difference
      N_E = log1p_exp(-E) - log1p_exp(-E - phi_k)

      ! Fermi functions for electrons
      f_s = 1.0 / (1.0 + exp(E))              ! f(E) in metal
      f_m = 1.0 / (1.0 + exp(E + phi_k))   ! f(E) in semiconductor

      ! Derivatives for electrons
      dN_dE = f_m - f_s
      dN_dphi = f_m

    else  ! CR_HOLE
      ! HOLES

      ! Hole occupancy difference
      N_E = log1p_exp(E) - log1p_exp(E + phi_k)

      ! Fermi functions for holes (note different signs)
      f_s = 1.0 / (1.0 + exp(-E))             ! 1-f(E) in metal = hole occupancy
      f_m = 1.0 / (1.0 + exp(-E - phi_k))  ! 1-f(E) in semiconductor = hole occupancy
      ! Derivatives for holes
      dN_dE = f_s - f_m     ! Opposite sign from electrons
      dN_dphi = -f_m        ! Opposite sign from electrons

    end if

    ! Integrand: f = T_wkb * N(E)
    f = T_wkb * N_E

    ! Derivative wrt energy (chain rule)
    ! df/dE = dT/dE * N + T_wkb * dN/dE
    dfdE = dT_dE * N_E + T_wkb * dN_dE

    ! Derivatives wrt parameters
    ! df/d(phi_eff) = dT/d(phi_eff) * N
    dfdp(1) = dT_dphi * N_E

    ! df/dF = dT/dF * N
    dfdp(2) = dT_dF * N_E

    ! df/d(m_tn): Set to zero (WKB doesn't return dT/dm)
    dfdp(3) = 0.0

    ! df/d(phi_k) = T_wkb * dN/d(phi_k)
    dfdp(4) = T_wkb * dN_dphi

    ! df/d(ci): Set to zero (carrier index is discrete)
    dfdp(5) = 0.0

    ! df/d(E_g): Set to zero (not computing this derivative)
    dfdp(6) = 0.0

  end subroutine

  subroutine tsu_esaki_current(phi_b, efield, m_tn, phi_k, ci, E_g, J_tn, dJ_tn_dF, dJ_tn_dpot)
    !! Calculate Tsu-Esaki tunneling current density through Schottky barrier
    !! Implements: J_tn = ±(m*/(2π²)) ∫₀^φ_eff T_wkb(E) N(E) dE
    !!
    !! Sign convention:
    !!   Electrons (ci=CR_ELEC): Negative prefactor (injection current)
    !!   Holes (ci=CR_HOLE):     Positive prefactor (extraction current)
    !!
    !! Barrier heights:
    !!   Electrons: phi_eff = phi_b (electron barrier to CB)
    !!   Holes:     phi_eff = E_g - phi_b (hole barrier to VB)
    !!
    !! This is the tunneling component only (under-barrier transport)
    !! Total current: J_total = J_TE + J_tn (ATLAS split model)
    !!
    !! All inputs/outputs in normalized units (kT/q normalization)

    real, intent(in)  :: phi_b     !! Schottky barrier height for electrons (with IFBL if applied)
    real, intent(in)  :: efield    !! Electric field magnitude (normalized)
    real, intent(in)  :: m_tn      !! Tunneling effective mass ratio (m*/m0)
    real, intent(in)  :: phi_k     !! Electrostatic potential at contact (from Poisson)
    integer, intent(in) :: ci      !! Carrier index (CR_ELEC or CR_HOLE)
    real, intent(in)  :: E_g       !! Band gap (normalized)
    real, intent(out) :: J_tn      !! Tunneling current density (normalized)
    real, intent(out) :: dJ_tn_dF  !! Derivative dJ_tn/dF for Jacobian
    real, intent(out), optional :: dJ_tn_dpot  !! Derivative dJ_tn/d(phi_k) for potential coupling

    real :: p(6)                   ! Parameter array for integrand
    real :: I                      ! Integration result
    real :: dIda, dIdb            ! Derivatives wrt bounds (not used)
    real :: dIdp(6)               ! Derivatives wrt parameters
    real :: err                   ! Integration error estimate
    integer :: ncalls             ! Number of function evaluations
    real :: prefactor             ! Normalization prefactor: ±m*/(2π²)
    real :: phi_eff               ! Effective barrier for this carrier

    ! Calculate effective barrier based on carrier type
    if (ci == CR_ELEC) then
      phi_eff = phi_b              ! Electrons: phi_bn (CB barrier)
    else  ! CR_HOLE
      phi_eff = E_g - phi_b        ! Holes: phi_bp = E_g - phi_bn (VB barrier)
    end if

    ! Check for degenerate cases
    if (phi_eff <= 0.0) then
      ! No barrier - no tunneling calculation needed
      J_tn = 0.0
      dJ_tn_dF = 0.0
      if (present(dJ_tn_dpot)) dJ_tn_dpot = 0.0
      return
    end if

    if (abs(efield) < 1e-10) then
      ! Zero field - no tunneling possible
      J_tn = 0.0
      dJ_tn_dF = 0.0
      if (present(dJ_tn_dpot)) dJ_tn_dpot = 0.0
      return
    end if

    ! Set up parameter array for integrand
    ! p(1) = phi_eff : effective barrier (carrier-dependent)
    ! p(2) = efield  : electric field
    ! p(3) = m_tn    : tunneling mass
    ! p(4) = phi_k   : contact potential
    ! p(5) = ci      : carrier index
    ! p(6) = E_g     : band gap
    p(1) = phi_eff
    p(2) = efield
    p(3) = m_tn
    p(4) = phi_k
    p(5) = real(ci)
    p(6) = E_g

    ! Perform integration (under-barrier only)
    ! Using adaptive quadrature with relative tolerance 1e-6
    !
      call quad(tsu_esaki_integrand, 0.0, phi_eff, p, I, dIda, dIdb, dIdp, &
                rtol=1e-6, err=err, max_levels=12, ncalls=ncalls)

    ! Apply normalization prefactor with carrier-dependent sign
    ! Electrons: negative (injection)
    ! Holes: positive (extraction)
    if (ci == CR_ELEC) then
      prefactor = m_tn / (2.0 * PI**2)  ! Negative for electron injection
    else  ! CR_HOLE
      prefactor = -m_tn / (2.0 * PI**2)  ! Positive for hole extraction
    end if

    ! Calculate tunneling current and field derivative
    J_tn = prefactor * I
    dJ_tn_dF = prefactor * dIdp(2)  ! dI/d(efield)

    ! Calculate potential derivative if requested
    if (present(dJ_tn_dpot)) then
      dJ_tn_dpot = prefactor * dIdp(4)  ! dI/d(phi_k)
    end if

  end subroutine

  subroutine test_tsuesaki_vs_te_comparison(phi_b, efield, m_tn, phi_k, ci, E_g, v_surf, n0b_base, delta_phi_b, n_actual)
    !! DEVICE VALIDATION: Tsu-Esaki Integral vs Robin BC Thermionic Current
    !!
    !! Test if Tsu-Esaki integral from phi_b → ∞ equals Robin BC thermionic current
    !! using the actual solved carrier density from the device.
    !!
    !! Physics hypothesis:
    !! Since T_wkb(E) → 1 for E > phi_b (over-barrier regime),
    !! the integral ∫_{phi_b}^∞ T_wkb(E) * N(E) dE
    !! should match J_TE = v_surf × (n - n0B_IFBL)
    !!
    !! This tests self-consistency between:
    !! - Tsu-Esaki integral formulation (with Fermi-Dirac occupancy N(E))
    !! - Robin BC formulation (with actual device carrier density)

    real, intent(in) :: phi_b         !! Schottky barrier height for electrons (normalized)
    real, intent(in) :: efield        !! Electric field magnitude (normalized)
    real, intent(in) :: m_tn          !! Tunneling effective mass ratio
    real, intent(in) :: phi_k         !! Contact potential (normalized)
    integer, intent(in) :: ci         !! Carrier index (CR_ELEC or CR_HOLE)
    real, intent(in) :: E_g           !! Band gap (normalized)
    real, intent(in) :: v_surf        !! Surface velocity (normalized)
    real, intent(in) :: n0b_base      !! Base injection density (normalized)
    real, intent(in) :: delta_phi_b   !! Barrier lowering from IFBL (normalized)
    real, intent(in) :: n_actual      !! Actual carrier density at contact (normalized)

    real :: p(6)                   ! Parameter array for integrand
    real :: I_overbarrier          ! Integration result for over-barrier regime
    real :: dIda, dIdb            ! Derivatives wrt bounds (not used)
    real :: dIdp(6)               ! Derivatives wrt parameters (not used)
    real :: err                   ! Integration error estimate
    integer :: ncalls             ! Number of function evaluations
    real :: prefactor             ! Normalization prefactor: ±m*/(2π²)
    real :: J_tn_overbarrier      ! Tsu-Esaki current in over-barrier regime
    real :: J_TE_net              ! Net thermionic emission current from Robin BC
    real :: n0B_IFBL              ! Field-enhanced injection density
    real :: upper_limit           ! Upper integration limit
    real :: relative_diff         ! Relative difference between methods
    real :: phi_bn_eff            ! Effective barrier
    real :: phi_eff               ! Carrier-dependent effective barrier

    print *
    print "(A)", "========================================================"
    print "(A)", "DEVICE VALIDATION"
    print "(A)", "Tsu-Esaki Integral vs Robin BC TE Current"
    print "(A,A,A)", "Carrier: ", CR_NAME(ci), ""
    print "(A)", "========================================================"
    print *

    ! Calculate effective barrier and field-enhanced n0B
    phi_bn_eff = phi_b - delta_phi_b

    ! Calculate carrier-dependent effective barrier
    if (ci == CR_ELEC) then
      phi_eff = phi_bn_eff              ! Electrons: phi_bn
      n0B_IFBL = n0b_base * exp(delta_phi_b)
    else  ! CR_HOLE
      phi_eff = E_g - phi_bn_eff        ! Holes: phi_bp = E_g - phi_bn
      n0B_IFBL = n0b_base * exp(-delta_phi_b)
    end if

    ! Print input parameters with carrier-aware labeling
    print "(A)", "Input Parameters:"
    print "(A,ES14.6,A)", "  phi_bn (e- barrier) = ", denorm(phi_bn_eff, "eV"), " eV"
    if (ci == CR_HOLE) then
      print "(A,ES14.6,A)", "  phi_bp (h+ barrier) = ", denorm(phi_eff, "eV"), " eV"
      print "(A,ES14.6,A)", "  E_g (band gap)      = ", denorm(E_g, "eV"), " eV"
    end if
    print "(A,ES14.6,A)", "  E_field             = ", denorm(efield, "V/cm"), " V/cm"
    print "(A,ES14.6)",   "  m_tn                = ", m_tn
    print "(A,ES14.6,A)", "  phi_k (potential)   = ", denorm(phi_k, "V"), " V"
    print "(A,ES14.6,A)", "  v_surf              = ", denorm(v_surf, "cm/s"), " cm/s"
    if (ci == CR_ELEC) then
      print "(A,ES14.6,A)", "  n (actual)          = ", denorm(n_actual, "cm^-3"), " cm^-3"
      print "(A,ES14.6,A)", "  n0B_base            = ", denorm(n0b_base, "cm^-3"), " cm^-3"
      print "(A,ES14.6,A)", "  n0B_IFBL            = ", denorm(n0B_IFBL, "cm^-3"), " cm^-3"
    else
      print "(A,ES14.6,A)", "  p (actual)          = ", denorm(n_actual, "cm^-3"), " cm^-3"
      print "(A,ES14.6,A)", "  p0B_base            = ", denorm(n0b_base, "cm^-3"), " cm^-3"
      print "(A,ES14.6,A)", "  p0B_IFBL            = ", denorm(n0B_IFBL, "cm^-3"), " cm^-3"
    end if
    print "(A,ES14.6,A)", "  delta_phi_b (IFBL)  = ", denorm(delta_phi_b, "eV"), " eV"
    print *

    ! Method 1: Robin BC thermionic current (using actual device carrier density)
    ! Electrons: J_TE_net = v_surf × (n - n0B_IFBL)
    ! Holes:     J_TE_net = v_surf × (p - p0B_IFBL)
    J_TE_net = v_surf * (n_actual - n0B_IFBL)

    print "(A)", "METHOD 1: Robin BC Thermionic Current"
    if (ci == CR_ELEC) then
      print "(A)", "  J_TE_net = v_surf × (n - n0B_IFBL)"
    else
      print "(A)", "  J_TE_net = v_surf × (p - p0B_IFBL)"
    end if
    print "(A,ES14.6,A,ES14.6)", "    v_surf              = ", v_surf, " (norm), ", denorm(v_surf, "cm/s"), " cm/s"
    print "(A,ES14.6,A,ES14.6)", "    carrier_actual      = ", n_actual, " (norm), ", denorm(n_actual, "cm^-3"), " cm^-3"
    print "(A,ES14.6,A,ES14.6)", "    carrier_0B_IFBL     = ", n0B_IFBL, " (norm), ", denorm(n0B_IFBL, "cm^-3"), " cm^-3"
    print "(A,ES14.6,A,ES14.6)", "    J_TE_net            = ", J_TE_net, " (norm), ", denorm(J_TE_net, "A/cm^2"), " A/cm^2"
    print *

    ! Method 2: Tsu-Esaki integral from phi_eff → ∞ (over-barrier regime)
    ! For E > phi_eff, T_wkb = 1.0, so integral is ∫_{phi_eff}^∞ N(E) dE

    ! Choose upper limit: large enough to capture all current
    ! Use 20*kT or phi_eff + 20, whichever is larger
    upper_limit = max(phi_eff + 20.0, 20.0)

    ! Set up parameter array (now 6 elements for new integrand)
    p(1) = phi_eff     ! effective barrier (carrier-dependent)
    p(2) = efield      ! electric field
    p(3) = m_tn        ! tunneling mass
    p(4) = phi_k       ! contact potential
    p(5) = real(ci)    ! carrier index
    p(6) = E_g         ! band gap

    ! Perform integration from phi_eff to upper_limit (over-barrier regime)
    ! In this regime, T_wkb(E) → 1.0 for E > phi_eff
    call quad(tsu_esaki_integrand, phi_eff, upper_limit, p, I_overbarrier, &
              dIda, dIdb, dIdp, rtol=1e-6, err=err, max_levels=12, ncalls=ncalls)

    ! Apply normalization prefactor with carrier-dependent sign
    if (ci == CR_ELEC) then
      prefactor = +m_tn / (2.0 * PI**2)  ! Negative for electrons
    else
      prefactor = -m_tn / (2.0 * PI**2)  ! Positive for holes
    end if
    J_tn_overbarrier = prefactor * I_overbarrier

    print "(A)", "METHOD 2: Tsu-Esaki Integral (Over-Barrier)"
    print "(A)", "  J = prefactor × ∫_{phi_eff}^∞ T_wkb(E) * N(E) dE"
    print "(A)", "  For E > phi_eff: T_wkb ≈ 1.0 (perfect transmission)"
    print "(A,ES14.6,A,ES14.6,A)", "  Integration limits  = ", denorm(phi_eff, "eV"), " to ", denorm(upper_limit, "eV"), " eV"
    print "(A,I8)",       "  Function calls      = ", ncalls
    print "(A,ES14.6)",   "  Integration error   = ", err
    print *
    print "(A)", "  Integration breakdown:"
    print "(A,ES14.6)",   "    Prefactor (m*/(2π²)) = ", prefactor
    print "(A,ES14.6)",   "    Prefactor (denorm)    = ", denorm(prefactor, "C*s/kg^-1/m^4"), " C·s/kg·m⁴"
    print "(A,ES14.6)",   "    Supply function I    = ", I_overbarrier
    print "(A,ES14.6,A)", "    J = prefactor × I    = ", denorm(J_tn_overbarrier, "A/cm^2"), " A/cm^2"
    print *

    ! Compare the two methods
    if (abs(J_TE_net) > 1e-30) then
      relative_diff = abs(J_tn_overbarrier - J_TE_net) / abs(J_TE_net)
      print "(A)", "COMPARISON:"
      print "(A,ES14.6,A)", "  J_Robin_BC (Method 1)   = ", denorm(J_TE_net, "A/cm^2"), " A/cm^2"
      print "(A,ES14.6,A)", "  J_integral (Method 2)   = ", denorm(J_tn_overbarrier, "A/cm^2"), " A/cm^2"
      print "(A,ES14.6)",   "  Relative difference     = ", relative_diff * 100.0, " %"
      print *

      if (relative_diff < 0.01) then
        print "(A)", "  ✓ EXCELLENT AGREEMENT (< 1%)"
      else if (relative_diff < 0.1) then
        print "(A)", "  ✓ GOOD AGREEMENT (< 10%)"
      else if (relative_diff < 0.5) then
        print "(A)", "  ~ MODERATE AGREEMENT (< 50%)"
      else
        print "(A)", "  ✗ POOR AGREEMENT (> 50%)"
        print "(A)", "  → Models may be inconsistent or integration limits insufficient"
      end if
    else
      print "(A)", "COMPARISON:"
      print "(A)", "  J_TE_net too small for meaningful comparison"
    end if

    print *
    print "(A)", "========================================================"
    print *

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

  subroutine schottky_tunnel_current_init(this, par, ci)
    !! Initialize schottky_tunnel_current variable for storing J_tn
    class(schottky_tunnel_current), intent(out) :: this
    type(device_params),            intent(in)  :: par
      !! device parameters
    integer,                        intent(in)  :: ci
      !! carrier index (CR_ELEC, CR_HOLE)

    type(grid_data1_real), pointer :: p1
    type(grid_data2_real), pointer :: p2
    type(grid_data3_real), pointer :: p3

    ! init base - only exists at vertices, units are A/cm^2
    call this%variable_init(CR_NAME(ci)//"jtn", "A/cm^2", g = par%g, idx_type = IDX_VERTEX, idx_dir = 0)
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

  subroutine calc_schottky_injection_init(this, par, ci, efield, n0b, delta_phi_b, pot, jtn_current, dens)
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
    type(potential), target, intent(in) :: pot
      !! electrostatic potential variable
    type(schottky_tunnel_current), target, intent(in) :: jtn_current
      !! tunneling current density variable
    type(density), target, intent(in) :: dens
      !! carrier density variable (for diagnostics)

    integer :: ict, i, iprov, iprov_delta, iprov_jtn, idx(par%g%idx_dim), dir, idep
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

    ! Store reference to electric field components, delta_phi_b variable, potential, and density
    this%efield => efield
    this%delta_phi_b_var => delta_phi_b
    this%pot => pot
    this%dens => dens

    ! Create n0b selector for all Schottky contact vertices
    ! We'll include all contacts but only Schottky ones will have non-zero values
    call this%n0b%init(n0b, [(par%transport_vct(ict)%get_ptr(), ict = 1, par%nct)])

    ! Create jtn selector for all Schottky contact vertices with tunneling enabled
    call this%jtn%init(jtn_current, [(par%transport_vct(ict)%get_ptr(), ict = 1, par%nct)])

    ! Provide n0B at Schottky contact vertices (provide for all, but only Schottky will be set)
    iprov = this%provide(n0b, [(par%transport_vct(ict)%get_ptr(), ict = 1, par%nct)])

    ! Provide J_tn at Schottky contact vertices with tunneling enabled
    iprov_jtn = this%provide(jtn_current, [(par%transport_vct(ict)%get_ptr(), ict = 1, par%nct)])

    ! Provide delta_phi_b at all transport vertices for output (only for electrons to avoid duplicate)
    ! Barrier lowering is carrier-independent, so we only need to compute it once
    if (this%ci == CR_ELEC) then
      iprov_delta = this%provide(delta_phi_b, [(par%transport_vct(ict)%get_ptr(), ict = 0, par%nct)])
    end if

    ! Set up dependencies and Jacobians for each E-field direction
    ! Separate Jacobians for n0B (thermionic/IFBL) and J_tn (tunneling)
    allocate(this%jaco_efield_n0b(par%g%dim))
    allocate(this%jaco_efield_jtn(par%g%dim))

    ! Create dependencies for each E-field component that Schottky contacts need
    do dir = 1, par%g%dim
      if (any(this%contact_normal == dir)) then
        ! This direction is used by at least one Schottky contact
        print *, "  Setting up E-field dependency for direction ", dir

        ! Create dependency for this E-field component at all contact vertices
        idep = this%depend(efield(dir), [(par%transport_vct(ict)%get_ptr(), &
                                         ict = 1, par%nct)])

        ! Create Jacobian for n0B: ∂n0B/∂E (IFBL only, attached to iprov)
        this%jaco_efield_n0b(dir)%p => this%init_jaco(iprov, idep, &
          st = [(this%st_dir%get_ptr(), ict = 1, par%nct)], &
          const = .false.)

        ! Create Jacobian for J_tn: ∂J_tn/∂E (tunneling, attached to iprov_jtn)
        this%jaco_efield_jtn(dir)%p => this%init_jaco(iprov_jtn, idep, &
          st = [(this%st_dir%get_ptr(), ict = 1, par%nct)], &
          const = .false.)
      end if
    end do

    ! Set up potential dependency and Jacobian for J_tn coupling
    ! ∂J_tn/∂phi_k (tunneling only, attached to iprov_jtn)
    print *, "  Setting up potential dependency for tunneling coupling"

    ! Create dependency on potential at all contact vertices
    idep = this%depend(pot, [(par%transport_vct(ict)%get_ptr(), ict = 1, par%nct)])

    ! Create Jacobian for J_tn: ∂J_tn/∂phi_k (attached to iprov_jtn)
    this%jaco_pot_jtn%p => this%init_jaco(iprov_jtn, idep, &
      st = [(this%st_dir%get_ptr(), ict = 1, par%nct)], &
      const = .false.)

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
    real :: phi_k, phi_bn_eff, v_surf
    real :: J_TE, J_tn, J_total, dJ_tn_dF, dJ_tn_dpot
    real :: n0B_IFBL, dn0B_eff_dE, dn0B_eff_dpot
    real :: n_actual_val
    real, allocatable :: tmp_n0b(:), tmp_jtn(:)
    real :: m_tun , m_tun_debug, m_tun_test


    ! Check if properly initialized
    if (.not. allocated(this%contact_normal)) return
    if (.not. allocated(this%jaco_efield_n0b)) return
    if (.not. associated(this%efield)) return

    allocate(tmp_n0b(this%n0b%n))
    allocate(tmp_jtn(this%jtn%n))
    tmp_n0b = 0.0
    tmp_jtn = 0.0

    ! Clear all Jacobians first
    do normal_dir = 1, this%par%g%dim
      ! Reset n0B Jacobians (thermionic/IFBL)
      if (associated(this%jaco_efield_n0b(normal_dir)%p)) then
        call this%jaco_efield_n0b(normal_dir)%p%reset()
      end if
      ! Reset J_tn Jacobians (tunneling)
      if (associated(this%jaco_efield_jtn(normal_dir)%p)) then
        call this%jaco_efield_jtn(normal_dir)%p%reset()
      end if
    end do

    ! Reset J_tn potential Jacobian
    if (associated(this%jaco_pot_jtn%p)) then
      call this%jaco_pot_jtn%p%reset()
    end if

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
          tmp_n0b(j) = 0.0  ! Non-Schottky contacts have zero injection
          tmp_jtn(j) = 0.0
          cycle
        end if

        normal_dir = this%contact_normal(ict)
        if (normal_dir <= 0 .or. normal_dir > this%par%g%dim) then
          tmp_n0b(j) = 0.0  ! Invalid normal direction
          tmp_jtn(j) = 0.0
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

        ! Calculate barrier lowering (IFBL)
        call schottky_barrier_lowering(this%par, ict, E_field, eps_r, delta_phi_b, d_delta_phi_dE)

        ! Store delta_phi_b in output variable (only electrons provide this variable)
        if (this%ci == CR_ELEC .and. associated(this%delta_phi_b_var)) then
          call this%delta_phi_b_var%set(idx, delta_phi_b)
        end if

        ! Get effective barrier height (including IFBL if enabled)
        phi_bn_eff = this%par%contacts(ict)%phi_b - delta_phi_b

        ! Get contact potential from Poisson solution
        if (associated(this%pot)) then
          phi_k = this%pot%get(idx)
          ! Debug: print potential information
          if (i == 1 .and. this%ci == CR_ELEC) then
            print "(A,I2,A,A,A)", "  Contact ", ict, " (", trim(this%par%contacts(ict)%name), ") potential fetch:"
            print "(A,3I5)", "    idx = ", idx
            print "(A,ES14.6,A)", "    phi_k (Poisson) = ", phi_k, " (normalized)"
            print "(A,ES14.6,A)", "    phi_k (denorm)  = ", denorm(phi_k, "V"), " V"
            print "(A,ES14.6,A)", "    phims (stored)  = ", this%par%contacts(ict)%phims, " (normalized)"
            print "(A,ES14.6,A)", "    phims (denorm)  = ", denorm(this%par%contacts(ict)%phims, "V"), " V"
            print "(A,ES14.6,A)", "    phi_b (barrier) = ", denorm(this%par%contacts(ict)%phi_b, "eV"), " eV"
            print "(A)", "    Expected: phi_k should equal phims at equilibrium"
          end if
        else
          phi_k = 0.0  ! Fallback if potential not available
          if (i == 1 .and. this%ci == CR_ELEC) then
            print "(A)", "  WARNING: Potential not associated, using phi_k = 0"
          end if
        end if

        ! Calculate surface velocity (thermionic emission velocity)
        v_surf = schottky_velocity(this%par, this%ci, ict)

        ! Calculate thermionic emission current (over-barrier transport)
        if (this%ci == CR_ELEC) then
          J_TE = n0b_base * exp(delta_phi_b) * v_surf
        else  ! CR_HOLE
          J_TE = n0b_base * exp(-delta_phi_b) * v_surf
        end if

        ! Calculate tunneling current if enabled (under-barrier transport)
        if (this%par%contacts(ict)%tunneling) then

          if (this%ci == CR_ELEC) then
            m_tun = this%par%contacts(ict)%m_tunnel_n
          else
            m_tun = this%par%contacts(ict)%m_tunnel_p
          end if
          call tsu_esaki_current(phi_bn_eff, E_field, m_tun, &
                                 phi_k, this%ci, this%par%smc%band_gap, &
                                 J_tn, dJ_tn_dF, dJ_tn_dpot)
        else
          J_tn = 0.0
          dJ_tn_dF = 0.0
          dJ_tn_dpot = 0.0
        end if

        ! NO FOLDING: Keep n0B as thermionic emission only (n0B_IFBL)
        ! Tunneling J_tn will be added as explicit boundary source in continuity equation
        ! This keeps n0B always positive and makes physics transparent

        if (this%ci == CR_ELEC) then
          ! Electrons: n0B_IFBL = n0B * exp(delta_phi_b) from IFBL
          n0B_IFBL = n0b_base * exp(delta_phi_b)
          n0b_val = n0B_IFBL  ! No folding - just thermionic

          ! Jacobian: d(n0B_IFBL)/dE (thermionic only, no tunneling term here)
          dn0B_eff_dE = n0B_IFBL * d_delta_phi_dE

        else  ! CR_HOLE
          ! Holes: n0B_IFBL = n0B * exp(-delta_phi_b)
          n0B_IFBL = n0b_base * exp(-delta_phi_b)
          n0b_val = n0B_IFBL  ! No folding

          ! Jacobian: d(n0B_IFBL)/dE (negative sign for holes)
          dn0B_eff_dE = -n0B_IFBL * d_delta_phi_dE
        end if

        ! Debug output for tunneling current analysis
        if (this%par%contacts(ict)%tunneling .and. mod(i, 10) == 1) then
          ! Enhanced debug for carrier-specific barriers
          if (i == 1) then

            if (this%ci == CR_ELEC) then
              m_tun_debug = this%par%contacts(ict)%m_tunnel_n
            else
              m_tun_debug = this%par%contacts(ict)%m_tunnel_p
            end if
            print "(A,I2,A,A,A,A,A)", "Contact ", ict, " (", trim(this%par%contacts(ict)%name), &
                                      ") ", CR_NAME(this%ci), " tunneling:"
            print "(A,F6.3,A)", "  m_tunnel = ", m_tun_debug, " m0 (carrier-specific)"
            print "(A,ES14.6,A)", "  phi_bn (e- barrier) = ", denorm(phi_bn_eff, "eV"), " eV"
            if (this%ci == CR_HOLE) then
              print "(A,ES14.6,A)", "  phi_bp (h+ barrier) = ", &
                    denorm(this%par%smc%band_gap - phi_bn_eff, "eV"), " eV"
              print "(A,ES14.6,A)", "  E_g (band gap)      = ", denorm(this%par%smc%band_gap, "eV"), " eV"
            end if
            print "(A,ES14.6,A)", "  E_field             = ", denorm(E_field, "V/cm"), " V/cm"
            print "(A,ES14.6,A)", "  J_tn                = ", denorm(J_tn, "A/cm^2"), " A/cm^2"
            print "(A,ES14.6,A)", "  J_TE                = ", denorm(J_TE, "A/cm^2"), " A/cm^2"
            print "(A,ES14.6)",   "  J_tn/J_TE ratio     = ", abs(J_tn / (J_TE + 1e-30))
            print "(A,A)",         "  Current sign        = ", merge("NEGATIVE (injection) ", &
                                                                     "POSITIVE (extraction)", this%ci == CR_ELEC)
          end if

          ! Run comparison test (only for first vertex of first contact)
          if (i == 1 .and. ict == 1) then
            ! Fetch actual carrier density at this vertex (same way as pot)
            if (associated(this%dens)) then
              n_actual_val = this%dens%get(idx)
            else
              n_actual_val = 0.0
            end if
            if (this%ci == CR_ELEC) then
              m_tun_test = this%par%contacts(ict)%m_tunnel_n
            else
              m_tun_test = this%par%contacts(ict)%m_tunnel_p
            end if
            call test_tsuesaki_vs_te_comparison(this%par%contacts(ict)%phi_b, E_field, &
                 m_tun_test, phi_k, this%ci, this%par%smc%band_gap, &
                 v_surf, n0b_base, delta_phi_b, n_actual_val)
          end if
        end if

        ! Set Jacobian entry for n0B: ∂n0B_IFBL/∂E (thermionic/IFBL only)
        if (associated(this%jaco_efield_n0b(normal_dir)%p)) then
          call this%jaco_efield_n0b(normal_dir)%p%set(idx, idx, dn0B_eff_dE)
        end if

        ! Set Jacobian entries for tunneling current J_tn
        ! These are attached to iprov_jtn, so chain rule works correctly:
        ! ∂F/∂E += ∂F/∂J_tn × ∂J_tn/∂E  and  ∂F/∂phi += ∂F/∂J_tn × ∂J_tn/∂phi
        ! Works for both electrons and holes
        if (this%par%contacts(ict)%tunneling) then
          ! ∂J_tn/∂E for field coupling
          if (associated(this%jaco_efield_jtn(normal_dir)%p)) then
            call this%jaco_efield_jtn(normal_dir)%p%set(idx, idx, dJ_tn_dF)
          end if
          ! ∂J_tn/∂phi_k for potential coupling
          if (associated(this%jaco_pot_jtn%p)) then
            call this%jaco_pot_jtn%p%set(idx, idx, dJ_tn_dpot)
          end if
        end if

        ! Store in temporary arrays at position j
        tmp_n0b(j) = n0b_val  ! Thermionic only
        tmp_jtn(j) = J_tn      ! Store tunneling current for continuity to use
      end do
    end do

    ! Materialize all Jacobians
    do normal_dir = 1, this%par%g%dim
      ! Materialize n0B Jacobians (thermionic/IFBL)
      if (associated(this%jaco_efield_n0b(normal_dir)%p)) then
        call this%jaco_efield_n0b(normal_dir)%p%set_matr(const = .false., nonconst = .true.)
      end if
      ! Materialize J_tn Jacobians (tunneling)
      if (associated(this%jaco_efield_jtn(normal_dir)%p)) then
        call this%jaco_efield_jtn(normal_dir)%p%set_matr(const = .false., nonconst = .true.)
      end if
    end do

    ! Materialize J_tn potential Jacobian
    if (associated(this%jaco_pot_jtn%p)) then
      call this%jaco_pot_jtn%p%set_matr(const = .false., nonconst = .true.)
    end if

    ! Set the n0B values (thermionic only)
    call this%n0b%set(tmp_n0b)

    ! Set the J_tn values (tunneling current for continuity to use)
    call this%jtn%set(tmp_jtn)

    deallocate(tmp_n0b, tmp_jtn)
  end subroutine

end module
