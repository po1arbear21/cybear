# Schottky Contact Implementation Notes
## Cybear General-Purpose Drift-Diffusion Simulator

**Author**: Chenyang Yu
**Date**: 2025-10-06
**Simulator**: Cybear DD - Supporting multiple device architectures

---

> **Note**: This document covers Schottky contact implementation for the Cybear general-purpose DD simulator. While examples include various device types (nanowires, 2D materials, perovskites), the physics and implementation are universally applicable to all semiconductor devices with metal-semiconductor junctions.

## Table of Contents
1. [Current Implementation Status](#current-implementation-status)
2. [Material Parameter Handling Issues](#material-parameter-handling-issues)
3. [Tunneling Physics for Ultrathin Contacts](#tunneling-physics-for-ultrathin-contacts)
4. [Implementation Roadmap](#implementation-roadmap)
5. [Device-Specific Tunneling Considerations](#device-specific-tunneling-considerations)
6. [Validation Strategy](#validation-strategy)
7. [Historical Development](#historical-development)

---

## Current Implementation Status

### ✅ Implemented Features

#### 1. **Thermionic Emission (TE) Model**
- Location: `src/schottky.f90`
- Robin boundary conditions in continuity equation
- Richardson constant support (A* = 112 A/cm²/K² for Si)
- Surface recombination velocity: v_surf = A*T²/(q*Nc)
- Zero-bias injection: n₀ = Nc × exp(-φ_Bn/kT)

#### 2. **Image Force Barrier Lowering (IFBL)**
- Schottky effect: Δφ_b = √(q|E|/(4π*ε))
- Field-dependent injection: n₀B = n₀ × exp(±Δφ_b/kT)
- Automatic contact normal direction detection
- Smoothing parameter ε = 1e-10 to avoid E=0 singularity
- Can be enabled/disabled per contact via `ifbl` flag

#### 3. **Electric Field Coupling**
- E-field calculated at vertices from potential gradients (`src/electric_field.f90`)
- Field components passed to Schottky injection calculation
- Full Jacobian matrix for E-field dependencies
- Proper coupling with continuity equation via `calc_schottky_injection`

#### 4. **Contact Infrastructure**
- CT_SCHOTTKY contact type (value = 3)
- Support for barrier height (phi_b) and Richardson constant (A_richardson)
- Proper parsing from configuration files
- Integration with device parameter system

### ⚠️ Issues Fixed Today (2025-10-06)

1. **Permittivity Retrieval** (FIXED)
   - Problem: Used vertex index to access edge data
   - Solution: Added `get_neighb` call to find connected edge
   - Location: `calc_schottky_injection_eval` line 450-463
   ```fortran
   ! NEW (correct):
   call this%par%g%get_neighb(IDX_VERTEX, 0, IDX_EDGE, normal_dir, idx, 1, edge_idx, edge_status)
   if (edge_status) then
     eps_r = this%par%eps(IDX_EDGE, normal_dir)%get(edge_idx)
   else
     eps_r = 11.7  ! Fallback to Si default
   end if
   ```

2. **Hardcoded Material Parameters** (FIXED)
   - Problem: eps_r = 11.7 hardcoded in `schottky_injection_mb_bias`
   - Solution: Added eps_r as function parameter
   - Note: Function marked as legacy (not used in main flow)

---

## Material Parameter Handling Issues

### Fixed Issues
- Permittivity now properly retrieved from device parameters via edge connectivity
- Richardson constant properly normalized with temperature scaling
- Fallback values for edge cases (though should ideally come from semiconductor parameters)

### Remaining Issues
1. Fallback permittivity should come from semiconductor parameters, not hardcoded
2. Missing tunneling effective mass in device parameters
3. Need interface doping concentration for TFE calculations

---

## Tunneling Physics for Ultrathin Contacts

### 🚨 **CRITICAL**: Current Implementation Lacks Tunneling

The current pure thermionic emission model is **insufficient** for modern nanoscale devices:

#### Universal Criteria (Any Device Type)
- Barrier thickness < 10 nm
- Doping concentration > 10¹⁸ cm⁻³
- Electric fields > 10⁶ V/cm
- Low temperatures < 200 K

#### Device-Specific Thresholds
| Device Type | Critical Dimension | When Tunneling Dominates |
|-------------|-------------------|--------------------------|
| **Si Nanowire FET** | d < 10nm | Always (quantum confinement) |
| **2D Material FET** | t < 1nm | Always (atomically thin) |
| **FinFET** | Fin width < 7nm | High doping or low T |
| **GAA FET** | d < 5nm | Always |
| **Perovskite VFET** | Channel < 50nm | Perforated contacts |
| **Schottky Diode** | Any | N_D > 10¹⁸ cm⁻³ |

### Required Transport Mechanisms

#### 1. **Thermionic-Field Emission (TFE)**

**Characteristic tunneling energy**:
```
E₀₀ = (qℏ/2) × √(N_D/(m*×ε_s))
```

**Regime determination**:
- kT/E₀₀ >> 1: Thermionic emission dominates
- kT/E₀₀ ≈ 1: Thermionic-field emission
- kT/E₀₀ << 1: Field emission dominates

**Padovani-Stratton TFE model**:
```
J_TFE = A*T² × exp(-qφ_eff/kT)
φ_eff = φ_b - ξ

where:
ξ = E₀₀ × coth(E₀₀/kT)  (energy position of tunneling maximum)
E₀ = E₀₀ × coth(E₀₀/kT)  (characteristic energy)
```

#### 2. **WKB Transmission Probability**

**General WKB formula**:
```
T(E) = exp(-2∫√(2m*(V(x)-E))/ℏ dx)
```

**Triangular barrier (Fowler-Nordheim)**:
```
T_FN = exp(-4√(2m*)φ_b^(3/2)/(3qℏE))
```

**Trapezoidal barrier (direct tunneling)**:
```
T_DT = exp(-4√(2m*q)φ_b × t_ox/(3ℏ(1+V/φ_b)))
```

#### 3. **Field Emission (FE)**

**Pure field emission current (T→0)**:
```
J_FE = (A*q³E²/8πhφ_b) × exp(-4√(2m*)φ_b^(3/2)/(3qℏE))
```

#### 4. **Direct Tunneling (DT)**

**For ultrathin barriers (< 3 nm)**:
```
J_DT = (q³m₀V²/8πℏ²m*d³) × exp(-4√(2m*q)φ_b × d/(3ℏ))
```

#### 5. **Tsu-Esaki Unified Tunneling Model**

**Comprehensive tunneling formulation with WKB transmission**:

The Tsu-Esaki model provides a unified framework for all tunneling regimes (TE/TFE/FE) through a single integral:

```
J = (4πqm*kT/h³) ∫₀^φb T(E) [ln(1+exp((EF-E)/kT)) - ln(1+exp((EF-qV-E)/kT))] dE
```

**Key features**:
- Integration limits from **0 to φ_b** (not to infinity)
- T(E) is the WKB transmission coefficient
- Naturally interpolates between all transport regimes
- Temperature-independent transmission probability

**Triangular barrier WKB transmission**:
```
T(E) = exp[-4√(2m*)(φb-E)^(3/2)/(3ℏqE)]
```

**Advantages over split TE/TFE models**:
1. **Physically consistent** - Single framework for all transport
2. **No artificial transitions** - Smooth interpolation between regimes
3. **Accurate at all biases** - Valid from equilibrium to high field
4. **Better for research** - More rigorous for publications

**Disadvantages**:
1. **Computationally intensive** - Requires numerical integration
2. **Convergence challenges** - Nested integrals (WKB + energy)
3. **Harder to debug** - Single integral masks individual physics

**Reference implementation (MATLAB)**:
```matlab
function J = tsu_esaki_J(E_V_per_cm, phi_b_eV, mstar, T, V_contact)
    % Physical constants
    q = 1.602e-19; kB = 1.381e-23; h = 6.626e-34; hbar = h/(2*pi);

    % Energy grid (0 to phi_b)
    Ez_eV = linspace(0, phi_b_eV, 1200);
    Ez_J = Ez_eV * q;

    % Occupancy difference
    arg1 = (EF - Ez_J)/(kB*T);
    arg2 = (EF - q*V_contact - Ez_J)/(kB*T);
    Nlog = log1p(exp(arg1)) - log1p(exp(arg2));

    % WKB transmission
    coeff = 4*sqrt(2*mstar)/(3*q*hbar);
    T_wkb = exp(-coeff * (phi_b_J - Ez_J).^(3/2) / E_SI);

    % Integration
    pref = 4*pi*q*mstar*kB*T/h^3;
    J = pref * trapz(Ez_J, T_wkb .* Nlog);
end
```

---

## Implementation Roadmap

### Phase 1: Essential TFE Model (IMMEDIATE PRIORITY)

```fortran
module schottky_tunneling_m
  use normalization_m, only: norm, denorm
  use math_m, only: PI

  implicit none

  type schottky_tunneling_params
    real :: E00           ! Characteristic energy (normalized)
    real :: m_tunnel      ! Tunneling effective mass ratio
    real :: N_interface   ! Interface doping concentration
    logical :: enable_tfe ! Enable TFE model
    logical :: enable_dt  ! Enable direct tunneling
    logical :: enable_fe  ! Enable field emission
  end type

contains

  function calc_E00(N_D, eps_r, m_eff) result(E00)
    !! Calculate characteristic tunneling energy
    real, intent(in) :: N_D    ! Doping concentration (normalized)
    real, intent(in) :: eps_r  ! Relative permittivity
    real, intent(in) :: m_eff  ! Effective mass ratio
    real :: E00

    ! E00 = (q*hbar/2) * sqrt(N_D/(m_eff*eps_0*eps_r))
    ! In normalized units:
    E00 = 0.5 * sqrt(N_D/(m_eff*eps_r))
  end function

  function calc_tfe_enhancement(E00, T, phi_b, E_field) result(f_tfe)
    !! Calculate TFE enhancement factor over pure TE
    real, intent(in) :: E00, T, phi_b, E_field
    real :: f_tfe
    real :: E0, xi, kT_E00

    kT_E00 = T/E00  ! Temperature ratio

    if (kT_E00 > 10.0) then
      ! Pure thermionic regime
      f_tfe = 1.0
    elseif (kT_E00 < 0.1) then
      ! Pure field emission
      f_tfe = exp(4.0*sqrt(2.0*phi_b**3)/(3.0*sqrt(E_field)))
    else
      ! TFE regime
      E0 = E00/tanh(E00/T)
      xi = E00*coth(E00/T)
      f_tfe = (E0/T)/sin(PI*T/E0) * exp((phi_b - xi)/T)
    end if
  end function

  subroutine calc_schottky_current_total(params, phi_b, E_field, T, V_bias, &
                                         J_total, dJ_dV, dJ_dE)
    !! Calculate total Schottky current including all mechanisms
    type(schottky_tunneling_params), intent(in) :: params
    real, intent(in)  :: phi_b, E_field, T, V_bias
    real, intent(out) :: J_total, dJ_dV, dJ_dE

    real :: J_TE, J_TFE, J_DT, J_FE
    real :: f_tfe, T_wkb

    ! Base thermionic emission
    J_TE = calc_thermionic_current(phi_b, T)

    ! TFE enhancement
    if (params%enable_tfe .and. params%E00 > 0.01*T) then
      f_tfe = calc_tfe_enhancement(params%E00, T, phi_b, E_field)
      J_TFE = J_TE * f_tfe
    else
      J_TFE = J_TE
    end if

    ! Direct tunneling (ultrathin barriers)
    if (params%enable_dt .and. phi_b < 0.1) then  ! < 100 meV
      J_DT = calc_direct_tunneling(phi_b, E_field, V_bias)
    else
      J_DT = 0.0
    end if

    ! Field emission (high field, low temp)
    if (params%enable_fe .and. E_field > 1e6 .and. T < 100) then
      J_FE = calc_field_emission(phi_b, E_field, params%m_tunnel)
    else
      J_FE = 0.0
    end if

    ! Total current
    J_total = J_TFE + J_DT + J_FE

    ! Derivatives for Newton solver
    call calc_current_derivatives(J_total, V_bias, E_field, dJ_dV, dJ_dE)
  end subroutine

end module
```

### Phase 2: WKB Implementation

```fortran
function calc_wkb_transmission(E, phi_b, E_field, d_barrier) result(T_wkb)
  !! WKB transmission through arbitrary barrier
  real, intent(in) :: E         ! Carrier energy
  real, intent(in) :: phi_b     ! Barrier height
  real, intent(in) :: E_field   ! Electric field
  real, intent(in) :: d_barrier ! Barrier width
  real :: T_wkb

  real :: x_tp  ! Classical turning point
  real :: gamma ! Tunneling exponent

  ! Find turning point where E = V(x)
  x_tp = min((phi_b - E)/(q*E_field), d_barrier)

  if (E > phi_b) then
    ! Over-barrier transport
    T_wkb = 1.0
  elseif (x_tp < d_barrier) then
    ! Triangular barrier
    gamma = 4.0*sqrt(2.0*m_eff)*(phi_b - E)**1.5/(3.0*hbar*q*E_field)
    T_wkb = exp(-gamma)
  else
    ! Trapezoidal barrier
    gamma = 4.0*sqrt(2.0*m_eff*q)/(3.0*hbar) * &
            (phi_b - E)**1.5 * (1.0 - (1.0 - q*E_field*d_barrier/(phi_b-E))**1.5) / &
            (q*E_field)
    T_wkb = exp(-gamma)
  end if
end function
```

### Phase 3: Integration with Continuity Equation

```fortran
! Modify continuity.f90 boundary conditions
if (par%contacts(ict)%type == CT_SCHOTTKY) then
  ! Get tunneling parameters
  call get_tunneling_params(par, ict, st_params)

  ! Calculate total current including tunneling
  call calc_schottky_current_total(st_params, phi_b, E_field, T, V_bias, &
                                   J_total, dJ_dV, dJ_dE)

  ! Update Robin BC coefficients
  v_surf_eff = J_total/(n - n0B)  ! Effective surface velocity
  call this%jaco_dens%set(idx1, idx1, A_ct * v_surf_eff)

  ! Add field derivative
  if (associated(this%jaco_efield)) then
    call this%jaco_efield%set(idx1, idx1, A_ct * dJ_dE)
  end if
end if
```

### Phase 4: Configuration Parameters

```ini
[schottky parameters]
  # Material parameters
  m_tunnel    = 0.3      : m0      # Tunneling effective mass
  m_dos       = 0.5      : m0      # DOS effective mass

  # Model selection
  enable_tfe  = true               # Thermionic-field emission
  enable_dt   = true               # Direct tunneling
  enable_fe   = true               # Field emission
  enable_btbt = false              # Band-to-band tunneling

  # Interface parameters
  N_interface = 1e19     : 1/cm^3  # Doping at interface
  d_barrier   = 2        : nm      # Effective barrier width
  trap_density = 1e12    : 1/cm^2  # Interface trap density

  # Numerical parameters
  E_mesh      = 100                # Energy grid points
  smooth_factor = 1e-10            # Field smoothing parameter
```

### Phase 5: Tsu-Esaki Tunneling Implementation

```fortran
module schottky_tunneling_tsu_esaki_m
  use device_params_m
  use normalization_m
  use quad_m  ! For tanh-sinh quadrature

  implicit none

contains

  function tsu_esaki_current_norm(E_field_norm, phi_b_norm, m_star, V_norm) result(J_norm)
    !! Normalized Tsu-Esaki current density
    !! All inputs/outputs in normalized units (kT/q normalization)

    real(dp), intent(in) :: E_field_norm  ! Normalized electric field
    real(dp), intent(in) :: phi_b_norm    ! Normalized barrier height
    real(dp), intent(in) :: m_star        ! Effective mass ratio (m*/m0)
    real(dp), intent(in) :: V_norm        ! Normalized voltage drop
    real(dp) :: J_norm                    ! Normalized current density

    ! Perform integration using adaptive quadrature
    call integrate_tsu_esaki_tanh_sinh(E_field_norm, phi_b_norm, m_star, V_norm, J_norm)

  end function

  subroutine integrate_tsu_esaki_tanh_sinh(E_field, phi_b, m_star, V, result)
    !! Core integration using tanh-sinh quadrature (recommended)
    use quad_m, only: quad_tanhsinh

    real(dp), intent(in) :: E_field, phi_b, m_star, V
    real(dp), intent(out) :: result
    real(dp) :: prefactor

    ! Normalization-aware prefactor
    ! In normalized units: J = (prefactor) * integral
    prefactor = 4.0_dp * PI * m_star  ! Simplified in normalized units

    ! Use tanh-sinh for smooth exponential integrands
    call quad_tanhsinh(integrand_func, 0.0_dp, phi_b, result, &
                      rtol=1e-6_dp, atol=1e-10_dp)

    result = prefactor * result

  contains

    function integrand_func(E) result(f)
      real(dp), intent(in) :: E
      real(dp) :: f
      real(dp) :: T_wkb, N_diff
      real(dp) :: E_smooth, coeff

      ! Smooth E_field to avoid division by zero
      E_smooth = sqrt(E_field**2 + 1e-10_dp)

      if (E < phi_b .and. E_smooth > 1e-20_dp) then
        ! Triangular barrier WKB transmission
        ! T = exp[-4√(2m*)/(3ℏq) * (φb-E)^(3/2) / E]
        ! In normalized units with proper scaling
        coeff = (4.0_dp/3.0_dp) * sqrt(2.0_dp * m_star * PI)
        T_wkb = exp(-coeff * (phi_b - E)**(1.5_dp) / E_smooth)
      else
        T_wkb = 1.0_dp  ! Above barrier or E→0 limit
      end if

      ! Occupancy difference using log1p for numerical stability
      ! N(E) = ln(1 + exp((EF-E)/kT)) - ln(1 + exp((EF-qV-E)/kT))
      ! In normalized units (kT=1):
      N_diff = log1p(exp(-E)) - log1p(exp(-E - V))

      f = T_wkb * N_diff

    end function

  end subroutine

  ! Alternative: Gauss-Legendre quadrature for comparison
  subroutine integrate_tsu_esaki_gauss(E_field, phi_b, m_star, V, result)
    !! Integration using Gauss-Legendre quadrature
    use gauss_m, only: gauss_legendre_nodes_weights

    real(dp), intent(in) :: E_field, phi_b, m_star, V
    real(dp), intent(out) :: result
    real(dp), allocatable :: nodes(:), weights(:)
    real(dp) :: E, T_wkb, N_diff, sum_integral
    integer :: n_points, i

    n_points = 100  ! Can be adaptive
    allocate(nodes(n_points), weights(n_points))

    ! Get Gauss-Legendre nodes and weights for [0, phi_b]
    call gauss_legendre_nodes_weights(n_points, 0.0_dp, phi_b, nodes, weights)

    sum_integral = 0.0_dp
    do i = 1, n_points
      E = nodes(i)
      ! Calculate T_wkb and N_diff (same as above)
      ! ... (code omitted for brevity)
      sum_integral = sum_integral + weights(i) * T_wkb * N_diff
    end do

    result = 4.0_dp * PI * m_star * sum_integral

    deallocate(nodes, weights)
  end subroutine

end module
```

#### Integration with Existing Schottky Module

```fortran
! Modified schottky.f90
subroutine schottky_injection_with_tunneling(par, ci, ict, E_field, ninj, dninj_dE)
  use schottky_tunneling_tsu_esaki_m

  type(device_params), intent(in) :: par
  integer, intent(in) :: ci, ict
  real, intent(in) :: E_field  ! Field magnitude at contact
  real, intent(out) :: ninj
  real, intent(out), optional :: dninj_dE

  real :: J_tunnel, v_surf, delta_J
  real :: E_perturb

  if (par%contacts(ict)%tunneling) then
    ! Tsu-Esaki tunneling current
    J_tunnel = tsu_esaki_current_norm( &
      E_field, &
      par%contacts(ict)%phi_b, &
      par%contacts(ict)%m_tunnel, &
      0.0_dp)  ! V_contact from BC

    ! Convert to injection density: n = J/(q*v)
    v_surf = calculate_surface_velocity(par, ict)
    ninj = J_tunnel / v_surf

    ! Calculate derivative for Jacobian (finite difference)
    if (present(dninj_dE)) then
      E_perturb = E_field * 1.001_dp
      delta_J = tsu_esaki_current_norm( &
        E_perturb, par%contacts(ict)%phi_b, &
        par%contacts(ict)%m_tunnel, 0.0_dp) - J_tunnel
      dninj_dE = delta_J / (v_surf * 0.001_dp * E_field)
    end if

  else
    ! Pure thermionic emission (existing code)
    call schottky_injection_mb(par, ci, ict, ninj)

    ! Apply IFBL if enabled
    if (par%contacts(ict)%ifbl) then
      call schottky_barrier_lowering(par, ict, E_field, eps_r, &
                                    delta_phi_b, d_delta_phi_dE)
      ninj = ninj * exp(delta_phi_b)
      if (present(dninj_dE)) then
        dninj_dE = ninj * d_delta_phi_dE
      end if
    end if
  end if

end subroutine
```

#### Configuration Parameters

```fortran
! In contact.f90, extend contact type
type contact
  ! ... existing fields ...

  ! Tunneling parameters
  logical :: tunneling = .false.     ! Enable Tsu-Esaki tunneling
  real :: m_tunnel = 1.0              ! Tunneling effective mass ratio (m*/m0)

  ! Can combine with existing IFBL
  logical :: ifbl = .false.          ! Image force barrier lowering
end type
```

```ini
[contact]
  name            = "SCHOTTKY"
  type            = "schottky"
  phi_b           = 0.7           : eV
  A_richardson    = 112           : A/cm^2/K^2

  # Tunneling configuration
  tunneling       = true           : Enable Tsu-Esaki tunneling
  m_tunnel        = 0.4            : Tunneling effective mass ratio

  # Can combine with IFBL (applied to TE component)
  ifbl            = false          : Image force barrier lowering
```

#### Critical Implementation Notes

1. **Normalization**: Cybear uses kT/q energy normalization
   - 1 energy unit = kT/q (thermal voltage)
   - Electric field normalized by L₀/V_T
   - Current density needs proper denormalization factor

2. **Integration Method**: Tanh-sinh recommended because:
   - Handles exponential decay at E→φb smoothly
   - Adaptive refinement for varying scales
   - Works well with log1p occupancy terms

3. **Numerical Stability**:
   - Use log1p for occupancy to avoid overflow
   - Smooth E_field with small epsilon (1e-10)
   - Cache integration results when E_field unchanged

4. **Performance Optimization**:
```fortran
type :: tsu_esaki_cache
  real :: last_E_field = -1.0
  real :: last_phi_b = -1.0
  real :: last_J = 0.0
  logical :: valid = .false.
end type

! Check cache before recalculating
if (abs(E_field - cache%last_E_field)/E_field < 1e-4) then
  J = cache%last_J
else
  J = tsu_esaki_current_norm(...)
  cache%last_E_field = E_field
  cache%last_J = J
end if
```

#### ⚠️ **CRITICAL: Current Density Sign Convention**

**Atlas Manual Convention**: The Atlas manual uses a **negative sign** for electron tunneling current density J_tn.

When applying boundary conditions in Cybear:
- If the interface normal **n̂ points OUT of the semiconductor**
- And the BC expects **current LEAVING the semiconductor as positive**
- Then use: **`J_tn_signed = -J_tn_magnitude`**

```fortran
! Sign convention check
normal_dir = get_schottky_contact_normal_dir(par, ict)
if (contact_normal_points_out) then
  ! Electrons flowing out → negative current by convention
  J_boundary = -J_tunnel
else
  J_boundary = J_tunnel
end if
```

This sign convention is critical for:
- Proper current continuity at interfaces
- Correct I-V characteristics (forward vs reverse bias)
- Convergence of the Newton solver

---

## Device-Specific Tunneling Considerations

### Applications Across Different Device Types

#### 1. **Nanowire FETs (Sub-10nm diameter)**
```
Quantum confinement effects:
- Subband formation modifies barrier shape
- Surface states dominate at Si/SiO₂ interface
- Radial tunneling through gate oxide
- Axial tunneling at source/drain Schottky contacts

Critical parameters:
- Wire diameter < 10nm → strong confinement
- Surface-to-volume ratio → interface trap dominance
- Gate-all-around → radial field enhancement
```

#### 2. **2D Material FETs (MoS₂, WSe₂, Graphene)**
```
Van der Waals contacts:
- No Fermi level pinning → tunable Schottky barrier
- Interlayer tunneling in vdW heterojunctions
- Edge contact vs surface contact geometries
- Thickness-dependent bandgap (1L vs bulk)

Tunneling considerations:
- Atomically thin channels (~0.7nm) → direct S-D tunneling
- Low DOS → enhanced TFE at moderate doping
- Anisotropic effective mass → direction-dependent tunneling
```

#### 3. **Perovskite Vertical FETs (Project A07)**
```
Perforated Source Electrode (80nm apertures):
```
Field enhancement at aperture edges:
β = 2-5 (geometric factor)
E_local = β × E_applied

Requires:
- Field emission at aperture edges
- 3D field distribution modeling
- Local barrier height variation
```

#### 2. **Ultra-short Channel (10-50nm)**
```
Source-drain tunneling:
- Becomes significant < 20nm
- Direct S-D tunneling in OFF state
- Band-to-band tunneling (BTBT) needed
```

#### 3. **Ion Migration Effects**
```
Time-dependent barrier:
φ_b(t) = φ_b0 + Δφ_ion(t)

where:
Δφ_ion = q × N_ion × x_ion / ε

Requires:
- Coupled ion-electron transport
- Dynamic barrier recalculation
- Hysteresis modeling
```

#### 4. **Mixed Conduction Mechanisms**
```
Total current:
J_total = J_drift + J_thermionic + J_tunneling + J_ionic

With transitions:
- Low field: Ohmic (J ∝ V)
- Medium field: SCLC (J ∝ V²)
- High field: Tunneling (J ∝ exp(-1/E))
```

### Implementation Priority for Advanced Devices

#### Universal Requirements (All Nanoscale Devices)
1. **MUST HAVE** (for accuracy):
   - TFE model with E₀₀ calculation
   - WKB transmission for arbitrary barriers
   - Field-dependent barrier lowering

2. **SHOULD HAVE** (for completeness):
   - Direct tunneling for ultrathin barriers
   - Interface trap states
   - Temperature-dependent effective mass

3. **DEVICE-SPECIFIC**:
   - **Nanowires**: Quantum confinement, subband structure
   - **2D Materials**: Interlayer tunneling, anisotropic transport
   - **Perovskites**: Ion migration, dynamic barriers
   - **Silicon**: Band-to-band tunneling, hot carriers

---

## Validation Strategy

### 1. **Unit Tests**
```fortran
! Test E00 calculation
N_D = norm(1e18, "1/cm^3")
E00_expected = norm(4.5e-3, "eV")  ! For Si at 1e18
E00_calc = calc_E00(N_D, 11.7, 0.26)
assert(abs(E00_calc - E00_expected) < 0.001)

! Test regime transitions
assert(is_thermionic(kT=0.026, E00=0.001))     ! kT/E00 = 26
assert(is_tfe(kT=0.026, E00=0.026))            ! kT/E00 = 1
assert(is_field_emission(kT=0.026, E00=0.26))  ! kT/E00 = 0.1
```

### 2. **I-V Characteristic Tests**
```python
# Expected behavior with tunneling
def test_iv_with_tunneling():
    # Low doping (1e16): pure thermionic
    J_1e16 = simulate(N_D=1e16, enable_tfe=False)
    J_1e16_tfe = simulate(N_D=1e16, enable_tfe=True)
    assert(J_1e16_tfe / J_1e16 < 1.1)  # < 10% difference

    # High doping (1e19): significant tunneling
    J_1e19 = simulate(N_D=1e19, enable_tfe=False)
    J_1e19_tfe = simulate(N_D=1e19, enable_tfe=True)
    assert(J_1e19_tfe / J_1e19 > 10)  # > 10x enhancement
```

### 3. **Temperature Dependence**
```
Plot: ln(J/T²) vs 1/T

Expected:
- Pure TE: Linear (slope = -qφ_b/k)
- With TFE: Curved (reduced slope at low T)
- Pure FE: Temperature independent
```

### 4. **Field Dependence**
```
Plot: ln(J/E²) vs 1/E (Fowler-Nordheim plot)

Expected:
- High field: Linear region (slope ∝ φ_b^(3/2))
- Low field: Deviation from linearity
```

### 5. **Benchmark Cases**

| Device Type | Test Case | Parameters | Expected J (A/cm²) |
|-------------|-----------|------------|-------------------|
| **Silicon** | n-type TE | N_D=1e16, φ_b=0.7V, T=300K | ~10⁻⁶ |
| **Silicon** | n-type TFE | N_D=1e19, φ_b=0.7V, T=300K | ~10⁻² |
| **Silicon** | Field emission | N_D=1e19, φ_b=0.7V, T=77K, E=1e6 V/cm | ~10¹ |
| **Si Nanowire** | GAA FET | d=5nm, φ_b=0.5V, T=300K | ~10⁻⁴ |
| **MoS₂ FET** | Monolayer | N_D=1e12/cm², φ_b=0.2V, T=300K | ~10⁻⁵ |
| **WSe₂ FET** | Bilayer | N_D=1e13/cm², φ_b=0.3V, T=300K | ~10⁻⁴ |
| **Perovskite** | Vertical FET | N_D=1e18, φ_b=0.4V, T=300K | ~10⁻³ |

### 6. **Comparison with Sentaurus Device**
```
Run identical structure in Sentaurus with:
- Thermionic emission model
- Nonlocal tunneling model
- Compare I-V curves within 20%
```

### 7. **Tsu-Esaki Model Validation**

#### Unit Tests for Tsu-Esaki Implementation
```fortran
! Test 1: Zero field limit → pure thermionic emission
E_field = 1e-10  ! Essentially zero
J_tsu = tsu_esaki_current_norm(E_field, phi_b, m_star, V)
J_te = thermionic_current_norm(phi_b, V)
assert(abs(J_tsu - J_te)/J_te < 0.01)  ! < 1% difference

! Test 2: High field limit → Fowler-Nordheim
E_field = 1e8  ! Very high field (normalized)
J_tsu = tsu_esaki_current_norm(E_field, phi_b, m_star, V)
J_fn = fowler_nordheim_current(E_field, phi_b, m_star)
assert(abs(J_tsu - J_fn)/J_fn < 0.1)  ! < 10% difference

! Test 3: Integration limits (0 to φb)
! Verify energy grid spans correct range
assert(E_min == 0.0)
assert(E_max == phi_b)

! Test 4: Sign convention check
J_magnitude = abs(J_tsu)
J_signed = apply_sign_convention(J_magnitude, normal_dir)
assert(J_signed < 0 if electrons_leaving_semiconductor)
```

#### Comparison with MATLAB Reference
```matlab
% Test configuration
E_field = 1e5;      % V/cm
phi_b = 0.7;        % eV
m_star = 0.4;       % Effective mass ratio
T = 300;            % K

% Run both implementations
J_matlab = tsu_esaki_J(E_field, phi_b, m_star*9.1e-31, T, 0);
J_fortran = denorm(tsu_esaki_current_norm(...), "A/cm^2");

% Should match within numerical tolerance
relative_error = abs(J_matlab - J_fortran)/J_matlab;
assert(relative_error < 1e-3);  % 0.1% tolerance
```

#### Field-Dependent I-V Curves
```python
# Generate validation plots
fields = [1e4, 1e5, 1e6]  # V/cm
for E in fields:
    V = np.linspace(0, 1.0, 100)
    J_te = [thermionic_only(v, E) for v in V]
    J_tsu = [tsu_esaki(v, E) for v in V]

    plt.semilogy(V, J_te, '--', label=f'TE only, E={E:.0e}')
    plt.semilogy(V, J_tsu, '-', label=f'Tsu-Esaki, E={E:.0e}')

# Expected: Tsu-Esaki shows higher current at high fields
```

#### Temperature Dependence Validation
```
Test: Arrhenius plot ln(J/T²) vs 1000/T

Temperature range: 77K to 400K
Expected behavior:
- Pure TE: Linear with slope = -qφb/k
- Tsu-Esaki: Deviation from linearity at low T (tunneling contribution)
- At T→0: Tsu-Esaki approaches finite value (pure tunneling)
```

---

## Code Quality Checklist

- [ ] Remove all debug print statements
- [ ] Add proper error handling for edge cases
- [ ] Implement unit tests for each function
- [ ] Add convergence monitoring for iterative solvers
- [ ] Document all physical assumptions
- [ ] Validate normalization consistency
- [ ] Profile performance bottlenecks
- [ ] Add configuration validation

---

## Historical Development

### Previous Implementation Milestones

1. **Step 1: Contact Type Infrastructure** (2025-08-29)
   - Added CT_SCHOTTKY constant
   - Integrated phi_b and A_richardson parameters
   - Updated region parsing and device initialization

2. **Step 2: Robin BC Implementation** (2025-08-29)
   - Created schottky.f90 module
   - Implemented thermionic emission model
   - Fixed stencil architecture for mixed BC types
   - Achieved proper I-V characteristics

3. **Step 3: Electric Field Integration** (2025-09-01)
   - Added electric field calculation
   - Implemented image force barrier lowering
   - Created field-dependent injection model
   - Added calc_schottky_injection equation

4. **Step 4: Material Parameter Fixes** (2025-10-06)
   - Fixed permittivity retrieval from edges
   - Removed hardcoded material parameters
   - Improved error handling and fallbacks

5. **Step 5: Tsu-Esaki Tunneling Specification** (2025-10-30)
   - Comprehensive Tsu-Esaki unified model theory
   - Detailed implementation specification for Cybear
   - Integration with existing TE and IFBL infrastructure
   - Critical sign convention documentation
   - Comparison with split TE/TFE approaches
   - Added MATLAB reference implementation
   - Specified tanh-sinh quadrature integration
   - Performance optimization strategies

---

## Next Steps

1. **Immediate** (This Week):
   - [ ] Implement Tsu-Esaki tunneling module (schottky_tunneling_tsu_esaki_m)
   - [ ] Add m_tunnel and tunneling flag to contact type
   - [ ] Integrate with existing schottky.f90
   - [ ] Verify sign convention for current density
   - [ ] Create unit tests for Tsu-Esaki integration
   - [ ] Compare with MATLAB reference implementation

2. **Short Term** (Next 2 Weeks):
   - [ ] Full WKB implementation
   - [ ] Field emission for perforated contacts
   - [ ] Benchmark against Sentaurus Device
   - [ ] Performance optimization

3. **Long Term**:
   - [ ] Band-to-band tunneling
   - [ ] Trap-assisted tunneling
   - [ ] Hot carrier effects
   - [ ] Full perovskite physics

---

## References

1. **Padovani-Stratton TFE Model**:
   F.A. Padovani and R. Stratton, "Field and thermionic-field emission in Schottky barriers," Solid-State Electronics, vol. 9, pp. 695-707, 1966.

2. **WKB Approximation**:
   S.M. Sze and K.K. Ng, "Physics of Semiconductor Devices," 3rd ed., Wiley, 2007, Ch. 3.

3. **Fowler-Nordheim Tunneling**:
   R.H. Fowler and L. Nordheim, "Electron emission in intense electric fields," Proc. R. Soc. Lond. A, vol. 119, pp. 173-181, 1928.

4. **Perovskite FETs**:
   Project A07 internal documentation, RWTH Aachen University.

5. **Sentaurus Device**:
   Synopsys Sentaurus Device User Guide, Version O-2018.06.

---

**Last Updated**: 2025-10-30
**Status**: In Active Development - Tsu-Esaki Tunneling Specified
**Priority**: HIGH - Critical for accurate ultrathin device simulation