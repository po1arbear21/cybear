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

#### 1. **WKB Transmission Probability**

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

#### 2. **Field Emission (FE)**

**Pure field emission current (T→0)**:
```
J_FE = (A*q³E²/8πhφ_b) × exp(-4√(2m*)φ_b^(3/2)/(3qℏE))
```

#### 3. **Direct Tunneling (DT)**

**For ultrathin barriers (< 3 nm)**:
```
J_DT = (q³m₀V²/8πℏ²m*d³) × exp(-4√(2m*q)φ_b × d/(3ℏ))
```

#### 4. **Tsu-Esaki Unified Tunneling Model**

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

**IMPORTANT - Normalization Note**:
The simulator's normalization module uses `2*PI/h` which equals `1/ℏ`, so the normalization implicitly uses ℏ (h-bar), not h directly. Therefore, the coefficient remains **(4/3)** in the WKB formula, not (8/3).

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

### Tsu-Esaki Direct Implementation

**Note**: We are skipping the Padovani-Stratton E00-based TFE model and going directly to the Tsu-Esaki unified approach. This provides:
- Single framework for all tunneling regimes (TE/TFE/FE)
- No artificial transitions between regimes
- Physically rigorous approach using WKB transmission

The Tsu-Esaki model naturally captures all transport through energy integration:
```

### WKB Implementation for Tsu-Esaki

```fortran
  function calc_wkb_transmission(E, phi_b, E_field, m_eff) result(T_wkb)
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

6. **Step 6: Tunneling Infrastructure Implementation** (2025-11-03)
   - ✅ Added `tunneling` and `m_tunnel` fields to `contact` type
   - ✅ Extended `region_contact` type with tunneling parameters
   - ✅ Implemented parameter parsing in `region.f90` with debug output
   - ✅ Added parameter transfer in `device_params.f90`
   - ✅ Verified compilation and parsing functionality
   - Decision: Skip E00/TFE approximate method, go directly to Tsu-Esaki
   - Rationale: Tsu-Esaki naturally includes TE→TFE→FE transitions via WKB

---

## Implementation Progress

### ✅ Completed (2025-11-03)
1. **Infrastructure for tunneling parameters**:
   - Added `tunneling` (boolean) and `m_tunnel` (real) fields to `contact` type
   - Extended `region_contact` type with same fields
   - Implemented parsing from .ini files with debug output
   - Added parameter transfer from regions to device contacts
   - **Tested**: Compiles successfully, ready for physics implementation

### 🚧 Current Implementation Plan

#### **Direct Tsu-Esaki Approach** (Revised - skipping E00/TFE)
**Rationale**: The Tsu-Esaki model with WKB transmission naturally captures all transport regimes:
- Low field/high T → T(E) ≈ 0 except near barrier top → Pure TE
- Medium field/T → Partial tunneling → TFE behavior
- High field/low T → T(E) ≈ 1 for all energies → Pure FE

No need for artificial regime switching based on E00!

## Next Steps

1. **Immediate** (Today):
   - [x] ~~Add m_tunnel and tunneling flag to contact type~~ ✅ DONE
   - [x] ~~Parse parameters from .ini files~~ ✅ DONE
   - [ ] Add Tsu-Esaki functions directly to `schottky.f90`:
     - `calc_wkb_transmission()` for triangular barrier
     - `tsu_esaki_integrand()` for quad interface
     - Use `quad` from quad_m for integration
   - [ ] Integrate into `schottky.f90`:
     - Modify `schottky_injection_mb` to check `tunneling` flag
     - If true: calculate injection from Tsu-Esaki current
     - If false: use existing TE+IFBL
   - [ ] Test with `schottky_diode_tunnel.ini`

### 🔍 Implementation Foundation - Codebase Integration Details

#### 1. **Mathematical Constants Location**
- **Module**: `lib/fortran-basic/src/util/math.f90` (line 29)
- **Import statement**: `use math_m, only: PI`
- **Definition**: `real, parameter :: PI = 4 * atan(1.0)`
- Also available: `PI_16` for quadruple precision if needed

#### 2. **Quadrature Integration Function (quad_m)**
- **Module**: `lib/fortran-basic/src/util/quad.f90`
- **Import statement**: `use quad_m, only: quad`
- **Function signature**:
  ```fortran
  subroutine quad(func, a, b, p, I, dIda, dIdb, dIdp, rtol, err, min_levels, max_levels, ncalls)
    procedure(integrand) :: func        ! Integrand function
    real, intent(in)  :: a, b           ! Integration bounds [a,b]
    real, intent(in)  :: p(:)           ! Parameters passed to integrand
    real, intent(out) :: I              ! Integral result
    real, intent(out) :: dIda, dIdb    ! Derivatives wrt bounds
    real, intent(out) :: dIdp(:)        ! Derivatives wrt parameters
    real, optional, intent(in)  :: rtol ! Relative tolerance (default 1e-13)
    real, optional, intent(out) :: err  ! Error estimate
    integer, optional, intent(in)  :: min_levels  ! Min recursion (default 2)
    integer, optional, intent(in)  :: max_levels  ! Max recursion (default 16)
    integer, optional, intent(out) :: ncalls      ! Function evaluations
  end subroutine
  ```
- **Required integrand interface**:
  ```fortran
  subroutine integrand(x, p, f, dfdx, dfdp)
    real, intent(in)  :: x       ! Integration variable (energy in our case)
    real, intent(in)  :: p(:)    ! Parameters [phi_b, F_field, m_tunnel, V_bias]
    real, intent(out) :: f       ! Function value: T(E) * [f_s(E) - f_m(E)]
    real, intent(out) :: dfdx    ! df/dE (derivative wrt energy)
    real, intent(out) :: dfdp(:) ! df/dp (derivatives wrt parameters)
  end subroutine
  ```
- **Example usage** from `src/current_integral.f90` (line 570):
  ```fortran
  call quad(integrand, eta(1), eta(2), [jj], res, dresdeta(1), dresdeta(2), &
            dresdjj1, max_levels = 8)
  ```

#### 3. **Integration Points in schottky.f90**
- **Line 10-15**: Add imports for `math_m` and `quad_m`
- **Line 150-191**: `schottky_injection_mb()` - Add tunneling check here
- **Line 267**: Insert Tsu-Esaki functions after `schottky_barrier_lowering`
- **Line 451-561**: `calc_schottky_injection_eval()` - E-field already available via `this%efield(normal_dir)%get(idx)`
- **Key available variables**:
  - `par%contacts(ict)%tunneling` - Enable flag (already parsed)
  - `par%contacts(ict)%m_tunnel` - Effective mass ratio (already parsed)
  - `par%contacts(ict)%phi_b` - Barrier height (normalized)
  - `E_field = this%efield(normal_dir)%get(idx)` - Electric field at contact

#### 4. **Minimal Implementation Strategy**
- **NO new files** - Add functions directly to `schottky.f90`
- **Two core functions needed**:
  1. `calc_wkb_transmission(E, phi_b, F, m_tunnel)` - Pure function for WKB
  2. `tsu_esaki_integrand` - Subroutine matching quad interface
  3. `tsu_esaki_current(phi_b, F, V_bias, m_tunnel)` - Wrapper calling quad
- **Use existing infrastructure**:
  - Contact fields `tunneling` and `m_tunnel` already present (Nov 3, 2025)
  - Richardson constant and normalization already handled
  - Sign convention already established (negative for injection)
- **Total code addition**: ~100-150 lines

#### 5. **Available Normalization Functions**
- **Module**: `lib/fortran-basic/src/util/normalization.f90`
- `norm(value, unit_string)` - Convert physical to normalized
- `denorm(value, unit_string)` - Convert normalized to physical
- Energy normalized by `kT/q`, field by `L₀/V_T`

#### 6. **Existing Physical Constants in par object**
- `par%T` - Temperature in Kelvin
- `par%smc%edos(CR_ELEC)` - Effective DOS for electrons (Nc)
- `par%smc%band_gap` - Bandgap (normalized)
- `par%g` - Grid object
- `par%eps()` - Permittivity array

### 📋 Detailed Next Step: Minimal Tsu-Esaki Implementation in schottky.f90

**File**: `src/schottky.f90` (MODIFY EXISTING - No new files!)

**Add to imports section** (around line 10-15):
```fortran
use math_m, only: PI
use quad_m, only: quad
```

**Add functions** (after line 267, following `schottky_barrier_lowering`):
```fortran
! Tsu-Esaki tunneling functions (minimal implementation)

pure function calc_wkb_transmission(E, phi_b, F, m_tunnel) result(T_wkb)
  real, intent(in) :: E, phi_b, F, m_tunnel
  real :: T_wkb
  ! WKB transmission for triangular barrier
  ! T(E) = exp[-4√(2m*)/(3ℏq) * (φb-E)^(3/2) / F]
  if (E >= phi_b) then
    T_wkb = 1.0
  elseif (abs(F) < 1e-10) then
    T_wkb = 0.0
  else
    T_wkb = exp(-(4.0/3.0) * sqrt(2.0*m_tunnel) * (phi_b - E)**1.5 / abs(F))
  endif
end function

! Integrand for quad - matches required interface
subroutine tsu_esaki_integrand(E, p, f, dfdE, dfdp)
  real, intent(in) :: E         ! Energy
  real, intent(in) :: p(:)      ! [phi_b, F, m_tunnel, V_bias]
  real, intent(out) :: f, dfdE
  real, intent(out) :: dfdp(:)
  real :: T_wkb, f_s, f_m

  T_wkb = calc_wkb_transmission(E, p(1), p(2), p(3))
  f_s = 1.0/(1.0 + exp(E))               ! Fermi-Dirac semiconductor
  f_m = 1.0/(1.0 + exp(E + p(4)))        ! Fermi-Dirac metal
  f = T_wkb * (f_s - f_m)

  ! Derivatives (simplified - set to 0 for basic implementation)
  dfdE = 0.0
  dfdp = 0.0
end subroutine
```

**Integration Point**: `src/schottky.f90`

Modify `schottky_injection_mb` (lines 150-191):
- Check `par%contacts(ict)%tunneling`
- If true: Use Tsu-Esaki to calculate injection
- If false: Use existing TE calculation

**Key Implementation Decisions**:
1. **Interface doping**: Use existing semiconductor doping at contact (no new parameter needed yet)
2. **Electric field**: Initially use zero-field approximation, later get from `calc_schottky_injection_eval`
3. **Energy grid**: Start with 100 points, increase if needed for convergence
4. **Temperature**: Always normalized to 1.0 in Cybear's units

2. **Short Term** (Next 2 Weeks):
   - [ ] Add field-dependence to Tsu-Esaki integration
   - [ ] Optimize integration (adaptive grid, caching)
   - [ ] Benchmark against Sentaurus Device
   - [ ] Validate sign convention

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

## Current Status Summary (2025-11-03)

### ✅ What's Working
- **Parameter Infrastructure**: `tunneling` and `m_tunnel` fields are parsed from .ini files
- **Compilation**: All changes compile successfully
- **Debug Output**: Can verify parameters are being read correctly

### 🚧 What's Next
**Create Tsu-Esaki physics module** (`schottky_tunneling_tsu_esaki.f90`):
1. WKB transmission probability function
2. Energy integration with Fermi-Dirac occupancy
3. Hook into existing `schottky.f90` infrastructure

### 🎯 Success Criteria
- When `tunneling = true`: Should see higher current than pure TE
- When `tunneling = false`: Should match existing behavior exactly
- Physical validation: Current enhancement should scale with doping and field

---

**Last Updated**: 2025-11-03
**Status**: In Active Development - Infrastructure Ready, Physics Implementation Next
**Priority**: HIGH - Critical for accurate ultrathin device simulation