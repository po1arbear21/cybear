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

---

## Next Steps

1. **Immediate** (This Week):
   - [ ] Add E₀₀ calculation function
   - [ ] Implement basic TFE model
   - [ ] Add tunneling parameters to device_params
   - [ ] Create tunneling unit tests

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

**Last Updated**: 2025-10-06
**Status**: In Active Development
**Priority**: HIGH - Critical for accurate ultrathin device simulation