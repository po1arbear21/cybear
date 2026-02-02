# EBIC Charged Surface Implementation Plan

## Current Status

### What's Implemented (Neutral Surface Model)

The current EBIC simulation includes **neutral surface recombination** matching Haney et al. 2016 (Sec. III):

```
R_surf = S × (n·p - nᵢ²) / (n + p + 2·nᵢ)
```

**Files:**
- `src/srh_recombination.f90:262-449` - Surface SRH implementation
- `src/continuity.f90:211-218` - Integration into continuity equation
- `devices/stem_ebic/stem_ebic_lateral.ini:62-63` - Configuration (S = 10⁶ cm/s)

**Physics captured:**
- Bulk SRH recombination (τ_n, τ_p)
- Surface SRH at x=0 and x=xmax (excluding contacts)
- Collection efficiency η < 100% emerges naturally from DD solution

### What's Missing (Charged Surface Model)

Haney Sec. IV describes charged surfaces with Fermi-level pinning:

```
ρ_surf = q × N_surf/2 × (1 - 2·f_surf)

f_surf = (n_s + p_surf) / (n_s + p_s + n_surf + p_surf)
```

**Not implemented:**
- Surface charge density in Poisson equation
- Fermi-level pinning at surface defect level
- Surface depletion fields (can extend ~150nm for 10¹⁵ cm⁻³ doping)

---

## Why Charged Surfaces Matter for This Setup

| Region | Doping | Surface Depletion Width | Lamella (300nm) |
|--------|--------|-------------------------|-----------------|
| p-side | 10¹⁷ cm⁻³ | ~15 nm | Negligible |
| n-side | 10¹⁵ cm⁻³ | ~150 nm | **50% of thickness!** |

The lightly-doped n-region could be entirely surface-depleted, significantly affecting EBIC collection.

---

## Implementation Plan

### Phase 1: Validation (No Code Changes)

**Goal:** Quantify current model accuracy before adding complexity.

- [ ] Run beam sweep with current neutral surface model
- [ ] Compare EBIC lineshape to Haney Fig. 3 (neutral surface)
- [ ] Check if max efficiency matches analytical formula:
  ```
  η = 1 - (S/2)/(μE + S/2) × (D/z_B)/(μE + D/z_B)
  ```
- [ ] Vary S_surf (10⁴, 10⁵, 10⁶ cm/s) and compare trends

### Phase 2: Charged Surface Implementation

**Goal:** Add Fermi-level pinning and surface charge to Poisson equation.

#### 2.1 New Parameters (device_params.f90)

```fortran
! In semiconductor type:
logical :: charged_surf = .false.
real    :: E_surf = 0.0        ! Surface defect level relative to midgap [eV]
real    :: N_surf = 0.0        ! Surface state density [cm⁻²]
```

#### 2.2 Surface Charge Calculation (new file or extend srh_recombination.f90)

```fortran
type :: surface_charge
  real, allocatable :: rho(:)   ! Surface charge at each surface vertex
  real :: E_surf                ! Defect energy level
  real :: N_surf                ! State density
end type

! Occupancy:
f_surf = (n_s + p_surf) / (n_s + p_s + n_surf + p_surf)

! Where:
n_surf = N_c × exp((E_surf - E_c) / kT)
p_surf = N_v × exp((E_v - E_surf) / kT)

! Charge density:
rho_surf = q × N_surf/2 × (1 - 2×f_surf)
```

#### 2.3 Poisson Equation Modification (poisson.f90)

Add surface charge as boundary source term:
```fortran
! At surface vertices (x=0 or x=xmax, not contacts):
F_poisson = div(ε·grad(φ)) + q·(p - n + N_D - N_A) + δ_surf × ρ_surf/ε
```

**Jacobian entries needed:**
- dρ_surf/dn_s (derivative w.r.t. surface electron density)
- dρ_surf/dp_s (derivative w.r.t. surface hole density)

#### 2.4 Configuration (stem_ebic_lateral.ini)

```ini
! Charged surface parameters
charged_surf = true
E_surf       = 0.0   : eV      ! Midgap pinning (relative to intrinsic level)
N_surf       = 1e11  : 1/cm^2  ! Surface state density
```

### Phase 3: Testing & Validation

- [ ] Compare to Haney Fig. 4 (max EBIC vs beam energy)
- [ ] Test three surface types:
  - n-type surface (E_surf = +0.38 eV): inverted EBIC
  - intrinsic surface (E_surf = 0): reduced max efficiency
  - p-type surface (E_surf = -0.38 eV): standard EBIC
- [ ] Verify position x₀ where surface field reverses direction
- [ ] Check convergence (surface charge couples Poisson ↔ continuity strongly)

---

## Open Questions

1. **Numerical stability**: Charged surfaces create strong field gradients at surface. May need:
   - Finer mesh near surface (dx < 1nm)
   - Damping in Newton solver for surface charge updates

2. **Both surfaces**: With 300nm lamella, both x=0 and x=xmax have charged surfaces. Do they have the same E_surf, N_surf?

3. **Contact interference**: Surface charge should be zero at contact boundaries (already handled by `par%ict%get(idx) == 0` check)

4. **2D vs 3D**: Haney's 2D model divides G_tot by diffusion length L_D. Our 2D model uses per-unit-depth. Verify unit consistency.

5. **Ion damage**: FIB creates Ga implantation damage. Should N_surf vary spatially?

---

## References

- Haney et al. 2016: `/home/yu/cybear/Haney2016_Depletion_Region_Surface_Effects_In_Electron_Beam_Induced_Current_Measurements.pdf`
- Current implementation: `src/srh_recombination.f90`, `src/continuity.f90`
- Device config: `devices/stem_ebic/stem_ebic_lateral.ini`

---

## Files to Modify/Create

| File | Action | Purpose |
|------|--------|---------|
| `src/semiconductor.f90` | Modify | Add charged_surf, E_surf, N_surf parameters |
| `src/device_params.f90` | Modify | Read new parameters from INI |
| `src/surface_charge.f90` | Create | Surface charge calculation + Jacobians |
| `src/poisson.f90` | Modify | Add surface charge source term |
| `src/device.f90` | Modify | Initialize and connect surface charge |
| `devices/stem_ebic/stem_ebic_lateral.ini` | Modify | Add charged surface config |
