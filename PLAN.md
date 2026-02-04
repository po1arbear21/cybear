# EBIC Charged Surface Implementation Plan

## Current Status: Phase 2 Complete

### What's Implemented

#### 1. Beam Generation (`src/beam_generation.f90`)

Three beam profile options:
- **LINE**: Delta function at grid node (default, use Y_BEAM for grid-independence)
- **POINT**: 2D localized at (beam_x, beam_y)
- **GAUSSIAN**: Depth-dependent cone with σ(x) = √(σ₀² + (α·x)²)

Hardcoded STEM parameters (200 keV in Si):
- Stopping power: 5.21 MeV/cm
- Energy per e-h pair: 3.6 eV
- Semi-convergence angle: 10 mrad
- Probe size: σ₀ = 0.1 nm

#### 2. Surface SRH Recombination (`src/srh_recombination.f90:262-449`)

```
R_surf = S × (n·p - nᵢ²) / (n + p + 2·nᵢ)
```

- Applied at x-boundaries (x=0, x=Lx), excluding contacts
- Area weighting via y-direction adjoint lengths
- Jacobians: dR/dn, dR/dp

#### 3. Charged Surface / Fermi-Level Pinning (`src/surface_charge.f90`)

Following Haney et al. 2016 Sec. IV:
```
ρ_surf = q × N_surf/2 × (1 - 2·f_surf)
f_surf = (n_s + p_trap) / (n_s + p_s + n_trap + p_trap)

where:
  n_trap = N_c × exp((E_surf - E_c) / kT)
  p_trap = N_v × exp((E_v - E_surf) / kT)
```

- **2D only** (enforced with error check)
- Surface vertices: x=0 and x=Lx, excluding contacts
- Coupled to Poisson via `jaco_scharge` in `poisson.f90`
- Jacobians: dρ_surf/dn, dρ_surf/dp

#### 4. Ramo-Shockley Current Collection (`src/ramo_shockley.f90`)

- Fundamental solutions computed once per device (Laplace with unit potential)
- Terminal current: I = ∫ J · ∇φ_fundamental dV
- Capacitance matrix for contact-contact coupling

#### 5. Equation Integration (`src/device.f90`)

```
sys_nlpe (Non-linear Poisson):
├─ calc_scharge  → ρ_surf(n, p)     [if charged_surf]
├─ poisson       → φ(ρ_vol, ρ_surf)
└─ calc_rho      → ρ_vol(n, p)

sys_dd (Drift-Diffusion per carrier):
├─ continuity    → ∂n/∂t + ∇·J - G + R = 0
├─ calc_cdens    → J(φ, n)
├─ calc_bgen     → G_beam           [if has_beam_gen]
├─ calc_srh      → R_bulk           [if srh]
└─ calc_surf_srh → R_surf           [if surf_recom]
```

---

## Known Issue: Surface SRH vs Charged Surface Inconsistency

**Problem**: Surface SRH and charged surface use different trap level assumptions.

| Model | Denominator | Assumption |
|-------|-------------|------------|
| Surface SRH | `n + p + 2·nᵢ` | Midgap trap (E_t = E_i) |
| Charged Surface | `n + p + n_trap + p_trap` | General trap (E_t = E_surf) |

**Impact**: When `E_surf ≠ 0`, the two models are inconsistent:
- Charged surface correctly uses E_surf-dependent n_trap, p_trap
- Surface SRH ignores E_surf and assumes midgap (n₁ = p₁ = nᵢ)

**Current workaround**: Use `E_surf = 0.0 eV` (midgap), which makes both models consistent.

**Proper fix**: Update `calc_surface_srh` to use:
```fortran
n1 = ni * exp((E_surf - E_i) / kT)  ! or equivalently Nc * exp((E_surf - Ec) / kT)
p1 = ni * exp((E_i - E_surf) / kT)  ! or equivalently Nv * exp((Ev - E_surf) / kT)
denom = n + p + n1 + p1
```

---

## Current Limitations

### Physical Limitations

1. **2D only for charged surface** - No 3D surface element handling
2. **Constant stopping power** - No depth-dependent dE/dz (ESTAR data not used)
3. **No multiple scattering** - Probe profile unrealistic for thick samples
4. **Static E_surf** - Not coupled to band bending at surface
5. **No recombination through surface states** - Only occupancy modeled, not tunneling/thermionic

### Numerical Limitations

6. **Grid-dependent beam** - Must specify Y_BEAM for grid-independence
7. **Contact surface states ignored** - Assumed ohmic at contact boundaries

---

## Device Configurations

### `devices/stem_ebic/stem_ebic_lateral.ini` (Primary test case)

```
Geometry:
  - Lamella: 300 nm (x) × 5000 nm (y)
  - Junction: y = 1500 nm
  - p-region: N_a = 1e18 cm⁻³ (y < 1500 nm)
  - n-region: N_d = 1e15 cm⁻³ (y > 1500 nm)

Surface Physics:
  - surf_recom = true, S_surf = 1e6 cm/s
  - charged_surf = true, E_surf = 0.0 eV, N_surf = 1e11 cm⁻²

Beam:
  - I_beam = 0.0825 nA/μm
  - beam_dist = "line"
  - Y_BEAM = 0, 1000, 2500, 5000 nm (mandatory grid nodes)
```

### Surface Depletion Impact

| Region | Doping | Surface Depletion Width | Lamella (300nm) |
|--------|--------|-------------------------|-----------------|
| p-side | 10¹⁸ cm⁻³ | ~1.5 nm | Negligible |
| n-side | 10¹⁵ cm⁻³ | ~150 nm | **50% of thickness!** |

The lightly-doped n-region may be significantly surface-depleted.

---

## Next Steps

### Phase 3: Validation

- [ ] Run beam sweep with current implementation
- [ ] Compare EBIC lineshape to Haney et al. 2016 Fig. 3-4
- [ ] Test surface types: E_surf = -0.38, 0.0, +0.38 eV
- [ ] Verify convergence behavior near surfaces

### Phase 4: Fix Surface SRH Consistency

- [ ] Add E_surf-dependent n1, p1 to `calc_surface_srh`
- [ ] Share n_trap, p_trap computation with `calc_surface_charge`
- [ ] Or: read E_surf in surface SRH init and compute locally

### Phase 5: Enhancements (Future)

- [ ] Depth-dependent stopping power (ESTAR tables)
- [ ] 3D charged surface support
- [ ] Collection efficiency tracking (carrier fate analysis)
- [ ] Transient solver for time-resolved EBIC

---

## File Summary

| File | Status | Purpose |
|------|--------|---------|
| `src/beam_generation.f90` | ✓ Complete | Beam generation (line/point/Gaussian) |
| `src/srh_recombination.f90` | ⚠ Needs fix | Bulk + surface SRH (midgap assumption) |
| `src/surface_charge.f90` | ✓ Complete | Charged surface (Fermi-level pinning) |
| `src/poisson.f90` | ✓ Complete | Poisson with surface charge coupling |
| `src/continuity.f90` | ✓ Complete | Carrier continuity with surface SRH |
| `src/ramo_shockley.f90` | ✓ Complete | Terminal current calculation |
| `src/device.f90` | ✓ Complete | Equation assembly and solver setup |
| `src/device_params.f90` | ✓ Complete | Parameter loading (charged_surf, E_surf, N_surf) |

---

## References

- Haney et al. 2016: "Depletion Region Surface Effects in Electron Beam Induced Current Measurements"
- Device config: `devices/stem_ebic/stem_ebic_lateral.ini`
- CLAUDE.md: Project overview and physics documentation
