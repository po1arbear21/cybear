# Schottky Barrier Tunneling Implementation
## Electron and Hole Contributions

### Overview
This document explains the complete implementation of Schottky barrier tunneling for both electrons and holes in the Tsu-Esaki model, as implemented in `src/schottky.f90`.

---

## Physical Model

### Energy Band Diagram
```
Metal                  Semiconductor
------                ----------------------
                      E_c (CB edge) ──────
                            │
E_F = 0 ────────            │ φ_bn (e- barrier)
                            │
                            │ E_g (band gap)
                            │
                      E_v (VB edge) ──────
                            │ φ_bp (h+ barrier)

where: φ_bp = E_g - φ_bn
```

### Key Physics Concepts

1. **Electron Tunneling (Metal → Conduction Band)**
   - Electrons tunnel FROM metal Fermi level TO semiconductor conduction band
   - Barrier height: φ_bn (Schottky barrier for electrons)
   - Energy reference: Metal Fermi level (E_F,metal = 0)
   - Current direction: Negative (electron injection into semiconductor)

2. **Hole Tunneling (Metal → Valence Band)**
   - Holes tunnel FROM semiconductor TO metal
   - Equivalently: Electrons tunnel FROM metal TO valence band states
   - Barrier height: φ_bp = E_g - φ_bn (valence band barrier)
   - Energy reference: Metal Fermi level (E_F,metal = 0)
   - Current direction: Positive (hole extraction from semiconductor)

---

## Mathematical Formulation

### General Tsu-Esaki Current
```
J = ±(m*/(2π²)) ∫₀^φ_barrier T_wkb(E) N(E) dE
```

### WKB Transmission Probability (Triangular Barrier)
```fortran
T_wkb(E) = exp[-(4/3) × √(2m*) × (φ - E)^(3/2) / |F|]
```
- Same for both carriers (carrier-agnostic)
- `φ` is the appropriate barrier (φ_bn for electrons, φ_bp for holes)
- `F` is the electric field magnitude
- `m*` is the tunneling effective mass

### Occupancy Difference N(E)

Energy E is measured from **metal Fermi level** (E_F,metal = 0).

#### Electrons (Conduction Band):
```
N_e(E) = f_metal(E) - f_semi,CB(E)
       = log(1 + exp(-E)) - log(1 + exp(-E - φ_k))
```
where:
- `φ_k` is the semiconductor electrostatic potential at contact (from Poisson)
- Conduction band edge: E_c = φ_bn + φ_k (in biased condition)

#### Holes (Valence Band):
```
N_h(E) = f_metal(E) - f_semi,VB(E)
       = log(1 + exp(-E)) - log(1 + exp(-E + E_g - φ_k))
```
where:
- Valence band edge: E_v = -E_g + φ_k (in biased condition)
- **Note**: `+E_g` term accounts for VB being E_g below metal Fermi level

### Integration Bounds

- **Electrons**: E ∈ [0, φ_bn] (under CB barrier)
- **Holes**: E ∈ [0, φ_bp] where φ_bp = E_g - φ_bn (under VB barrier)

### Prefactor Signs

Critical for correct current direction:

- **Electrons**:
  ```fortran
  prefactor = -m*/(2π²)  ! Negative = injection
  ```

- **Holes**:
  ```fortran
  prefactor = +m*/(2π²)  ! Positive = extraction
  ```

---

## Implementation Details

### Modified Functions

#### 1. `tsu_esaki_integrand(E, p, f, dfdE, dfdp)`
**Location**: `src/schottky.f90:413-523`

**Changes**:
- Added parameter `p(5) = ci` (carrier index)
- Added parameter `p(6) = E_g` (band gap)
- Conditional occupancy calculation:
  ```fortran
  if (ci == CR_ELEC) then
    E_arg_m = E + phi_k            ! CB: E_c = phi_k
  else  ! CR_HOLE
    E_arg_m = E - E_g + phi_k      ! VB: E_v = -E_g + phi_k
  end if
  ```

#### 2. `tsu_esaki_current(phi_b, efield, m_tn, phi_k, ci, E_g, J_tn, dJ_tn_dF)`
**Location**: `src/schottky.f90:525-621`

**Changes**:
- Added input parameter `ci` (carrier index)
- Added input parameter `E_g` (band gap)
- Effective barrier calculation:
  ```fortran
  if (ci == CR_ELEC) then
    phi_eff = phi_b              ! Electrons: φ_bn
  else
    phi_eff = E_g - phi_b        ! Holes: φ_bp = E_g - φ_bn
  end if
  ```
- Carrier-dependent prefactor:
  ```fortran
  if (ci == CR_ELEC) then
    prefactor = -m_tn / (2.0 * PI**2)  ! Injection
  else
    prefactor = +m_tn / (2.0 * PI**2)  ! Extraction
  end if
  ```

#### 3. `calc_schottky_injection_eval(this)`
**Location**: `src/schottky.f90:806-980`

**Changes**:
- Updated call to `tsu_esaki_current`:
  ```fortran
  call tsu_esaki_current(phi_bn_eff, E_field, m_tunnel, &
                         phi_k, this%ci, this%par%smc%band_gap, &
                         J_tn, dJ_tn_dF)
  ```
- Enhanced debug output showing carrier type and barrier heights

---

## Physics Validation

### Expected Behavior

1. **Electron Tunneling (n-type Schottky)**:
   - At forward bias: φ_k < 0 → barrier lowering → increased tunneling
   - Negative current (electrons injected)
   - J_tn magnitude increases with field

2. **Hole Tunneling (p-type Schottky)**:
   - At forward bias: φ_k > 0 → barrier lowering → increased tunneling
   - Positive current (holes extracted)
   - J_tn magnitude increases with field

3. **Barrier Height Relationship**:
   - For symmetric contact: φ_bn + φ_bp ≈ E_g
   - Smaller barrier → more tunneling
   - Typical: φ_bn = 0.3-1.0 eV, E_g = 1.1 eV (Si)

### Critical Sign Checks

| Carrier | Barrier | φ_k (fwd bias) | Current Sign | Prefactor |
|---------|---------|----------------|--------------|-----------|
| Electron| φ_bn    | Negative       | Negative     | -m*/(2π²) |
| Hole    | φ_bp    | Positive       | Positive     | +m*/(2π²) |

---

## Usage Example

### Contact Configuration (INI file)
```ini
[contact.source]
type = schottky
phi_b = 0.35           # eV (electron barrier)
tunneling = true
m_tunnel = 0.26        # Tunneling mass (relative to m0)
```

For this configuration:
- **Electron barrier**: φ_bn = 0.35 eV
- **Hole barrier**: φ_bp = E_g - 0.35 = 0.76 eV (for Si, E_g=1.11 eV)
- **Expected**: More hole tunneling than electron (lower barrier)

### Debug Output
```
Contact 1 (elec) vertex 1:
  J_TE  = -2.345E+02 A/cm^2
  J_tn  = -1.234E+01 A/cm^2
  J_tot = -2.469E+02 A/cm^2
  J_tn/J_TE ratio = 5.26E-02
  E_field = 1.234E+06 V/cm
  phi_k   = -0.123 V
  phi_bn_eff = 0.320 eV

Contact 1 (hole) vertex 1:
  J_TE  = 3.456E+01 A/cm^2
  J_tn  = 4.567E+01 A/cm^2
  J_tot = 8.023E+01 A/cm^2
  J_tn/J_TE ratio = 1.32E+00
  E_field = 1.234E+06 V/cm
  phi_k   = -0.123 V
  phi_bn_eff = 0.320 eV
  phi_bp_eff = 0.790 eV (VB barrier)
```

---

## Testing Recommendations

1. **Unit Tests**:
   - [ ] Verify opposite prefactor signs for electrons/holes
   - [ ] Check barrier height: φ_bp = E_g - φ_bn
   - [ ] Validate occupancy difference at various biases

2. **Physics Tests**:
   - [ ] I-V curves for n-type and p-type Schottky diodes
   - [ ] Tunneling contribution vs field strength
   - [ ] Temperature dependence (kT normalization)

3. **Benchmarks**:
   - [ ] Compare with ATLAS Tsu-Esaki model
   - [ ] Validate against experimental Schottky diode data
   - [ ] Check current continuity in device simulation

---

## References

1. ATLAS Device Simulation Framework, Silvaco
2. Tsu, R. and Esaki, L., "Tunneling in a finite superlattice", Appl. Phys. Lett. 22, 562 (1973)
3. Sze, S.M., "Physics of Semiconductor Devices", 3rd Ed., Wiley (2006)

---

## Notes

- All calculations in **normalized units** (kT/q normalization)
- Electric field includes image-force barrier lowering (IFBL) if enabled
- Tunneling mass can differ from transport effective mass
- Integration uses adaptive quadrature (rtol=1e-6, max_levels=12)

---

**Implementation Date**: 2025-11-19
**Author**: Claude Code (Anthropic)
**Validated**: Build successful, awaiting physics validation
