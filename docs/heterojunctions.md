# Heterojunctions

This document specifies how cybear models heterojunctions: how the user describes
multi-material devices, what physical convention the simulator applies at material
interfaces, how heterojunctions are represented on the mesh, and what the v1
implementation does and does not support.

**Status (2026-05-13).** Implementation in progress on the `heterojunction-v2`
branch. Phase A (multi-material catalog `device_params%mat(:)`, unified
`[semiconductor]` INI parser, region-to-material resolution, v1 uniformity
enforcement for `dist`/`stab`) is being re-implemented against the post-`dd-arch-sync`
codebase. **Not yet committed:** per-vertex band-edge fields and the
Scharfetter–Gummel `dpot` modification — until Phase B lands, multi-material
INIs parse and resolve correctly but produce identical physics to single-material
runs (no Δχ contribution to the SG flux). The schema described below is the
contract; the parser will enforce it.

## Overview

A heterojunction is a junction between two materials with different band
structures. Cybear represents this by allowing band-structure parameters
(conduction-band edge `Ec`, valence-band edge `Ev`, bandgap `Eg`, CB-edge
offset `dEc`, density-of-states `Nc`/`Nv`) to vary spatially across the device.

Targets enabled by this feature:

- Si/SiGe channels in BiCMOS and HBT structures.
- GaAs/AlGaAs HEMTs.
- Perovskite VFETs with multi-layer stacks (Project A07).
- 2D-material vdW heterostructures (MoS2/WSe2 and related).

Single-material devices (the entire pre-heterojunction test suite) remain
fully supported and produce bit-identical output.

## Band offset convention

Each `[semiconductor]` block declares a CB-edge offset `Ec_offset` (stored on
`type semiconductor` as `dEc`) on a user-chosen reference scale. The simulator
only ever consumes **differences**:

    Delta Ec = dEc(material_1) - dEc(material_2)

The choice of zero is a per-INI user convention. The reference can be any
material's CB edge, vacuum level, or just `0` on whichever material the user
calls reference. Three writers of the same heterojunction stack can pick three
different zeros and obtain bit-identical simulation results.

The valence-band offset is derived: `Delta Ev = Delta Ec - Delta Eg`. There is
no separate user input for VB offset; materials characterized primarily by VB
offset are converted by the user before being written into the INI.

This convention does **not** apply the Anderson electron-affinity rule.
`Ec_offset` is meant to carry experimentally measured or DFT-computed band
offsets directly. Anderson's rule is rarely accurate in practice (it ignores
interface dipoles, strain, charge transfer, and Bardeen pinning); the direct-
offset convention sidesteps it.

## Mesh representation

Band-structure quantities are stored **per vertex** on the simulation mesh:

| Quantity            | Field                               | Type                          |
| ------------------- | ----------------------------------- | ----------------------------- |
| Conduction-band edge `Ec(v)` | `device_params%band_edge_v(1)` | `grid_data_real`, vertex      |
| Valence-band edge `Ev(v)`    | `device_params%band_edge_v(2)` | `grid_data_real`, vertex      |
| Effective DOS `Nc(v)`        | `device_params%edos_v(1)`      | `grid_data_real`, vertex      |
| Effective DOS `Nv(v)`        | `device_params%edos_v(2)`      | `grid_data_real`, vertex      |
| Bandgap `Eg(v)`              | `device_params%band_gap_v`     | `grid_data_real`, vertex      |
| CB-edge offset `dEc(v)`      | `device_params%dEc_v`          | `grid_data_real`, vertex      |

Vertex storage matches the existing consumer pattern: density, quasi-Fermi
level, and contact boundary conditions are all evaluated at vertices, and the
Scharfetter-Gummel current along an edge requires only the *difference* of
adjacent-vertex band-edge values. There is therefore no separate per-edge
band-edge field.

A single vertex on a material interface holds **one** (`Ec`, `Ev`, `Nc`, `Nv`,
`dEc`) tuple, taken from the vertex's owning region per the existing
region-loop convention (see `init_doping` in `device_params.f90` for the
template). The discontinuity is realized across the *adjacent edge*, not at
the interface vertex itself. **Mesh refinement near interfaces is the user's
responsibility.**

## INI schema

Heterojunction devices are described by:

1. One or more `[semiconductor]` blocks (the material catalog).
2. `[transport]` regions that reference a material by name.

This is the same pattern the simulator already uses for `[contact]` blocks
(stored in `device_params%contacts(:)` with `contact_map` for name lookup).
Existing single-material INIs continue to work unchanged — the new schema
only *extends* the legacy single-`[semiconductor]` form.

### Material catalog

Each `[semiconductor]` block declares one material. The optional `name` key
is the handle used by `[transport]` regions; it is required when more than
one block is present:

```ini
[semiconductor]
name      = Si
E_gap     = 1.12        ! eV
Ec_offset = 0.0         ! eV    (this material is the reference)
N_c0      = 2.86e19     ! cm^-3 at 300 K (scaled by T**1.5 internally)
N_v0      = 3.10e19     ! cm^-3 at 300 K
electrons = true    ! read once from mat_default; controls ci0/ci1 globally
holes     = true
dos   = parabolic   ! uniformity-checked across blocks (v1)
dist  = maxwell     ! uniformity-checked (v1 deferral)
stab  = sg          ! uniformity-checked (v1 deferral)
mob_n = 1450        ! cm^2/Vs
mob_p = 450         ! cm^2/Vs
! Caughey-Thomas params (low-field mobility per material; saturation params
! are still global in v1; see "Deferred" below).
mob_min_n = 92
mob_min_p = 54
N_ref_n   = 1.3e17
N_ref_p   = 2.4e17
alpha_n   = 0.91
alpha_p   = 0.81

[semiconductor]
name      = SiGe20
E_gap     = 0.97
Ec_offset = -0.05       ! eV    (50 meV below Si CB; illustrative)
N_c0      = ...
N_v0      = ...
```

### Transport regions

Each `[transport]` region adds a single optional `material` key that
references a `[semiconductor]` block by name:

```ini
[transport]
xmin = 0
xmax = 50
material = Si

[transport]
xmin = 50
xmax = 100
material = SiGe20
```

### Resolution rules

| INI shape | Behavior |
| --------- | -------- |
| 1 `[semiconductor]` block, no `name =` | Anonymous default. `mat(1)%name%s = "default"`. Transport regions resolve to `mat_default = 1`. Every existing single-material INI works unchanged via this path. |
| 1 `[semiconductor]` block, with `name =` | Named, still the lone default. Transport regions may specify `material = <name>` redundantly. |
| N `[semiconductor]` blocks, all with `name =` | Multi-material mode. First block becomes `mat_default`. Transport regions without `material =` resolve to default; with `material =` resolve via `mat_map`. |
| N `[semiconductor]` blocks, ≥1 missing `name =` | **Hard-error:** `multiple [semiconductor] sections require every block to declare name=` |

`Ec_offset` defaults to `0.0` when absent (single-material INIs never evaluate
any ΔEc). For heterojunction simulations, `Ec_offset` should be supplied
explicitly on every block on a consistent user-chosen reference scale.

## Discretization

The Scharfetter-Gummel argument along an edge between vertices `v1` and `v2`
is

    dpot = - ch * (pot(v2) - pot(v1)) + (be(v2) - be(v1))

where `ch` is the per-carrier sign (`-1` for electrons, `+1` for holes), `pot`
is the electrostatic potential, and `be` is the per-vertex band-edge field
(`Ec` for electrons, `Ev` for holes). When `be` is uniform this collapses to
the legacy single-material expression `dpot = -ch * (pot(v2) - pot(v1))`,
guaranteeing bit-identical regression behavior. The reference for the
augmented form is Selberherr §4.3 / Yu & Dutton.

DOS normalization in the SG flux uses per-vertex `Nc(v)`/`Nv(v)` rather than
a scalar:

    n(v) = dens(v) / edos(v)

The `dens/edos` denormalization is preserved in the final flux, which is the
carrier-charge-conserving convention.

The Maxwell-Boltzmann shortcut for intrinsic carrier density becomes
vertex-local:

    ni(v) = sqrt(Nc(v) * Nv(v)) * exp(- 0.5 * Eg(v))

## Assumptions (v1)

1. **Static band offsets, no interface dipole.** `Delta Ec` follows directly
   from each material's user-declared `Ec_offset` value. No spontaneous-
   polarization corrections, no fixed interface charge, no Bardeen-style
   pinning. Users wanting to model these effects must fold them into the
   `Ec_offset` values themselves (e.g., shift one material's offset by the
   measured dipole contribution).
2. **Per-vertex band edges are well-defined at the interface.** Each
   interface vertex inherits one material's parameters; the discontinuity is
   carried by the adjacent edge.
3. **Stabilization method is uniform across the device.** `smc%stab`
   (Scharfetter-Gummel vs. enhanced-diffusion vs. exact integration) is read
   once from `mat_default`. Per-material `stab` values trigger a hard-error
   during the v1 uniformity check.
4. **Distribution function is uniform.** Maxwell-Boltzmann vs. Fermi-Dirac is
   a global choice (read from `mat_default`, uniformity-checked). Mixing
   degenerate and non-degenerate materials in one device is unsupported.
5. **No interface trap states.** Surface and interface traps are not
   modeled. The existing `incomplete_ionization` machinery applies only to
   bulk dopants.
6. **No bias-dependent band-edge shifts.** `Ec`, `Ev`, `Nc`, `Nv` are
   static. Bandgap-narrowing-with-doping (Slotboom) and bandgap-with-T are
   not implemented (also true today in single-material mode).

## Deferred / not supported in v1

Each entry below is a real limitation. They are tracked separately and will
be addressed in follow-up work as use cases demand.

### Per-material stabilization method

`smc%stab` (one of `STAB_SG`, `STAB_ED`, `STAB_EXACT`) is a global
device-wide setting in v1. Cross-material upwind logic in
`current_density.f90` would add branching in the hot path and has no
compelling test case yet.

**Workaround:** pick the strictest-needed method globally (typically
`STAB_ED` or `STAB_EXACT` for any device that mixes degenerate and
non-degenerate carriers).

### Per-material distribution function

`smc%dist` is global. Crossing the degenerate / non-degenerate boundary
mid-device is an open research question, not a v1 deliverable.

**Workaround:** use the more demanding distribution globally (`DIST_FERMI`).

### Per-material velocity-saturation parameters

Caughey-Thomas low-field mobility (`mob`, `mob_min`, `mob_max`, `N_ref`,
`alpha`) is wired per-material via the existing `mob0(:,:,ci)` per-cell
field. **Saturation parameters** (`v_sat`, `beta`) inside
`mobility.f90`'s field-dependent kernel still read from the global `smc`
(populated from `mat_default`).

**Workaround:** equilibrium and low-bias heterojunction simulations are
unaffected. For high-bias work, set `v_sat`/`beta` to the limiting
material's values until the kernel is refactored.

### Per-material incomplete ionization

The Altermatt-Schenk parameters (`ii_E_dop0`, `ii_N_crit`, `ii_g`,
`ii_dop_th`) and the Poole-Frenkel / tunneling extensions
(`ii_pf`, `ii_pf_a`, `ii_tau_tun`, `ii_m_tun`) are read once from
`mat_default`. The single-material assumption is wired through
`ionization.f90` and the ohmic-contact charge-neutrality solve in
`contact.f90`.

**Workaround:** disable `incomp_ion` in heterojunction tests, or use
materials with similar dopant ionization energies (e.g. Si and SiGe) where
the global value is a reasonable approximation for both.

### Interface trap states and fixed interface charge

There is no data structure for interface charge in v1.

**Workaround:** approximate with a thin highly-doped transport region at the
interface.

### Polar interfaces (nitrides, ferroelectric stacks)

Anderson alone is wrong for GaN/AlGaN, ferroelectric/semiconductor stacks,
and polar 2D heterostructures (spontaneous + piezoelectric polarization
charge at the interface). These devices need an additional fixed interface
charge contribution that v1 does not provide.

**Workaround:** none in v1. Out of scope; revisit with a separate plan.

### Quantum confinement at the heterojunction

Cybear's quantum-confinement support (subband structure in nanowires and
2D systems) is orthogonal to this feature. A combined heterojunction +
quantum-confinement workflow is not validated.

## Equilibrium current floor

A heterojunction equilibrium simulation in cybear does **not** reach the
machine-precision zero (~1e-30 A) that single-material runs achieve. The
realistic floor is ~1e-20 to 1e-21 A. The cause is structural to the
discretization, not a bug:

- In single-material bulk at equilibrium, every interior edge has
  `dpot = 0` and `n(1) = n(2)` to bit-precision (Newton converges all
  bulk values to identical floats). The SG kernel returns
  `Bern(0)·n − Bern(0)·n = n − n = 0` exactly, and continuity propagates
  this through all edges.
- At a Si/SiGe interface edge, `dpot ≠ 0` (carries the `Δbe` term) and
  `n(1) ≠ n(2)` (different per-vertex `Nc` normalization). The SG kernel
  evaluates `Bern(-dpot)·n(1) − Bern(dpot)·n(2)`, which is *algebraically*
  zero at flat quasi-Fermi but has a floating-point cancellation residual
  of order `eps · |Bern·n| ≈ 1e-14` in normalized units (1e-20 to 1e-21 A
  after denormalization). Continuity then propagates this residual to all
  edges in the chain, giving the same magnitude at the contacts.

Eliminating this floor requires the discretization used by Sentaurus
Device, Sesame, and similar TCAD tools — **the quasi-Fermi potential is a
primary Newton variable** rather than a derived quantity. With
`iref_e(contact) = V_applied` enforced as a Dirichlet BC at the metal
contacts, every vertex stores its `iref` as a literal float. At V=0
equilibrium, Newton drives `iref(v) = 0` identically across the device,
and the Sesame-form SG flux

    J_n ∝ M(ψ_n(2), ψ_n(1)) · [exp(EFn(2)) − exp(EFn(1))]

vanishes exactly because `exp(0) − exp(0) = 0` in floats. Cybear's current
Newton system uses `(pot, dens)` as primary variables; switching to
`(pot, iref)` is a substantial refactor and is **not** in v1 scope.

For the v1 implementation, the per-vertex band-edge fields and SG
modification are correct (verified by detailed balance derivation); the
1e-21 A floor is the FP precision limit of the dens-based SG kernel at a
non-zero `dpot` interface edge.

## Validation

A device is considered correctly heterojunction-modeled when:

1. **Single-material regression.** `fargo run schottky_test`,
   `fargo run debug`, and `fargo run pn_1e18_1e17` produce bit-identical
   output (1e-12 relative tolerance) with the new code path.
2. **Synthetic uniform Ec.** Setting both regions to the same material
   produces fluxes identical to the single-material baseline.
3. **Si / SiGe equilibrium test (`devices/het/si_sige.ini`).** At V=0:
   - `Ec(x)` shows a step `Delta Ec = Ec_offset(Si) - Ec_offset(SiGe20)` at
     the interface vertex (i.e. the user-declared band-offset values).
   - The quasi-Fermi level is flat across the device (max deviation
     < 1e-6 eV).
   - Integrated electron density matches the analytical N+/N step
     prediction.
4. **Newton convergence.** The new test converges in <= 15 iterations at
   V=0 (parity with single-material at the same bias). If not, the SG
   `dpot` sign is the prime suspect; bisecting by zeroing the
   `(be(v2) - be(v1))` term must recover the legacy behavior in the
   homogeneous limit.

## References

- S. M. Sze, *Physics of Semiconductor Devices*, ch. 1 (heterojunction band
  alignment; Anderson rule and its limitations).
- S. Selberherr, *Analysis and Simulation of Semiconductor Devices*,
  §4.3 (Scharfetter-Gummel with position-dependent band edges).
- Z. Yu, R. W. Dutton, "SEDAN III - A Generalized Electronic Material
  Device Analysis Program" (DOS-discontinuity convention at
  heterojunction edges).
