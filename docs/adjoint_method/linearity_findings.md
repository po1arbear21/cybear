# EBIC Linearity Test — Findings

Results from the forward-only linearity suite that gates the EBIC adjoint-method
study. Test spec: `linearity_test.md` (same directory). Driver:
`src/linearity_test.f90`. Device + run INIs: `devices/linearity_test/`.
Launched as a single fargo target: `fargo run linearity_test`.

## Purpose

The adjoint method for EBIC collection maps is exact only for the linearized
drift-diffusion system around the dark operating point u₀. Before implementing
adjoint infrastructure, we needed to verify empirically whether
F: G → I_EBIC is linear in I_beam at the experimental operating point I_exp.
If it isn't, one precomputed collection map cannot describe the response and
the adjoint approach fails for that geometry.

## What was done

**Test 1 (homogeneity scaling)** was implemented. At a fixed beam position r_A,
I_beam is swept across five decades (10⁻² to 10² × I_exp) and the terminal
currents recorded per decade. Linearity is confirmed if I_EBIC(k) / I_beam(k)
is constant across decades.

Tests 2 (superposition) and 3 (small-signal ratios) from the spec are not yet
implemented.

The driver runs each (device, run) pair as a separate subprocess via the
`--suite` mode so that Fortran module-singleton state doesn't leak between
cases, and aborts on first failure.

## Test configurations

All cases short-circuit bias (V=0), nonzero I_exp = 0.0825 nA/μm (matches
target STEM beam current), same I_beam decade sweep.

| case             | geometry        | contacts                | SRH   | surf recom | beam profile | r_A              |
|------------------|-----------------|-------------------------|-------|------------|--------------|------------------|
| lateral          | 2D cross-section, 300 nm × 5 μm  | ohmic, y-ends  | 2e-11 s | off | line  | 2000 nm (0.5 μm n-side of junction at 1500) |
| lateral_surf     | same            | ohmic, y-ends  | 2e-11 s | **on**  | line  | 2000 nm |
| lateral_schottky | same but 400 nm thick | Schottky (0.65 eV), y-ends | 5.9e-6 s | off | line | 2000 nm |
| flat_ohmic       | 2D top-view, 14 × 8.5 μm, no slits | ohmic corners (P top-right, N bottom-left) | 2e-11 s | off | **point** | 7000 nm (at junction y=7) |
| flat_schottky    | same top-view, no slits | Schottky (0.42 eV) corners | 2e-11 s | off | point | 7000 nm |
| slit             | 2D top-view with FIB voids, same 14 × 8.5 μm  | Schottky (0.42 eV), full-edge spans | 2e-11 s | off | point | 7000 nm |

The flat_ohmic, flat_schottky, and slit variants were designed as an
isolation matrix to separate (geometry × contact type × SRH) effects.

## Results — isolation matrix

Four variables were explored on the top-view (point-beam) geometries:
slit voids present/absent, contact type (ohmic/Schottky), SRH on/off.

| slit voids | contacts  | SRH | Test 1 result                    |
|------------|-----------|-----|----------------------------------|
| no         | ohmic     | on  | **linear** (flat_ohmic)          |
| no         | Schottky  | on  | **nonlinear** (flat_schottky+SRH) |
| no         | Schottky  | off | **linear** (flat_schottky−SRH)   |
| yes        | Schottky  | on  | **nonlinear** (slit)              |
| yes        | Schottky  | off | **linear** (slit−SRH)             |

The cross-section (line-beam) lateral variants were all linear regardless of
contact type or SRH state.

## Conclusion

**Linearity fails when and only when SRH recombination AND Schottky contacts
are both present on the point-beam top-view geometry.** Removing either
ingredient restores linearity. The slit voids are not required — `flat_schottky`
(no voids) has the same nonlinearity as `slit`. The slit geometry is incidental
to the mechanism.

### Why lateral (line beam) stays linear regardless

Line profile spreads the beam generation over the full x-column at y=r_A
(~60 grid vertices in the lateral geometry). Local injection density Δn at
each vertex is ~60× lower than if the same G_tot were deposited at one point.
Estimated Δn/n₀ ≈ 10⁻³ at I_exp for the lateral cases — deep in the
low-injection regime where SRH is linear in Δn and Schottky BCs see a
negligible density perturbation. Neither ingredient matters here.

### Why point beam + SRH + Schottky is nonlinear

Point profile concentrates all G_tot into a single grid vertex. With the
refined mesh at the junction (~100 nm × 2 nm at y=r_A), local injection
density reaches Δn/n₀ ≈ 0.2 already at I_exp, and ≈ 2 at I_exp × 10. This
crosses the SRH transition region (Δn ~ n₀), where the effective lifetime
shifts from τ_minority (low injection) toward τ_n + τ_p (high injection).
At the same time, the Schottky Robin boundary condition has an exponential
dependence on local ψ, which couples amplified carrier-density variations
back into the potential. Together these produce a bias-point-dependent
collection efficiency that is not captured by a single linear operator.

Removing SRH eliminates the density-dependent sink in the continuity
equation; removing Schottky eliminates the exponential density-to-current
coupling at the contacts. Either removal decouples the mechanism.

Notably, neither is nonlinear *by itself* at the operating point. The
nonlinearity requires them together. This pairing is specific to devices
where the photogenerated minority carriers must both (i) survive SRH on
their way to the contact, and (ii) leave the device via a Schottky-BC
contact — exactly the topology of the FIB slit lamella in the SISPAD
draft.

## Caveats on interpretation

- Test 1 alone does not quantify Δu/u₀ at I_exp; it only detects whether
  the response *ratio* varies across injection levels. Test 3 (small-signal
  ratios, not yet implemented) would give a direct yes/no on whether the
  linearization is accurate at I_exp specifically.
- The point-beam model is conservative. Cybear's built-in Gaussian profile
  (σ₀ = 0.1 nm at entrance, growing as σ(x) = √(σ₀² + (αx)²) with α = 10 mrad)
  is closer to a real STEM beam. At 500 nm lamella depth σ ≈ 5 nm; that still
  concentrates most generation in one grid cell of a 100 nm mesh, but distributes
  enough into neighbors to substantially reduce peak Δn/n₀. The nonlinearity
  magnitude we observed is therefore an upper bound.
- The point-beam run on the slit uses `lamella_t = 500 nm` but the beam's y=7
  line crosses two vacuum slit voids; cybear's generation model uses a fixed
  Si stopping power over the full lamella_t regardless, overestimating G_tot
  by ~8%. This is a constant multiplier and does not affect linearity, only
  the absolute collection-efficiency numbers.
- The dark-current offset is not subtracted from I_EBIC before the ratio is
  computed. On devices with any geometric or contact asymmetry at V=0 (e.g.,
  the slit with two dissimilar voids), low-decade ratios are contaminated by
  the dark offset rather than being a clean high-injection diagnostic.
- PASS/FAIL reporting in the driver currently only checks subprocess exit
  status (Newton convergence), not whether the ratio is actually flat. The
  linearity call has to be made by eyeballing the printed table.

## Implications for the adjoint EBIC method

The Donolato adjoint formulation is exact only when
I_EBIC = c⊤ A⁻¹ s holds for the linearized system at u₀. Given:

- Lateral cross-section devices: linear at I_exp → adjoint valid.
- Top-view point-beam devices with Schottky + SRH: **nonlinear at and above
  I_exp** → standard Donolato adjoint is not valid at the simulated
  operating point for the slit-W/Pt production device.
- Top-view point-beam devices with ohmic *or* without SRH: linear → adjoint
  valid for those variants, but they aren't the physical production device.

Three things need resolution before concluding whether the adjoint is
workable for the paper's production device:

1. **Beam profile realism**: re-run slit and flat_schottky with
   `beam_dist = "gaussian"` to see whether a more physical STEM-like
   injection density ever stays inside the linear regime at I_exp.
2. **Direct small-signal measurement**: implement Test 3, compute
   max(Δn/n₀) and max(Δψ/V_T) at I_exp for each device, and compare against
   the spec's < 0.1 criterion. That is the honest yes/no on linearization
   validity at the operating point, independent of Test 1's edge artifacts.
3. **Injection-level window**: even if the adjoint fails at the simulated
   I_exp, it may be valid at ~0.1 × I_exp. If the experiment can be
   performed at lower beam current without SNR loss, that's a path forward.

If all three checks fail, the adjoint approach should not be used for the
slit-W/Pt paper and the study should revert to one forward solve per beam
position. The findings here would themselves be publishable as a clean
isolation of why the SRH-Schottky interaction breaks linearization in this
device class.

## Status

- Test 1 implemented, run on 6 cases, isolation matrix completed.
- Test 2 and Test 3: stubs only.
- Adjoint infrastructure: not started (blocked on linearity verdict
  at I_exp for the production device — needs Gaussian profile and/or Test 3
  before proceeding).
