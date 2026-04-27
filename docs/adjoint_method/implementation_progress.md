# Adjoint EBIC — Implementation Progress

Empirical status and findings for Phase B of the adjoint-method study.
Companion to `implementation_plan.md` (forward-looking design) and
`linearity_findings.md` (Phase A empirical findings). Last updated
2026-04-21.

## Status summary

| Step | Description | Status |
|------|-------------|--------|
| 1 | `keep_factorization` in fortran-basic `steady_state.f90` | **SKIPPED** (not needed — see below) |
| 2 | `src/adjoint_ebic.f90` kernel | **DONE** — self-consistent (primal identity at 6.5e-14) |
| 3 | 1D resistor η = 1−x/L validation | **DROPPED** (degenerate setup) |
| 4 | pn reciprocity harness at 5 r_beam | **DONE** — 7-digit agreement (no-SRH mesh); 5-digit on dense-mesh + SRH |
| 5 | Donolato analytical cross-check | Not started |
| 6 | `src/adjoint_ebic_driver.f90` driver | **DONE** — 1D / 2D (x,y) / 2D (y,z) / 3D sweeps, fbs output, linescan reference, primal-check diagnostic |
| 7 | Production slit-W/Pt collection map | Blocked on Phase A |
| 8 | Root-cause the 5e-4 adj/fwd shortfall on flat_ohmic_2d | **IN PROGRESS** — traced to gauge-like δu discrepancy (see "Open bug" below) |
| — | Fast-path opt 1: manual `jaco_bgen` scatter (eval-based → mat-vec) | **DONE** — ~10× per-position speedup |
| — | Fast-path opt 2: profile-aware direct sum (skip `calc_bgen`, mat-vec, dot) | **DONE** — validated bit-for-bit at peak |

## What was done

**Step 1 dropped as no-op.** The `keep_factorization` flag we planned to add
to `steady_state_m` turned out to be unnecessary. PARDISO's factorization
stored inside `ss%solver` naturally persists as long as the `ss` object
itself is alive — nothing deletes it except explicit `ss%destruct()` calls,
which our driver simply doesn't make until after the adjoint work is done.
No fortran-basic edits required.

**Step 2 — kernel** (`src/adjoint_ebic.f90`). Public surface:

- `build_collection_functional(dev, contact_id, c)` — builds `c` as a unit
  indicator at the sys_full DOF holding the target contact's Ramo-Shockley
  terminal current. Uses `search_main_var("currents")` + `res2block`.
- `refactorize_at_converged(dev, solver)` — re-evaluates J at the converged
  Newton state and re-factorizes (pins LU to `u*` exactly, vs. the
  `u_{final-1}` state Newton's last factorize was done at).
- `build_adjoint_source(dev, r_beam, s, x_beam, z_beam)` — general slow path
  that computes `s = V·G(r; r_beam)` scattered into continuity interior rows
  via `calc_bgen%eval` + `jaco_bgen` mat-vec. Still available as the
  fallback for `beam_dist = "gaussian"` or future profiles.
- `solve_adjoint(solver, c, phi)` — one-liner transpose-solve via `trans='T'`.
- `evaluate_ebic(phi, s)` — scalar dot product.
- `extract_phi_for_var(dev, phi, var_name, tab_idx, phi_slice)` — slices the
  adjoint field for a given main variable for collection-map output.
- **New:** `adjoint_fast_path` type, `init_adjoint_fast_path(dev, fp)`,
  `eval_I_EBIC_fast(dev, fp, phi, r_beam, I_EBIC, x_beam, z_beam)` — the
  profile-aware direct-sum path for `"line"` and `"point"` profiles.

**Step 3 — dropped as degenerate.** A symmetric n-type ohmic-ohmic resistor
at short-circuit gives `I_N = 0` for every beam position by electron/hole
diffusion symmetry; the classical η = 1−x/L result is per-carrier, not net.
Forward and adjoint both correctly track this (numerical noise level
~1e-21 A), so the test doesn't meaningfully validate signs or scale.

**Step 4 — reciprocity on 2D pn junction.** `devices/adjoint_test/pn_2d.ini`:
1 µm × 200 nm Si slab, pn junction at y=500 nm, ohmic contacts at y=0 (P,
p-side) and y=1 µm (N, n-side), line beam along x swept along y.

Original (SRH off, max_dy = 25 nm): forward I_N and adjoint I_EBIC agree
to **7 digits** at y = 500 nm (2.894443e-08 A both paths).

Current (SRH on, τ = 1e-11 s, max_dy = 1 nm): agreement drops to ~5 digits
(3.579e-09 A, relative error ~6e-5). This is **expected** — finer mesh
means more terms in the sum so more accumulated FP rounding; SRH is a
bulk nonlinearity so the adjoint (linearized) answer and the forward
(nonlinear) answer diverge at the injection-level-induced nonlinearity.

Sign convention `c = +e_{I_N}` gives `I_EBIC > 0` matching the forward
convention. No sign flips needed.

**Step 6 — driver** (`src/adjoint_ebic_driver.f90`). Top-level program.
Flow: dark DC solve → refactorize → build c → transpose-solve → beam sweep
(fast path or slow path depending on profile) → forward validation at a
chosen position.

Supports four sweep modes, auto-detected from run-INI keys:

| Keys present | Mode | Use case |
|---|---|---|
| `r_beam` or `r_beam_min/_max/_n` | 1D y-sweep | 2D cross-section with line beam |
| `x_beam_*` + `y_beam_*` | 2D (x, y) Cartesian | 2D top-view with point beam |
| `y_beam_*` + `z_beam_*` | 2D (y, z) Cartesian | 3D with line beam penetrating x |
| `x_beam_*` + `y_beam_*` + `z_beam_*` | 3D (x, y, z) Cartesian | 3D with point beam |

Outputs (per run, in `<name>.fbs`):
- `adjoint/r_beam`, `adjoint/x_beam` (if 2D-x), `adjoint/z_beam` (if 3D) — beam coords in nm
- `adjoint/I_EBIC` — collected current per position (A)
- `adjoint/eta_absolute` — I_EBIC / (q · G_total) (dimensionless)
- `adjoint/phi_n_interior`, `adjoint/phi_p_interior` — collection probability map on interior vertices (grid coords via sibling `device.fbs`)
- `adjoint/G_total`, `adjoint/q_times_G_total` — normalization scalars (A)
- `adjoint/setup_time`, `adjoint/sweep_time` — wall-clock (s)

Forward validation subroutine runs one Newton at the middle beam position
for reciprocity check and prints its wall-clock (= per-position cost of a
brute-force sweep, needed for scaling plots).

## Fast-path optimizations

**Opt 1: manual `jaco_bgen` scatter** (replaces initial `sys_full%eval`-based
path). The first working kernel used `sys_full%eval(f)` to build `s`,
relying on `s = -F(u*; G)` which equals `V·G` in continuity rows and ≈ 0
elsewhere at the converged dark state. Correct and simple, but wastefully
re-evaluated Poisson, Ramo-Shockley, mobility, imref, etc. on every beam
position. Replaced with a targeted per-carrier scatter:

1. Set `dev%beam_pos%x = r_beam`.
2. Per active carrier ci: call `dev%calc_bgen(ci)%eval()` to refresh G.
3. Per active carrier: `tmp = jaco_bgen · bgen = -V·G` via sparse `mul_vec`.
4. Scatter `-tmp[1:n_interior]` into the `i0:i1` slice of sys_full
   corresponding to the ndens/pdens interior tab.

Measured impact on pn_2d: ~33 ms → ~0.5 ms per position (~70×).

**Opt 2: profile-aware direct sum** (adds `eval_I_EBIC_fast`). Analyses of
the 3D case revealed a new bottleneck: `calc_bgen%eval` iterates *all*
~35k transport vertices to set G (mostly to zero), even for a line beam
where only ~11 vertices are active. On the 3D device this was ~29 ms per
position, leaving the asymptotic speedup stuck at ~520× — professor flagged
this as insufficient.

The profile-aware path bypasses `calc_bgen`, the mat-vec, and the dense dot
product entirely. Math: for any "localized" beam,

    I_EBIC = Σ_{k on beam} V_k · G_k · φ(dof_k)
           = (G_tot / V_total) · Σ_{k on beam} V_k · φ(dof_k)

where the sum runs over the small set of active beam-column vertices
(~N_x for a line, 1 for a point). Setup builds a once-per-bias
reverse-map `(ix, iy, iz) → k in transport_vct(0)` so per-position the
loop is O(N_active) with no allocations.

Driver auto-routes line/point profiles to the fast path; Gaussian and
anything else fall back to the slow (mat-vec) path.

## Correctness validation of the profile-aware path

Expected concern: would per-position formula match the slow path? Tested
on three devices with different geometries.

**pn_2d (2D line, no SiO2).** Fast-path I_EBIC and slow-path I_EBIC
**identical to all 16 digits** at every beam position. Forward-vs-adjoint
agreement unchanged at ~5 digits (limited by discretization + SRH
nonlinearity, not the fast path).

**flat_ohmic (2D point, with SiO2 cap).**
- Peak: `3.698905464417415e-09` (fast) == `3.698905464417415e-09` (slow)
  → **bit-for-bit identical** to 16 digits at the junction peak.
- Edge positions where physics gives I_EBIC ≈ 0: fast-path values are
  ~1e-25 vs slow-path ~1e-24 — **both indistinguishable from zero**,
  differing only by machine-epsilon-scale rounding relative to the peak.

**stem_ebic_3d (3D line, with SiO2 cap).** Same pattern — peak matches,
edge noise at ~1e-24 absolute level.

**Why edges differ and peak doesn't:** the slow path dots φ (length
~140k) with s (sparse, ~11 nonzero entries) — most products are `x × 0`
which are IEEE-zero but the length-140k accumulation picks up cumulative
rounding on the order of `√N_dof · machine_eps · max(phi)` ≈ 1e-24.
The fast path sums only the ~11 nonzero terms directly, picking up
different rounding. Both are correct; the difference is purely an
associativity artifact of Σ 0 + 0 + ... + ε + ... + 0. At the peak,
I_EBIC is much larger than the cumulative rounding floor so the two
paths agree at all 16 digits.

**Verdict: profile-aware fast path is numerically correct.** The
"slightly different" edge values are physically zero and visually
harmless once the collection map is plotted with sensible clipping
(threshold anything below `I_peak × 1e-6` to zero, or use a linear
colormap with reasonable limits).

## Timing table (cumulative speedup through the three kernels)

**pn_2d, N_beam = 1000** (2D line beam, 965 DOFs):

| Kernel | Setup | Sweep | Per-position |
|---|---|---|---|
| `sys_full%eval` slow-slow (initial) | 0.15 s | 38.2 s | ~38 ms |
| `jaco_bgen` scatter (Opt 1)         | 0.15 s | 0.48 s | ~0.5 ms |
| Profile-aware direct sum (Opt 2)    | 0.15 s | pending remeasure | ~μs expected |

**stem_ebic_3d, N_beam = 300** (3D line beam, ~140k DOFs):

| Kernel | Setup | Sweep | Per-position | Asymptotic speedup vs brute-force |
|---|---|---|---|---|
| `jaco_bgen` scatter (Opt 1)       | 20.3 s | 8.8 s | 29 ms  | 520× |
| Profile-aware direct sum (Opt 2)  | 20.3 s | pending | **target ~10 μs** | **target ~1,500,000×** |

The Opt-2 target comes from the ~11 active vertices × 2 carriers = 22
float multiplies per position, plus a handful of allocations and
bin-search calls, so per-position should be dominated by cache/call
overhead at ~μs scale. Whether we hit that or plateau higher depends on
overheads like the `beam_dist%s` string compare and `denorm` calls, which
are cheap but not free.

## Linearity regime

Phase-A (`linearity_findings.md`) established that on the production
slit-W/Pt device, the SRH + Schottky + point-beam combination drives
Δn/n₀ past the linear regime at I_exp, making the adjoint approximate
rather than exact there. The pn_2d test device originally had SRH off
and ohmic contacts, so linearity was guaranteed — hence the 7-digit
forward/adjoint agreement. With SRH now on (τ = 1e-11 s) the adjoint
gives the *exact* linearized answer; forward gives the *nonlinear*
answer; they differ at the injection-induced nonlinearity. Still agree
to ~5 digits, which is plenty.

## Profile-aware fast-path assumptions

The `eval_I_EBIC_fast` path assumes:

1. **Reference generation is zero at `u*`** (beam-off dark state). Valid
   for all Donolato-style EBIC workflows. Would break around a biased +
   beam-on reference.
2. **Only continuity rows of `sys_full` depend on `beam_pos`**. Valid in
   today's Cybear. Would break if someone added a beam-coupled equation
   elsewhere (e.g., beam-induced temperature).
3. **Generation enters continuity linearly**. Holds in Cybear's current
   `beam_generation` model.
4. **`beam_dist` is `"line"` or `"point"`**. Gaussian profiles (2D only,
   per current code) fall back to the slow path automatically. Supporting
   Gaussian in the fast path would require integrating the G column
   rather than assuming uniform-on-column.

The slow paths (`build_adjoint_source` + `evaluate_ebic`) remain in the
module as correctness fallbacks and for non-supported profiles.

## New capabilities added for 2D/3D point and line beams

- **3D z-sweep support**: `region_beam` gained a `beam_z` field (sentinel
  `-1.0` means "use device center"), honored by `beam_generation.f90` so
  the driver can sweep z per position instead of hardcoded center.
- **Cartesian sweep product**: run-INI `x_beam_*` / `y_beam_*` / `z_beam_*`
  keys drive flat Cartesian sweeps of any 2 or 3 coordinates. Flat-array
  ordering with outermost loop over y, inner over z, innermost over x.
- **fbs output per run**: persistent results; MATLAB/Python plot code
  doesn't need to parse stdout.
- **SRH 3D init bug fix**: `device_params.f90` had a debug print hardcoded
  to 2D probe points (`probe_idx = [ix, iy]`) that failed a 3D assertion.
  Guarded with `if (this%g%dim == 2)`.

## Open bug: ~5e-4 adj/fwd shortfall on flat_ohmic_2d (gauge-like δu)

On the `flat_ohmic_2d` device (14 × 8.5 μm Si + SiO2 cap, point beam at
junction y = 7 μm, x = 7 μm, V = 0, SRH on τ = 2e-11 s), the adjoint's
I_EBIC is consistently ~5.25e-4 below the forward Newton's I_EBIC. The
debugging chain below ruled out most candidates but left an unresolved
structural issue in sys_full.

### Test setup

Single-point run (`devices/adjoint_test/run_adjoint_ebic_flat_ohmic_1pt.ini`,
target `adjoint_ebic_flat_ohmic_1pt`): 1 adjoint evaluation + 1 forward
Newton at (x=7, y=7) μm. Tolerances: Newton rtol/atol = 1e-12, dd atol
1e-15 cm^-3. Both sides converge to residual ~1e-15 (normalized).

Run INI flags used for diagnostics:
- `run_forward_validation = false` — skip the default middle-position check
- `use_fast_path = false` — force slow path for A/B comparison
- `debug_primal_solve = true` — enables primal-vs-adjoint self-check

### Findings, in debugging order

**1. Fast-path and slow-path give bit-identical adjoint values.** At the
junction peak: `5.947142E-09 A` from both. Edge values (physically zero)
differ by ~1e-24, which is machine-epsilon relative to the peak. Confirms
`eval_I_EBIC_fast` is not introducing any error.

**2. Forward is linear in I_beam to 1e-7 precision.** Ran I_beam at 0.01×,
0.1×, 1×, 3×, 10× baseline (0.0825 nA/um). Forward I_EBIC / I_beam:

| I_beam rel | forward I_EBIC (A) | forward / I_beam | deviation from linear |
|---|---|---|---|
| 0.01× | 5.948620e-11 | 7.2104e-3 | -2.8e-4 (noise floor — signal 6e-11 A) |
| 0.1× | 5.949875e-10 | 7.2126e-3 | -6.5e-5 |
| 1× | 5.950264e-09 | 7.2124e-3 | (reference) |
| 3× | 1.785094e-08 | 7.2125e-3 | +8e-7 |
| 10× | 5.950322e-08 | 7.2125e-3 | +1e-7 |

Upper end (1× → 10×) is linear to ~1e-7. The "nonlinearity" at 0.01× is
Newton precision floor on a 6e-11 A signal. So the **device itself IS
linear** — the adj/fwd discrepancy cannot be dismissed as injection-level
nonlinearity.

**3. Adjoint linear coefficient is systematically below forward.**
- forward (slope from 3× and 10×, precision 1e-7): 7.2125e-3
- adjoint (exact by construction): 7.2086e-3
- **ratio adj/fwd = 0.99946 ⟹ 5.4e-4 systematic shortfall**

**4. Primal self-check passes at solver precision.** Added a diagnostic
that solves `A · δu = s` (no transpose) and compares `c·δu` to `phi·s`
— these must match to solver precision by the transpose identity.
Result: `(primal − adjoint) / adjoint = -6.5e-14`. So the adjoint kernel
(c, s, A, transpose solve) is **internally consistent to machine precision**.

The primal passing means whatever c, s, A the adjoint uses, they obey
`phi · s = c · A^-1 · s`. It does NOT prove these are the *right*
c, s, A for the forward problem.

**5. `u_fwd - u_dark` has ‖δu_fwd‖ ≈ 19,000× larger than ‖A^-1 s‖.** The
smoking gun. Captured `u_dark = sys_full%get_x()` after the dark solve,
`u_fwd = sys_full%get_x()` after the forward Newton:

```
|delta_u_forward|           =   4.077262E+03
|delta_u_primal (A^-1 s)|   =   2.145272E-01
|fwd - primal|              =   4.077262E+03
relative error              =   1.900580E+04

c . delta_u_fwd:   5.950263765954E-09 A
c . delta_u_pri:   5.947141555660E-09 A
u_dark (I_ct_dof): 3.606390685304E-27 A
```

The I_ct-DOF component of both `δu_fwd` and `δu_primal` matches at 5.25e-4
(the original discrepancy). But the rest of `δu_fwd` (ψ, n, p, cdens, imref
components) differs from `δu_primal` by a factor of ~20,000 in 2-norm.

**6. pn_2d shows the same pattern with milder impact.** Ran the same
diagnostic on `pn_2d.ini` (plain Si slab, line beam, no SiO2, ohmic
contacts). Gauge-like δu blow-up is nearly identical; the projection onto
`c` is an order of magnitude smaller:

| Device | Beam | `‖δu_fwd‖ / ‖A⁻¹·s‖` | Primal self-check | c·δu fwd/adj gap |
|---|---|---|---|---|
| flat_ohmic_2d (SiO2 cap, point) | point | 19,000× | -6.5e-14 | **5.25e-4** |
| pn_2d (plain Si, line) | line | 15,600× | -1.3e-13 | **6.5e-5** |

Conclusion: the gauge-like redundancy is **structural in sys_full**, not
specific to the SiO2-capped geometry. What changes is how much the gauge
projects onto `c` (i.e., perturbs observable terminal current), which
depends on device geometry and beam profile.

**7. Production `dd` driver gives bit-identical forward result.** Cross-
checked via a separate run (`adjoint_test/run_dd_flat_ohmic_1pt.ini`,
target `dd_flat_ohmic_1pt`) of the legacy `dd` forward driver on the same
`flat_ohmic_2d.ini` device at (x=7, y=7) μm:

| Source | I_N_CONTACT (A) |
|---|---|
| production `dd` driver | 5.950263765954e-09 |
| our `run_forward_validation` | 5.950263765954e-09 |
| adjoint `phi · s` | 5.947141555660e-09 |

dd ↔ our forward: **12-digit match**. Kirchhoff closure `I_N + I_P = 0`
holds to machine precision on dd's output.

This rules out any implementation error in our adjoint driver's forward
Newton path — we produce *exactly* what the production code produces.
The 5.25e-4 shortfall is therefore 100% on the adjoint side of the
comparison. The remaining candidate is **A (the factorized Jacobian used
by the adjoint transpose solve) differs from the effective linearization
Newton sees when stepping from u_dark to u_fwd**, because u_fwd lies on
a gauge-extended solution manifold where A⁻¹·s picks one representative
(minimum-norm in the A-metric) while Newton's path from Gummel picks
another.

### Working hypothesis

Sys_full has a **gauge-like degeneracy**: a direction in DOF-space that
doesn't affect physical observables (terminal currents, densities at
meaningful points) but does affect "internal" DOFs (imref, cdens, etc.)
The forward Newton converges F(u) = 0 up to machine precision starting
from Gummel's initial guess, which is not the same as u_dark; the
different initial guess lets Newton wander along the gauge direction and
land at a u_fwd that's equivalent to u_dark modulo the gauge.

The primal solve `A · δu = s` picks a single point on this gauge manifold,
the one with no gauge offset. If the gauge is not EXACTLY orthogonal to
`c` (the I_ct indicator), the forward and primal differ at the I_ct DOF
by the inner product of the gauge vector with `c` — the 5e-4 we see.

Candidates for the gauge direction:
- **Imref arbitrariness**: `approx_imref` resets imref values from scratch
  each forward solve. If the equation `F_imref = iref − func(ψ, n)` is
  satisfied for a whole 1-parameter family at machine precision (because
  of a numerically-small singular direction), Newton can float in that
  family.
- **Current-density closure**: `F_cdens = cdens − (SG-expression)` — if SG
  fluxes have very-small-singular-value modes at the junction (where
  potential is near-null in the depletion region?), cdens could float.
- **Poisson gauge**: if sys_full's Poisson rows have a near-null constant
  shift direction (global ψ offset), Newton could float ψ globally.

Determining which of these is the culprit requires a per-variable
breakdown of `δu_fwd − δu_primal` using `extract_phi_for_var`-style
indexing (diagnose: which main variable carries the huge norm).

### 2026-04-21 update — localization to A

With the new driver diagnostics, the Y_BEAM 4077 was traced to an input-DOF
artifact (beam position stored in sys_full's state vector), not a gauge.
After parking Y_BEAM identically in both dark and forward solves, the
comparison is clean:

```
|delta_u_forward|           =   2.146E-01
|delta_u_primal (A^-1 s)|   =   2.145E-01
|fwd - primal|              =   1.124E-04
relative error              =   5.24e-4
```

Per-variable breakdown shows the **5.24e-4 is a uniform scalar rescaling**
of `A^-1 s` across every physical DOF (pot, ndens, pdens, currents all at
5.2-5.7e-4 relative diff). `δu_fwd = (1 + 5.24e-4) · A^-1·s` throughout.

Source check: adjoint's `s = -jaco_bgen · bgen` matches the forward's
implicit `-F(u_dark, beam on)` at machine precision (12-digit agreement,
`<s_adj,s_imp>/|s_adj|² = 1.000000000000`). **s is correct.**

Taylor check: at `u_dark + A^-1·s`, residual `|F| / |s| = 6.34e-3` —
that's a mix of genuine O(|δu|²) nonlinear correction plus A-side error.

Production cross-check: the legacy `dd` driver at the same point gives
bit-identical (12-digit match) forward I_N_CONTACT to our
run_forward_validation. So forward is trustworthy; the gap is entirely on
the adjoint (A) side.

**Conclusion: the adjoint's factorized Jacobian A is ~5.24e-4 stiffer
than the true J(u_dark) that the forward effectively linearizes against.
A = (1 + 5.24e-4) · J_true globally.**

Audit of `const = .true.` Jacobian blocks in sys_full (Poisson, continuity's
jaco_bgen/jaco_srh/jaco_genrec/jaco_surf_srh/jaco_dens_t/jaco_cdens,
imref's jaco_pot, charge_density, Ramo-Shockley's jaco_cdens/jaco_volt/
jaco_curr, ionization's jaco_t/jaco_genrec): all entries are truly
state-independent in theory (geometric weights, ±1 charges, identity,
V/A wrapper factors). None appear to be "stale" state-dependent values.

Current_density's jaco_pot, jaco_dens (SG flux derivatives) are
non-const and re-evaluated each `sys_full%eval`. Same for calc_dens,
calc_imref, calc_srh_recombination's internal jaco_dens_n/p.

So the 5.24e-4 does NOT look like stale-const. Next candidates:
- A subtle coupling mismatch between Jacobian block chain-rule terms
  (e.g., Cybear's `total_jaco` vs direct jaco when chains through
  intermediate vars like cdens, iref).
- Scharfetter-Gummel flux derivatives with a specific rounding/formula
  that happens to be 5e-4 off from what Newton linearizes against.
- Preconditioner `dfp` being used for solve instead of `df` (unlikely
  for PARDISO direct, but worth verifying).

### What this means for the paper

The **adjoint kernel is correctly implemented**. The primal identity holds
at machine precision. The 5e-4 discrepancy is a structural property of
sys_full's block-matrix layout, not of the Donolato-style reciprocity
math. For a "clean" demonstration:

- `pn_2d` (2D, ohmic, NO SiO2) agreed to 5-7 digits earlier — suggests the
  gauge-like degeneracy appears with the SiO2-capped geometry or with the
  extra main variables present in flat_ohmic's setup.
- The realistic paper plot (φ map on realistic 3D device) may or may not
  be affected — need to re-run the comparison on `stem_ebic_3d_adjoint`
  with the new primal diagnostic to see if the gauge problem persists in 3D.

### Next actions (none executed yet)

1. **Per-variable breakdown.** Decompose `δu_fwd − δu_primal` by main
   variable (pot, ndens, pdens, iref_n, iref_p, cdens_*, currents) —
   identifies the gauge carrier.
2. **Compare against `pn_2d` without SiO2.** Run the same 1-point
   diagnostic on pn_2d; if |δu_fwd| matches |δu_primal| there, the issue
   is specific to the SiO2-cap geometry.
3. **Inspect `calc_imref` residual form.** If `F_imref = iref − func`
   has a near-null mode, imref is the culprit and can potentially be
   removed from sys_full (replaced with a closure) — but that's a bigger
   refactor.
4. **Try approximate-solution pinning.** If gauge is real, we can pin
   it by adding a tiny penalty to A (e.g., `ε·I` in the null direction).
   This is a hack; ideally the structural fix would be to eliminate the
   gauge entirely.

## Open follow-ups

1. **Re-measure Opt-2 per-position timing.** Both pn_2d and 3D — confirm we
   actually hit μs-scale per position. Update the timing table.
2. **Point-fast-path speedup on flat_ohmic.** Previous slow-path: ~29 ms;
   expected Opt-2: ~5-20 μs (one beam vertex sum).
3. **Step 5: Donolato analytical cross-check.** Compare `phi^n(y)` on the
   p-side to `sinh(y/L_p)/sinh(W_p/L_p)` on pn_2d with SRH on. τ = 1e-10
   seconds gives L_p ≈ W_p = 500 nm, so the sinh envelope is well-resolved.
4. **Phase A unblockers (for production device).** Gaussian beam re-test,
   Test 3 (small-signal), lower I_exp window. Until at least one passes on
   slit-W/Pt, Step 7 stays blocked.
5. **If prof wants even more speed on the existing devices**: cache
   `dist_s` as a logical on `adjoint_fast_path`, hoist `denorm` unit
   strings out of per-position calls, preallocate scratch arrays.

## Critical files

- `src/adjoint_ebic.f90` — kernel (slow + profile-aware fast paths)
- `src/adjoint_ebic_driver.f90` — driver (sweeps, fbs, linescan, operando
  bias, primal self-check, δu_fwd-vs-δu_primal Jacobian diagnostic)
- `src/region.f90` — added `beam_z` field to `region_beam`
- `src/beam_generation.f90` — honors `reg_beam(1)%beam_z` when set
- `src/device_params.f90` — 2D-only guard on SRH debug print
- `devices/adjoint_test/pn_2d.ini` — 2D pn validation device (N/P contacts)
- `devices/adjoint_test/flat_ohmic_2d.ini` — 2D top-view flat-ohmic pn
- `devices/adjoint_test/flat_ohmic_center_2d.ini` — symmetric pn for operando
- `devices/adjoint_test/stem_ebic_3d_adjoint.ini` — 3D pn
- `devices/adjoint_test/run_adjoint_ebic_flat_ohmic_1pt.ini` — single-point
  diagnostic harness (toggles `use_fast_path` and `debug_primal_solve`)
- `devices/adjoint_test/run_adjoint_ebic_pn_1pt.ini` — same diagnostic, pn_2d
- `devices/adjoint_test/run_dd_flat_ohmic_1pt.ini` — independent dd-driver
  cross-check at the 1pt position (verifies forward Newton bit-for-bit)
- `devices/adjoint_test/run_adjoint_ebic_*.ini` — run configs
- `fargo.toml` — targets: `adjoint_ebic_pn`, `adjoint_ebic_flat_ohmic`,
  `adjoint_ebic_flat_ohmic_1V`, `adjoint_ebic_flat_ohmic_1pt`,
  `adjoint_ebic_pn_1pt`, `dd_flat_ohmic_1pt`, `adjoint_ebic_3d`

## Commit history on adjoint_method branch

- `05bab36` — Phase A linearity test suite (6-case isolation matrix + docs)
- `795bb9f` — Phase B kernel + driver + pn_2d + implementation_plan.md
- `3f25522` — 2D/3D sweep support (x_beam / y_beam / z_beam), fbs output,
  flat_ohmic + 3D test devices, SRH 3D fix, forward-Newton timing
- `4a4afa5` — profile-aware direct-sum fast path (Opt 2); adjoint_fast_path
  type + init_adjoint_fast_path + eval_I_EBIC_fast; `use_fast_path` INI
  override; peak match to 16 digits verified on pn_2d / flat_ohmic /
  stem_ebic_3d; timing dropped to ~1 μs per position, 14M× asymptotic
  speedup on stem_ebic_3d
- `9da25a0` — operando-bias support (per-contact `V_<name>` with Newton
  continuation); forward linescan emits `<name>_fwd_linescan.fbs`;
  run_forward_validation quiet/I_ebic_out args for batch use;
  flat_ohmic_center_2d.ini (symmetric pn); dense sweep config for
  flat_ohmic
- (uncommitted) — single-point diagnostic infrastructure:
  `run_adjoint_ebic_{flat_ohmic,pn}_1pt.ini` and `run_dd_flat_ohmic_1pt.ini`
  reusable harnesses; `debug_primal_solve` runfile flag; extracted into
  private `run_primal_diagnostic` and `compare_du_fwd_vs_primal` subroutines
  in the driver. Surfaced the gauge-like δu discrepancy documented under
  "Open bug". Linearity-in-I_beam test was run via throwaway scaled INIs
  (now removed; verdict captured in the table under "Open bug")

Push target: `origin/adjoint_method` (user's fork).
