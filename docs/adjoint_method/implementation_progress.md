# Adjoint EBIC — Implementation Progress

Empirical status and findings for Phase B of the adjoint-method study.
Companion to `implementation_plan.md` (forward-looking design) and
`linearity_findings.md` (Phase A empirical findings). Last updated
2026-04-20.

## Status summary

| Step | Description | Status |
|------|-------------|--------|
| 1 | `keep_factorization` in fortran-basic `steady_state.f90` | **SKIPPED** (not needed — see below) |
| 2 | `src/adjoint_ebic.f90` kernel | **DONE** — validated |
| 3 | 1D resistor η = 1−x/L validation | **DROPPED** (degenerate setup) |
| 4 | pn reciprocity harness at 5 r_beam | **DONE** — 7-digit agreement (no-SRH mesh); 5-digit on dense-mesh + SRH |
| 5 | Donolato analytical cross-check | Not started |
| 6 | `src/adjoint_ebic_driver.f90` driver | **DONE** — 1D / 2D (x,y) / 2D (y,z) / 3D sweeps, fbs output |
| 7 | Production slit-W/Pt collection map | Blocked on Phase A |
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
- `src/adjoint_ebic_driver.f90` — driver (1D / 2D / 3D sweeps, fbs output)
- `src/region.f90` — added `beam_z` field to `region_beam`
- `src/beam_generation.f90` — honors `reg_beam(1)%beam_z` when set
- `src/device_params.f90` — 2D-only guard on SRH debug print
- `devices/adjoint_test/pn_2d.ini` — 2D pn validation device (N/P contacts)
- `devices/adjoint_test/flat_ohmic_2d.ini` — 2D top-view flat-ohmic pn
- `devices/adjoint_test/stem_ebic_3d_adjoint.ini` — 3D pn (doping + I_beam fixed)
- `devices/adjoint_test/run_adjoint_ebic_*.ini` — run configs
- `fargo.toml` — `adjoint_ebic_pn`, `adjoint_ebic_flat_ohmic`, `adjoint_ebic_3d` targets

## Commit history on adjoint_method branch

- `05bab36` — Phase A linearity test suite (6-case isolation matrix + docs)
- `795bb9f` — Phase B kernel + driver + pn_2d + implementation_plan.md
- `3f25522` — 2D/3D sweep support (x_beam / y_beam / z_beam), fbs output,
  flat_ohmic + 3D test devices, SRH 3D fix, forward-Newton timing
- (uncommitted) — profile-aware direct-sum fast path; adjoint_fast_path
  type + init_adjoint_fast_path + eval_I_EBIC_fast; updated progress doc

Push target: `origin/adjoint_method` (user's fork).
