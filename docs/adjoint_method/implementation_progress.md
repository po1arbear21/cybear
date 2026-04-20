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
| 4 | pn reciprocity harness at 5 r_beam | **DONE** — 7-digit agreement |
| 5 | Donolato analytical cross-check | Not started |
| 6 | `src/adjoint_ebic_driver.f90` driver | **DONE** |
| 7 | Production slit-W/Pt collection map | Blocked on Phase A |
| — | Fast-path optimization (manual jaco_bgen scatter) | **DONE**, awaiting re-validation |

## What was done

**Step 1 dropped as no-op.** The `keep_factorization` flag we planned to add
to `steady_state_m` turned out to be unnecessary. PARDISO's factorization
stored inside `ss%solver` naturally persists as long as the `ss` object
itself is alive — nothing deletes it except explicit `ss%destruct()` calls,
which our driver simply doesn't make until after the adjoint work is done.
No fortran-basic edits required.

**Step 2 — kernel** (`src/adjoint_ebic.f90`). Five public subroutines:

- `build_collection_functional(dev, contact_id, c)` — builds `c` as a unit
  indicator at the sys_full DOF holding the target contact's Ramo-Shockley
  terminal current. Uses `search_main_var("currents")` + `res2block`.
- `refactorize_at_converged(dev, solver)` — re-evaluates J at the converged
  Newton state and re-factorizes (pins LU to `u*` exactly, vs. the
  `u_{final-1}` state Newton's last factorize was done at).
- `build_adjoint_source(dev, r_beam, s)` — computes `s = V·G(r;r_beam)`
  scattered into continuity interior rows.
- `solve_adjoint(solver, c, phi)` — one-liner transpose-solve via
  `trans='T'`.
- `evaluate_ebic(phi, s)` — scalar dot product.
- `extract_phi_for_var(dev, phi, var_name, tab_idx, phi_slice)` — slices
  the adjoint field for a given main variable (e.g., "ndens" interior) for
  collection-map output.

**Step 3 — dropped as degenerate.** A symmetric n-type ohmic-ohmic resistor
at short-circuit gives `I_SRC = 0` for every beam position by electron/hole
diffusion symmetry; the classical η = 1−x/L result is per-carrier, not net.
Forward and adjoint both correctly track this (numerical noise level
~1e-21 A), so the test doesn't meaningfully validate signs or scale.

**Step 4 — reciprocity on 2D pn junction.** `devices/adjoint_test/pn_2d.ini`:
1 µm × 200 nm Si slab, pn junction at y=500 nm, ohmic contacts at y=0 (P,
p-side) and y=1 µm (N, n-side), line beam along x swept along y.

Five-position test result (forward Newton vs. adjoint at matched bias):
- Forward `I_N` at y = 500 nm: `2.894443e-08 A`
- Adjoint `I_EBIC` at y = 500 nm: `2.894443e-08 A`
- **Agreement: 7 digits** (rel. error ≈ 1e-7).

Beam-sweep shape matches classical pn-diode EBIC:

| y [nm] | region | I_EBIC [A] | mechanism |
|---|---|---|---|
| 100 | deep-p | 6.2e-9 | electron diffusion across p quasi-neutral |
| 400 | p-edge of depletion | 2.5e-8 | electron near junction, mostly collected |
| 500 | junction center | 2.9e-8 (peak) | both carriers separated by field |
| 600 | n-edge of depletion | 2.8e-8 | holes near junction, mostly collected |
| 900 | deep-n | 7.1e-9 | hole diffusion across n quasi-neutral |

Device is deliberately short (1 µm) and SRH-off, so the profile is
geometry-limited (diffusion to contacts without recombination). Peak at
junction ≈ q · G_tot → collection efficiency η ≈ 1 (every generated pair
collected).

Sign convention `c = +e_{I_N}` gives `I_EBIC > 0` matching the forward
convention. No sign flips needed.

**Step 6 — driver** (`src/adjoint_ebic_driver.f90`). Top-level program
running: dark DC solve → refactorize → build c → transpose-solve → beam
sweep → forward validation at a chosen position. Accepts `r_beam` as either
an explicit list or `r_beam_min` / `r_beam_max` / `r_beam_n` linspace.
Prints I_EBIC + η_absolute + η_relative for each position plus setup/sweep
timing.

## Fast-path optimization

The first working kernel used `sys_full%eval(f)` to build `s`, relying on
`s = -F(u*; G)` which equals `V·G` in continuity rows and ≈ 0 elsewhere at
the converged dark state. Correct and simple, but wastefully re-evaluates
Poisson, Ramo-Shockley, mobility, imref, etc. on every beam position even
though none of those depend on `beam_pos`.

Replaced with a targeted scatter:

1. Set `dev%beam_pos%x = r_beam`.
2. Per active carrier ci: call `dev%calc_bgen(ci)%eval()` to refresh G.
3. Per active carrier: compute `tmp = jaco_bgen · bgen = -V·G` on interior
   rows, 0 on contact rows (via sparse `mul_vec`).
4. Scatter `-tmp[1:n_interior]` into the `i0:i1` slice of sys_full
   corresponding to the ndens/pdens interior tab.

Mathematically identical to the eval-based path (both produce the same `s`
up to Newton residual floor, which in the fast path is zero by construction
in non-continuity rows). Expected ~10× per-position speedup.

## Timing on pn_2d (slow path, before fast-path swap)

| N_beam | Setup | Sweep | Per-position | Sweep/setup |
|---|---|---|---|---|
| 101  | 0.25 s | 3.34 s | ~33 ms | 13 |
| 1000 | 0.15 s | 38.2 s | ~38 ms | 256 |

Setup cost (Newton + refactorize + adjoint transpose-solve) is ~150 ms,
dominated by Newton. Per-position cost was ~35 ms, dominated by the
wasteful parts of `sys_full%eval`. After the fast-path swap, per-position
should drop to ~2-5 ms.

Speedup vs. brute-force forward (estimated: one Newton ~150 ms per beam
position):

| Regime | Brute-force | Adjoint (slow) | Adjoint (fast, proj.) |
|---|---|---|---|
| 1000 positions | ~150 s | ~38 s | ~3-5 s |
| Asymptotic ratio (T_Newton / T_per_pos) | 1× | ~4× | ~30-50× |

On this small 965-DOF device the setup is cheap, so the adjoint advantage
is modest. On a realistic 3D device with ~1M DOFs, factorization dominates
Newton; the expected adjoint advantage there is 100-1000× depending on
how many beam positions are needed.

## Linearity regime

Phase-A (`linearity_findings.md`) established that on the production
slit-W/Pt device, the SRH + Schottky + point-beam combination drives
Δn/n₀ past the linear regime at I_exp, making the adjoint approximate
rather than exact there. The pn_2d test device has SRH off and ohmic
contacts, so linearity is guaranteed — which is why we see 7-digit
forward/adjoint agreement. When SRH is re-enabled or Schottky contacts
are used, we'd expect agreement to degrade to the Δu/u₀ order for
high-injection beam positions. This is the physics, not a bug: the adjoint
gives the *exact* answer for the linearized problem, which just isn't the
same as the nonlinear forward at high injection.

## Caveats on the fast-path

The fast path assumes:

1. **Reference generation is zero at `u*`**. Valid for all Donolato-style
   EBIC workflows (we always linearize around a beam-off state). Would
   break if someone set up an adjoint around a biased + beam-on reference.
2. **Only continuity rows of `sys_full` depend on `beam_pos`**. Valid in
   today's Cybear; would break if someone added a beam-coupled equation
   elsewhere (e.g., beam-induced temperature).
3. **Generation enters continuity linearly**. Holds in Cybear's current
   `beam_generation` model. A nonlinear beam coupling would require
   returning to the full-eval path.

The eval-based path remains in git history as a safe fallback.

## Open follow-ups

1. **Re-run pn_2d with the fast path.** Verify 7+ digit reciprocity match
   (should in fact improve to ~9 digits since fast-path `s` has no Newton
   noise floor). Measure actual per-position timing.
2. **Step 5: Donolato analytical cross-check.** Compare `phi^n(y)` on the
   p-side to the closed form `sinh((y+W_p)/L_p)/sinh(W_p/L_p)`. Requires
   turning SRH on in `pn_2d.ini` (otherwise L_p → ∞ and the formula reduces
   to the trivial linear-ramp that we already see). Needs the user's call
   on τ_p values.
3. **SRH inheritance test.** Flip `srh = true` in `pn_2d.ini` with a
   moderate τ. Confirm forward/adjoint still agree to ~7 digits. Visually
   confirm `phi^n(y)` develops exponential envelope (instead of linear
   ramp) with decay length ~L_n.
4. **Phase A unblockers (for production device).** Gaussian beam re-test,
   Test 3 (small-signal), lower I_exp window. Until at least one of these
   passes on slit-W/Pt, Step 7 stays blocked.

## Critical files

- `src/adjoint_ebic.f90` — kernel (fast-path implementation in
  `build_adjoint_source`)
- `src/adjoint_ebic_driver.f90` — driver
- `devices/adjoint_test/pn_2d.ini` — 2D pn validation device
- `devices/adjoint_test/run_adjoint_ebic_pn.ini` — run config
- `fargo.toml` — `adjoint_ebic_pn` target

## Commit history on adjoint_method branch

- `05bab36` — Phase A linearity test suite (6-case isolation matrix + docs)
- `795bb9f` — Phase B kernel + driver + pn_2d + implementation_plan.md
- (uncommitted) — fast-path optimization in `build_adjoint_source`

Push target: `origin/adjoint_method` (user's fork).
