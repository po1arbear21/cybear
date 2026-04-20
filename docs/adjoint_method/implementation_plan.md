# Adjoint-method EBIC — Implementation Plan

Design doc for replacing Cybear's "one forward solve per beam position" EBIC
sweep with "one adjoint solve per bias point" using Donolato-style
reciprocity. Sibling docs: `linearity_test.md` (Phase-A spec),
`linearity_findings.md` (Phase-A results).

## 1. Context

For a linearized drift-diffusion system `A · δu = s(r_beam)` around the dark
short-circuit state `u* = (ψ*, n*, p*)`, the EBIC terminal current at
collecting contact k is

    I_EBIC(r_beam) = c_k⊤ A⁻¹ s(r_beam)
                   = (A⁻⊤ c_k)⊤ s(r_beam)
                   = φ_k⊤ s(r_beam).

Solve `A⊤ φ_k = c_k` **once** per bias point. Every beam position is then a
dot product between the stored adjoint `φ_k` and the per-position source
`s(r_beam)`. The scaling advantage is enormous: one LU factorization
amortized over thousands of beam positions replaces thousands of independent
nonlinear solves.

Equivalently in volume-integral form — standard in the adjoint-EBIC
literature — the terminal current perturbation is

    δI_k(r_beam) ≈ −e Σ_i tr_vol_i [ν*_{k,i} + κ*_{k,i}] G_E(r_i; r_beam)

where `(ν*, κ*)` are the electron- and hole-continuity components of `φ_k`.
The ψ-component drops out because the beam generation only appears in
continuity rows, not in Poisson.

The headline deliverable is the spatial collection-probability map
`φⁿ(r) + φᵖ(r)` — a single figure per bias point that encodes EBIC response
for every beam position.

## 2. Dependency on Phase A

The adjoint formulation is exact only when the forward map `F: G → I_EBIC`
is linear around `u*`. Phase A (linearity validation) addressed this
empirically — results in `linearity_findings.md`. Headline:

> Linearity fails iff **SRH + Schottky + point beam** are all present on
> the top-view geometry. Line-beam variants are linear regardless.

The production slit-W/Pt geometry (point beam + Schottky + SRH) is
**non-linear at I_exp as currently simulated**. Three open unblockers
before the adjoint result on that device is physically meaningful:

1. Re-run with `beam_dist = "gaussian"` — point beam is a pessimistic
   upper bound on Δn/n₀.
2. Implement Phase-A Test 3 to directly measure `max(Δn/n₀)` and
   `max(Δψ/V_T)` at I_exp — gives a clean yes/no.
3. Lower-I_exp injection-level window if the experiment can tolerate it.

**Phase B proceeds regardless** on simpler geometries (2D-extruded pn
junction) because:
- The adjoint infrastructure is device-agnostic — it's wired once and
  works for every device that builds a `sys_full` Jacobian.
- The 1D cases have analytical solutions, so signs and prefactors are
  pinned *independently* of the forward solver (which could share a sign
  bug).
- If Phase-A unblockers succeed, the adjoint is drop-in for the
  production device.
- If they fail, the 1D results still stand as a methods paper separate
  from the slit-W/Pt SISPAD draft.

## 3. Inherited infrastructure

Each finding verified against the current source on the `adjoint_method`
branch.

### F1. Transpose solve is already supported
- `lib/fortran-basic/src/util/solver/solver_base.f90:101, 253` — abstract
  `solve_` / `solve_vec` take optional `trans ∈ {'N','T','C'}`.
- `pardiso.f90:254` — sets `iparm(12) = 2` for transpose, reuses the
  existing LU.
- `mumps.f90:239` — sets `ICNTL(9) = 2` for transpose, same.

`φ = A⁻⊤ c` is `solver%solve(c, phi, trans='T')`. No algorithmic work.

### F2. Jacobian does NOT persist past Newton convergence (to be fixed)
- `steady_state.f90:148` — `this%solver%destruct()` wipes the LU.
- The `dfdx` pointer is local to the Newton loop (`steady_state.f90:314`)
  and re-evaluated every iteration (`:428`).

Fix: add a `keep_factorization` option + accessor (Step 1).

### F3. Terminal current `I_ct` is a DOF of `sys_full`
- `device.f90:62` — `type(current), allocatable :: curr(:)`.
- `device.f90:349` — Ramo-Shockley equation `ramo_curr` attached to
  `sys_full`.
- `ramo_shockley.f90:224` — `I_ct` is the main variable of that residual.
- `ramo_shockley.f90:253` — `jaco_cdens` carries
  `curr_fact · surf · (ramo_shape(idx2) − ramo_shape(idx1))`
  (Ramo-Shockley weighting, displacement-current aware).
- `ramo_shockley.f90:265` — `jaco_curr` is identity.

**Structural consequence:** the collection functional `c_k` is a unit
indicator vector at the DOF offset of `dev%curr(k)%x`. No hand-assembly
of ∂I/∂u from edge-flux derivatives. See appendix for the comparison
against a Dirichlet-probe `c`.

### F4. Source `s` builder is already in place
- `continuity.f90:330` — `F = dn/dt · V + div(j) · V − G_beam · V = 0`
  implies `dF/dG = −V`.
- `continuity.f90:334` — writes `−tr_vol(idx)` into `jaco_bgen`.
- `continuity.f90:416` — `jaco_bgen%matr%mul_vec(bgen%get(), ...)` gives
  the per-row source contribution for free.

The linearized RHS is `s = −J_bgen · G(r; r_beam)` scattered into n- and
p-continuity rows. Sign convention is dictated by the existing forward
EBIC path — we match it.

### F5. SG edge-flux derivatives are available but not needed
- `current_integral.f90:75–92` returns `djdn(2)` and `djddpot`.

Only relevant if `c` were assembled from a non-Ramo-Shockley contact
definition — not the chosen path, but available for cross-checks.

## 4. Workflow (per bias point)

Three steps; steps 1 and 2 run once per bias, step 3 per beam position.

1. **Nonlinear DC solve at target bias.** Gummel warm-up → Full Newton.
   Save `J` and its factorization at the converged state via
   `keep_factorization`.
2. **Per collecting contact k** (usually one, see D2):
   build c_k via `build_collection_functional`, solve `J⊤ x_k* = c_k`.
   One triangular sweep, no re-factorization (transpose flag).
3. **Per beam position r_beam:** form `G_E(r; r_beam)`, evaluate
   `δI_k = −e Σ_i tr_vol_i [ν*_{k,i} + κ*_{k,i}] G_E_i`. Equivalent to
   `φ_k⊤ s(r_beam)` with `s = −jaco_bgen · G`.

The matrix form and the volume-integral form are the same equation —
different bookkeeping.

## 5. Implementation order

### Step 1 — Jacobian persistence (fortran-basic)
`lib/fortran-basic/src/analysis/steady_state.f90` — add
`keep_factorization` option. When set, skip `solver%destruct()` on
teardown. Expose `solver` and the `dfdx` pointer via accessors. This is
the **only** edit outside `cybear/src/`.

### Step 2 — Adjoint kernel (cybear)
`src/adjoint_ebic.f90` (new):
- `build_collection_functional(dev, contact_id) → c` — unit indicator at
  `dev%curr(contact_id)%x` DOF in `sys_full`.
- `build_adjoint_source(dev, r_beam) → s` — set `beam_pos`, refresh `bgen`
  via `beam_generation%eval`, apply `jaco_bgen%matr%mul_vec`.
- `solve_adjoint(solver, c) → φ` — one line, `trans='T'`.
- `evaluate_ebic(φ, s) → I_EBIC` — scalar dot product.
- `extract_collection_map(φ, dev) → (φⁿ + φᵖ)(r)` — spatial field for
  output via existing `c_vars` plumbing.

### Step 3 — (dropped)
A uniform ohmic-ohmic resistor at short-circuit is degenerate for
forward-vs-adjoint validation: electron and hole diffusion give equal
and opposite contributions to the terminal current, so net I_EBIC = 0
for every beam position by symmetry. The η = 1 − x/L identity is a
per-carrier diffusion probability, not the net current that the forward
or adjoint solver computes. Validation collapses into Step 4 on the pn
junction, which has non-degenerate I_EBIC and supports direct
forward-vs-adjoint reciprocity.

### Step 4 — pn reciprocity harness
`src/adjoint_ebic_test.f90` (new). `devices/adjoint_test/pn_2d.ini` at 5
beam positions: deep-p, depletion p-edge, depletion center, depletion n-edge,
deep-n.
Forward path via existing EBIC flow; adjoint via Step 2 kernel. Report
per-point relative error and pass/fail at tolerance D5.

### Step 5 — Donolato analytical cross-check
Compare `φⁿ(x) + φᵖ(x)` on 1D pn to
`sinh((x + W_p) / L_p) / sinh(W_{qn,p} / L_p)` on p-side (symmetric on
n-side). Expected match within discretization error ~mesh/L_{n,p}.
Methods-section figure, not a pass/fail gate.

### Step 6 — Driver wire-up
New top-level program `src/adjoint_ebic_driver.f90` (or similar) rather
than extending `stem_ebic.f90`. Shared code with stem_ebic is `dev`
initialization and `beam_generation` — both cleanly exportable. Fargo
target: `adjoint_ebic` or `adjoint_ebic_test`.

### Step 7 — Conditional: production device
Only if Phase-A unblockers succeed. Extend the driver to dump the full
`(φⁿ + φᵖ)(r)` collection map on
`devices/stem_ebic/stem_ebic_slit_pt.ini`. SISPAD-draft headline plot.

## 6. Verification ladder

Each step gates the next.

1. **pn reciprocity.** Forward vs adjoint I_EBIC agree at 5 r_beam within
   D5 tolerance. This is the gateway validation — pins signs/prefactors
   and catches any bug in the kernel, driver, or Jacobian-reuse path.
   If this fails, no downstream result is trustworthy.
2. **pn Donolato analytical cross-check.** Discretization-error bound
   for the methods section.
3. **Low-injection audit (D6).** From one forward solution at r_beam in
   depletion, compute `max(δn/n*)`, `max(δp/p*)`, `max(δψ/V_T)`. Report.
4. **Production slit-W/Pt collection map.** No analytical truth;
   sanity-check against ~5 forward spot-checks at the Step 1 tolerance.

## 7. Design decisions

**D1. Validation device ordering: pn junction → (conditional) slit-W/Pt.**
The 2D-extruded pn diode at `devices/adjoint_test/pn_2d.ini` gives
non-degenerate I_EBIC and supports direct forward-vs-adjoint reciprocity
— this pins signs/prefactors and catches kernel bugs. Slit deferred to
Step 7.

**D2. Collecting contact: single configurable `contact_id`.** On a
two-contact device Kirchhoff gives `I_SRC = −I_DRN`, so one adjoint per
run suffices. Driver default: first `ohmic` or `schottky` contact in
the device INI, or user-specified.

**D3. Output.** Both forms every run:
- Scalar `I_EBIC(r_beam_k)` at a prescribed beam-position list — direct
  diff against the existing forward EBIC output.
- Spatial field `φⁿ(r) + φᵖ(r)` — the collection-probability map.

**D4. Sign convention.** Pinned empirically by Step 4: adjoint I_EBIC at
a fixed r_beam must match the sign and magnitude of the forward I_EBIC
at the same position. Placeholder proposal: `c = +e_{I_collect}`.
Confirm or flip on first reciprocity run.

**D5. Reciprocity tolerance: target < 1e-10, report actual.** Default
solver settings first; if we land at 1e-9, report 1e-9 rather than
tighten to meet a round number.

**D6. Low-injection check: post-hoc only.** No runtime assertion;
`max(δn/n*)` in depletion computed and reported for the methods
section.

## 8. Appendix: Ramo-Shockley indicator vs Dirichlet probe

Two ways to build the collection functional `c_k`:

**(A) Classical Donolato / Dirichlet probe:** set a +1 test BC at
contact k, zero elsewhere. Assemble `b_k` as the corresponding RHS
contribution.

**(B) Ramo-Shockley indicator (chosen):** `c_k = e_{I_ct(k)}` — a unit
vector at the current-variable DOF index of `dev%curr(k)%x`.

Advantages of (B):

- **Displacement-current aware.** The Ramo-Shockley row in the forward
  Jacobian carries `curr_fact · surf · Δramo_shape` weights that include
  displacement-current contributions. Dirichlet-probe `b_k` is a pure-DC
  construct and would miss these on transient / AC applications.
- **Sign is pinned by the forward code.** Whatever the forward EBIC
  driver reports as `I_EBIC(k)` is exactly what the adjoint computes —
  no independent sign convention to reconcile.
- **Robust at Schottky contacts.** With TE/TFE boundary conditions,
  "probe +1 Dirichlet" is ambiguous — the contact row is not a
  Dirichlet row. Ramo-Shockley is agnostic to the BC type at the
  collecting contact.

The equivalent volume-integral evaluation form
`δI_k ≈ −e Σ_i tr_vol_i [ν*_i + κ*_i] G_E_i` is unchanged; only the
assembly of `c_k` differs.
