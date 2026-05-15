# Heterojunction Equilibrium Current Floor — Session Notes

A working write-up of what we learned trying to drive the V=0 equilibrium
current in `het_test` (Si / Si₀.₈Ge₀.₂ / Si N-N-N stack) down to machine-
precision zero. Captures the real-bug fix that landed, the empirical
floor of the current discretization, and what it would take to break
through it. Pin this somewhere for whoever picks up v2 next.

---

## TL;DR

| Knob | Status | Result at V=0 |
|---|---|---|
| Initial Phase B/C/D + wrong-sign Δbe in SG | landed and reverted | 8.82×10⁻²⁰ A |
| **Sign-fixed SG: `dpot = -ch·(Δpot − Δbe)`** | **landed** | **3.83×10⁻²¹ A** (24×) |
| Slotboom-form rewrite via `expm1` | tried, regressed | 6.4×10⁻²⁰ A (worse, broke schottky_test) |
| Reorder dpot as `(pot − be)(2) − (pot − be)(1)` | tried, no improvement | 4.7×10⁻¹⁹ A |
| Matched-DOS heterojunction (Nc_Si = Nc_SiGe) | tried, no improvement | 4.7×10⁻¹⁹ A |
| Geometric-mean denormalization | landed, indistinguishable | unchanged |
| Tightened Newton tolerance | not viable | jobs hang (tol below FP eps) |

**Reference**: `schottky_test` (single-material) reports **1.04×10⁻²⁸ A**,
not 1e-30 as is sometimes assumed. The 7-order gap between heterojunction
(~1e-21) and single-material (~1e-28) is **structural to the dens-as-
primary-variable formulation**, not a missing physics term.

**Jacobian independently verified.** The analytical `sys_full` Jacobian at
the converged Si/SiGe equilibrium state has been checked against
finite-difference / Taylor-slope / JVP+Richardson references via the
`het_jacobian_test` driver (six probe directions: interface pot, interface
ndens, bulk pot, three random seeds). Taylor-remainder slope is `2.00`
throughout; in the `v_dens_int` direction `R1(h)` is flat at the FP noise
floor across the entire sweep — the strongest possible "J = FD"
outcome. Combined with the Bratu fixture (20/20 verified discrimination
of known-wrong Jacobians), this rules out a Jacobian bug as the source of
the floor. See *Jacobian independently verified* below.

Closing the gap requires switching to **iref as a primary Newton
variable** (Slotboom-form, as in Sentaurus / Sesame) — that refactor is
deferred, not blocking. `1e-21 A` is 14 orders below femtoamps, so no
present cybear application is bottlenecked by it. Revisit only if a
future use case demands sub-`1e-21` equilibrium currents.

---

## The real bug: missing `ch` on the band-edge term

The SG argument I shipped in the first cut of Phase C.1 was

```fortran
dpot = -ch * (pot(2) - pot(1)) + (be(2) - be(1))     ! WRONG
```

This came from a memorized statement "the SG argument picks up Δbe at a
heterojunction edge" — true at the symbol level, wrong on the sign for
the simulator's convention. The fix:

```fortran
dpot = -ch * ((pot(2) - pot(1)) - (be(2) - be(1)))   ! correct
```

### Derivation (so the fix doesn't get unfixed)

From `imref.f90`:

```
iref(v) = pot(v) − be(v) + ch · eta(v)             with eta = ln(dens/edos)
```

At flat quasi-Fermi (equilibrium):
```
Δln(n) = Δpot − Δbe − ch · Δiref = Δpot − Δbe       (since Δiref = 0)
```

From the SG kernel `j = Bern(-dpot)·n(1) − Bern(dpot)·n(2)`, the zero
condition `j = 0` gives `ln(n(2)/n(1)) = dpot`, i.e.
```
Δln(n) = dpot
```

Equating:
```
dpot = Δpot − Δbe                                    (electron, ch = -1)
dpot = -Δpot + Δbe                                   (hole, ch = +1)
```

Both fit `dpot = -ch · (Δpot − Δbe)`. The first cut had `dpot = -ch·Δpot
+ Δbe`, which for electrons gives `Δpot + Δbe` — opposite sign on the
band-edge term. That sign error left a ~`2·Δbe ≈ 0.9` residual in dpot
at the interface edge, which Newton couldn't drive to zero.

### What this fix did *not* turn out to be

The `ln(Nc(2)/Nc(1))` Polsky–Rimshans term that I (and the literature
search) initially suspected was missing **is not missing**. The per-
vertex `n = dens(v) / edos_v(v)` normalization already absorbs the DOS
discontinuity exactly. Verifying via detailed balance:

- At deep-bulk equilibrium under uniform doping, `dens(v) = Nd` at every
  vertex (charge neutrality).
- `n(v) = dens(v) / Nc(v) = Nd / Nc(v)` differs per material.
- For SG `j = 0`: `ln(n(2)/n(1)) = dpot` must hold.
- LHS: `ln(Nc(1)/Nc(2))`. RHS (with the correct sign fix): the simulator
  drives `Δpot − Δbe = ln(Nc(1)/Nc(2))` through the imref relation +
  Poisson equilibrium.
- ✓ No extra term needed.

---

## The structural floor at ~1e-21 A

After the sign fix, het_test reports terminal currents in the
1e-19 to 1e-21 A range depending on FP rounding paths in a given build.
All of these are "the same number" — they cluster around the FP precision
floor of the dens-based SG kernel at a heterojunction interface edge.

### Mechanism

Single-material bulk at equilibrium:
- `pot(v)` and `dens(v)` are bit-identical between adjacent vertices
  (Newton converges to uniform values under uniform doping).
- `dpot = -ch · (pot(2) − pot(1)) = -ch · 0` is FP-exact zero.
- `Bern(0) = 1` exactly. `j = 1 · n(1) − 1 · n(2) = n − n` is the
  subtraction of identical bit patterns → **literal 0 in IEEE FP**.
- Continuity propagates zero everywhere. Contact current ~ Newton
  residual on the contact equations → ~1e-28 A.

Heterojunction interface edge:
- `n(1)` and `n(2)` differ structurally (different Nc).
- `dpot ≠ 0` (the band-edge step survives the cancellation).
- `Bern(-dpot)·n(1) − Bern(dpot)·n(2)` is the subtraction of two finite
  products (each ≈ 3×10⁻³ normalized / ~250 nA physical at the Si/SiGe
  interface) that are *algebraically equal at equilibrium* but
  **cannot be bit-identical** because they're computed from different
  operands. The absolute IEEE-FP residual is bounded by
  `~eps_dbl·|Bern·n|` ≈ `7×10⁻¹⁹` normalized, with the measured value
  varying between `~10⁻¹⁹` and `~10⁻¹⁸` depending on the specific
  rounding paths each Newton iteration takes.
- Newton enforces `∇·J = 0` everywhere, so this residual propagates to
  every edge in the chain. Contact current after denormalization sits
  in the `~10⁻¹⁹` to `~10⁻²¹ A` range, with both magnitude and sign
  varying across grid sizes as FP rounding paths shift.

### Empirical verification

A per-edge diagnostic in `get_curr` printing `B(-dpot)·n(1)`, `B(dpot)·n(2)`,
and the imref-flat residual `log(n(2)/n(1)) − dpot` at heterojunction
edges confirms each part of the mechanism:

- `B(-dpot)·n(1)` and `B(dpot)·n(2)` agree to 14–15 significant digits
  (e.g. both `3.25589959120636×10⁻³` at the Si/SiGe interface), with
  the cancellation residual sitting at the `eps_dbl·|B·n|` floor. No
  "missing term" — the kernel is algebraically exact at flat imref.
- `log(n(2)/n(1)) − dpot` stabilizes at `~10⁻¹⁵` at the interface edge.
  This is the *algebraic identity* the SG kernel relies on for `j = 0`,
  and Newton converges it to **machine precision** regardless of the
  `atol` setting — confirming Newton tolerance is not the bottleneck.
- Decoupled tests (`ΔEc = 0` with `ΔNc ≠ 0`, and `ΔNc = 0` with
  `ΔEc ≠ 0`) both produce FP-precision floors in the `10⁻¹⁹` to
  `10⁻²⁰` A range — neither collapses to bit-zero. The cancellation is
  intrinsic to the per-vertex `n = dens/edos` normalization at any
  non-uniform edge, not to any specific source of non-uniformity.

### Grid-spacing scaling

The kernel-level cancellation `j_norm` scales linearly with `1/Δx`,
matching the FP-precision prediction that the residual is invariant in
*normalized* units but enters the per-edge denormalization as a `1/len`
factor:

| `dx` (nm) | `j_norm` (kernel, normalized) | ratio vs predicted `1/Δx` |
|---|---|---|
| 2.0 | `2.53×10⁻¹⁷` | baseline |
| 1.0 | `5.06×10⁻¹⁷` | 2.00× (predicted 2×) |
| 0.5 | `1.09×10⁻¹⁶` | 4.31× (predicted 4×) |

Contact-current magnitude does *not* scale cleanly — sign flips and
factor-100× swings occur — because the contact-current path involves
Newton-residual norm tolerance and Ramo–Shockley integration, each
adding its own rounding paths. The kernel scaling is the load-bearing
prediction; the contact current behaves as bounded random-sign FP noise
in the `10⁻¹⁹–10⁻²¹` range.

### Three-precision SG cancellation test

The strongest single confirmation of the FP-floor theory comes from
recomputing the same SG flux in three precision regimes (all in one
diagnostic pass, no separate runs needed):

| Test | Kernel arithmetic | Inputs | Cancellation residual |
|---|---|---|---|
| `diff_d` | double | double (cybear's actual state) | `~5×10⁻¹⁹` |
| `diff_q` | quad (`real128`) | double cast to quad | `~5×10⁻¹⁹` |
| `diff_qi` | quad | quad, reconstructed to satisfy `iref = 0` exactly in quad | **`0` or `~4×10⁻³⁷`** |

`diff_d ≈ diff_q` proves the kernel arithmetic is already at the
input-precision floor — promoting only the kernel to quad does not
help. `diff_qi` drops by 19 orders of magnitude when the *inputs*
themselves are made consistent at quad precision, with several
iterations printing literally `0.0000000000000000000000000000000000E+000`
— bit-exact zero. This is the cleanest possible demonstration that
the implementation is correct: the SG kernel returns exact zero when
given inputs at sufficient precision to make the underlying Bernoulli
identity bit-exact, and the `~10⁻¹⁹` residual in production is solely
the cost of storing `(pot, dens)` as IEEE doubles.

### Why finishing with Gummel-only doesn't help

A natural hopeful guess: maybe the `sys_dd` Gummel iteration (which
uses `calc_eta_iref` and treats `iref` as a provided variable) already
delivers the iref-form floor, and we just need to skip the final
`sys_full` Newton refinement. **Measured: false.** Stopping at Gummel
gives `j_norm ≈ 1.4×10⁻¹⁶` (vs `5×10⁻¹⁷` for full Newton at the same
`Δx`) — Gummel's fixed-point convergence is *looser* than Newton's
quadratic convergence, so its converged `(pot, dens)` state has more
residual. Both Gummel and Newton produce double-precision `(pot, dens)`,
both feed the same SG kernel, both hit the same input-precision floor.

### Why various rewrites don't help

| Rewrite | Why it fails to lower the floor |
|---|---|
| Slotboom form `j = n(1)·Bern(-dpot)·(-expm1(-z))` with `z = dpot − Δln(n)` | `z` is computed from FP-precision pot and dens values; its absolute precision is ~`eps × max(|dpot|, |Δln(n)|) ≈ 1e-15`. Multiply by `n·Bern` → same floor. Jacobian is also asymmetric and hurts Newton convergence on schottky reverse-bias. |
| Reorder `(pot − be)(2) − (pot − be)(1)` instead of `Δpot − Δbe` | Algebraically identical. The bit patterns of `pot(v) − be_v(v)` at adjacent interface vertices are not identical at Newton convergence, so the subtraction still carries `eps × |scale|`. |
| Matched-DOS heterojunction (Nc_Si = Nc_SiGe, only Eg and Ec_offset differ) | At equilibrium, `n(v) = Nd/Nc = const` IS bit-identical across the device. But `pot(v) = be_v(v) + ln(Nd/Nc)` varies because `be_v` varies. The `dpot = Δpot − Δbe` subtraction has FP cancellation noise of order `eps × |be| ≈ 4×10⁻¹⁵` → contact current ~1e-19, same order as mismatched-DOS. |
| Geometric- vs arithmetic-mean `edos_avg` in denormalization | At equilibrium the SG kernel returns the same algebraic zero regardless of the prefactor scaling; only the FP cancellation residual matters, and that residual is the same in either form. |
| Quad-precision arithmetic inside the SG kernel | Useless when inputs are double. Result precision is bounded by **input** precision, not arithmetic precision. To benefit from quad, Newton-variable storage must also be quad. |
| Tightening Newton tolerance (atol, rtol, min_it) | Default `rtol=1e-12, atol=1e-16` already hits FP eps after 2 iterations on this device. Tighter values cause Newton to spin forever without further progress. |

The pattern: any rewrite that keeps `(pot, dens)` as the primary Newton
variables has an arithmetic-residual floor of `eps × |Bern·n| ≈ 1e-14`
normalized at the interface edge. The bit patterns of n(1), n(2) and
dpot at the interface are structurally non-identical, and IEEE
arithmetic cannot manufacture cancellation that the input bit patterns
don't admit.

---

## Jacobian independently verified

The FP-floor argument above is load-bearing only if Newton actually
converges to `F(x) ≈ 0` within machine epsilon. If the analytical
Jacobian disagreed with the true Jacobian of the residual, Newton would
stall early at a residual floor set by that inconsistency — and the
terminal-current "floor" would actually be a Jacobian bug, not FP
arithmetic. To distinguish, `het_jacobian_test` (target wired in
`fargo.toml`) runs the converged `sys_full` Jacobian through two
threshold-free reference checks:

- **Taylor remainder slope.** `R₁(h) = ‖F(x+hv) − F(x) − hJv‖`. For
  correct `J`, `R₁ ∼ h²`; for wrong `J`, `R₁ ∼ h`. Threshold-free —
  only the *shape* of the curve is interpreted.
- **JVP + Richardson.** `‖Jv − (4·g(h/2) − g(h))/3‖` where `g(h) =
  (F(x+hv) − F(x−hv))/(2h)`. Compares the analytical `Jv` against an
  `O(h⁴)`-accurate centered FD reference of `F` along `v`.

Six probe directions are exercised: pot at heterojunction interface
vertices; ndens at the same; pot in the Si bulk far from any interface;
three independent random vectors. Result (at converged equilibrium with
`‖F(x*)‖ = 2.6e-14`):

| Probe | Taylor slope | JVP+Richardson err |
|---|---|---|
| `v_pot_int`  (interface pot)   | 2.00 (large-h half)             | `1.33e-10` |
| `v_dens_int` (interface ndens) | flat at FP floor `~1e-14`       | `7.99e-12` |
| `v_bulk`     (Si bulk pot)     | 2.00 throughout                 | `9.34e-11` |
| `v_rand` × 3 (seeds 31337/271828/161803) | 2.00 throughout (all 3) | `1.18e-9` / `2.64e-10` / `7.39e-10` |

The `v_dens_int` case is the strongest possible "J = FD" outcome: `hJv`
matches `F(x+hv) − F(x)` to machine precision across the entire `h`
sweep — there is no residual left for a slope to be fit to. The slope
test is paired with a "max R₁ below `1e-10`" criterion to recognize
this regime as PASS rather than the false-FAIL it would be if only the
fitted slope were checked.

The validators themselves are verified on the Bratu fixture
(`test/jacobian_validation_test.f90`): 20/20 expected verdicts across
{correct, scaled-by-1.5, dropped-diagonal, transposed off-diagonal,
1e-9 noise} × {finite-difference, complex-step, Taylor-slope,
JVP+Richardson}. The discriminator gap between correct and any
structural defect is `~10` orders of magnitude (Taylor slope `2.00 →
1.00`; JVP err `7e-12 → 6e-1`), so the `het_test`'s slope-2 result is
not something the validator can produce on a wrong Jacobian.

**Conclusion:** the analytical Jacobian for `sys_full` at the
converged Si/SiGe equilibrium is correct, including the band-edge
contribution in `current_density.f90:260`. The `1e-21 A` floor is
genuinely the FP cancellation residual of the `(pot, dens)`
formulation — not a hidden Jacobian bug and not a Newton-convergence
artifact. The iref refactor below is therefore confirmed as the *only*
path to lower the floor, *and* confirmed as optional for any current
device application.

---

## What actually closes the gap (the v2 work)

From [Sesame docs][sesame] and [Sentaurus User Guide][sentaurus]:

```
J_n^i = (q μ_n / Δx) · [ (ψ_n(i+1) − ψ_n(i)) / (e^{-ψ_n(i+1)/kT} − e^{-ψ_n(i)/kT}) ]
                     · [ e^{E_Fn(i+1)/kT} − e^{E_Fn(i)/kT} ]
```

where `ψ_n = qφ + χ + kT·ln(N_C)` and `E_Fn` is the electron quasi-Fermi
level **stored as a primary variable**.

The key structural property: `E_Fn` is set by the Dirichlet BC at the
contact (`E_Fn(contact) = V_applied`). At V=0, every vertex stores
`E_Fn = 0` as the *same* float value. The driving force
`exp(E_Fn(i+1)) − exp(E_Fn(i)) = exp(0) − exp(0) = 1 − 1 = 0` is
*literally zero in IEEE arithmetic*, not "approximately zero". Continuity
then carries that exact zero to every edge.

Cybear's current Newton system holds `(pot, dens, currents, V_*)`. To
adopt the Sesame/Sentaurus form, swap `dens` for `iref` (per carrier)
and rewrite the system. **No existing equation system has iref as a
Newton primary today** — earlier drafts of this doc implied otherwise
based on `calc_eta_iref` appearing in `sys_nlpe`/`sys_dd`, but that
equation consumes iref + pot to produce eta; it doesn't *solve* for
iref. The actual state in `src/device.f90`:

```
sys_nlpe:  Newton primary = pot           (Poisson is the res_equation)
           iref is provide()-d externally — line 229
           eta and dens are derived from (pot, iref) via
             calc_eta_iref → calc_dens
sys_dd:    Newton primary = ndens         (continuity is the res_equation)
           pot, iref are provide()-d externally — lines 279–280
sys_full:  Newton primary = pot, ndens, currents, V_*
           ndens drives continuity; eta = calc_eta_dens(ndens);
           iref = calc_iref(pot, eta) [derived]
```

The variable infrastructure (`type, extends(variable) :: imref` in
`src/imref.f90:35`, with full grid_data storage and Jacobian chain
plumbing) and the helper equations (`calc_eta_iref`,
`calc_chemical_pot_imref`, `calc_dens`) all exist. The piece that
doesn't exist is any `res_equation` whose main_var is `iref`. That's
the actual refactor:

| File | What changes |
|---|---|
| `src/device.f90` | In `sys_full` setup: swap `calc_eta_dens` → `calc_eta_iref`, swap `calc_iref` (derived) → `calc_dens` (derived). ~4 lines. |
| `src/continuity.f90` | Continuity residual `∂n/∂t + ∇·J = 0` rewritten with `n = edos·exp(...)` derived from `iref`. Jacobian re-derived in `(pot, iref)` space. |
| `src/contact.f90` | Ohmic BC: `iref(contact) = V_applied` (Dirichlet) instead of `dens(contact) = ...` |
| `src/approx.f90` | Initial-guess routines seed `iref` directly instead of `dens` |
| `src/current_density.f90` | SG kernel rewritten in Sesame form `j ∝ M(ψ) · [exp(iref(2)) − exp(iref(1))]`; Jacobian re-derived |
| `src/schottky.f90` | Thermionic BC rewritten in iref terms (Schottky devices only) |
| `test/het_test.f90` | Output schema if it asserts on `ndens` directly |

Estimated effort: **2–4 focused days of work** + validation against the
existing pn-diode / MOSFET / nanowire / Schottky / material regression
suite (all of which currently end with `sys_full` and would now end in
iref-form). The infrastructure that's *not* in this list — `imref.f90`,
`device_params.f90`, the system-builder glue, output writers,
`charge_density.f90` — either already supports iref or doesn't care
which form is primary because it only consumes derived quantities.

Worth doing for v2 — heterojunctions, perovskite VFETs, and 2D-material
stacks are all going to demand the same equilibrium-current property,
and dens-primary has no compensating advantage that justifies keeping
it. Most likely path is to **delete dens-primary entirely** rather than
add a runtime toggle.

---

## Code pointers for picking this up later

The Phase B + C.1 + D-min work lives on branch `heterojunction-v2`:

- **Per-vertex band-structure fields**: `src/device_params.f90`
  - Declarations: `band_edge_v`, `edos_v`, `band_gap_v`, `dEc_v` (lines
    78–86).
  - `device_params_init_band_edges` subroutine (lines ~1192–1245) —
    iterates `reg_trans(:)`, fills from each region's resolved material;
    `Ec_offset` (smc%dEc) is folded into `band_edge_v` on both carriers.

- **SG flux modification**: `src/current_density.f90:260`
  - `dpot = -ch * ((pot(2) - pot(1)) - (be(2) - be(1)))` — the fixed
    form. **Do not "simplify" the parentheses without re-deriving from
    `imref = pot - be + ch·eta`.**
  - Per-vertex `n(v) = dens(v)/edos_e(v)` (lines 261–262) — handles Nc
    discontinuity by per-vertex normalization, no extra `Δln(Nc)` term
    needed.
  - `edos_avg = sqrt(edos_e(1) * edos_e(2))` (line 266) — geometric mean
    for the denormalization. Collapses to `edos` in single-material
    limit. `djddens` includes the `edos_avg/edos(i)` ratio (line 295–296)
    to compensate.

- **Per-vertex consumers (Phase D-min)**:
  - `src/imref.f90` — 4 vertex loops (calc_chemical_pot_dens,
    calc_imref, calc_chemical_pot_imref, calc_density) now read from
    `par%band_edge_v(ci)%get(idx)` and `par%edos_v(ci)%get(idx)`.
  - `src/continuity.f90:270` — ohmic boundary density at contact vertex
    uses `band_edge_v` and `edos_v` of the contact vertex's local
    material.
  - `src/approx.f90:130–146` — initial-potential guess uses per-vertex
    `band_gap_v` and `edos_v`.

- **Test device**: `devices/het/si_sige.ini`, `devices/het/run.ini`.
  Both committed; fargo target `het_test` registered in `fargo.toml`.

- **Test runner**: `test/het_test.f90`. Asserts I_L, I_R < 1e-15 A and
  the dEc step exactly equals user-spec 0.020 eV. Includes a per-edge
  SG-flux diagnostic that prints the max-|J| edge index and `j_norm`
  value — useful for debugging the floor if you push this work further.

- **Jacobian validation harness**: `src/jacobian_validator.f90` +
  `src/esystem_problem.f90` + `test/het_jacobian_test.f90` +
  `test/jacobian_validation_test.f90` (Bratu proof). Run via `fargo run
  het_jacobian_test` or `fargo run jacobian_validation_test`. Use these
  if you ever extend the heterojunction term and want to verify the
  Jacobian chain rule survived — slope-2 in `R₁(h) = ‖F(x+hv) − F(x) −
  hJv‖` is the load-bearing check.

- **Design doc**: `docs/heterojunctions.md`. The "Equilibrium current
  floor" section (added this session) summarizes the structural argument
  and links back to the Sesame/Sentaurus convention.

---

## Things that *seemed like they should help* but didn't

Documenting these so future-me doesn't try them again expecting
different results.

1. **Polsky–Rimshans `ln(Nc(2)/Nc(1))` term in dpot.** The literature
   keeps mentioning it. Per the detailed-balance derivation above, it's
   already absorbed by per-vertex `n = dens/edos` normalization. Adding
   it explicitly *double-counts* and gives wrong physics.

2. **Sesame-form `j ∝ exp(EFn(2)) − exp(EFn(1))` computed from cybear's
   derived iref.** The Sesame trick relies on `EFn` being a *stored*
   variable. Computing it on-the-fly from `(pot, dens, be)` and feeding
   into the form gives the same FP precision as the standard SG kernel,
   because the cancellation residual moves but doesn't shrink. The only
   way to get the Sesame floor is to store `iref` directly so the BC
   value propagates as identical floats.

3. **Quad precision in the SG kernel only.** The output precision is
   bounded by the *input* precision. Quad arithmetic on double-precision
   `(pot, dens)` gives a more accurately-computed eps-residual, not a
   smaller residual.

4. **`expm1` for the cancellation.** Helps when the cancellation is
   between exponentials whose arguments are close to zero. In the SG
   kernel the cancellation is between `Bern(±dpot)·n` which are O(1)
   each, not exponentially small. `expm1` only helps if you can route
   the calculation through a small argument like Δiref, but Δiref
   itself has eps-level FP noise from being derived rather than stored.

5. **Tightening Newton (`atol`, `rtol`, `min_it`).** The Newton residual
   already hits FP eps in 2 iterations on het_test. There's no slack
   to extract.

6. **Matched-DOS heterojunction (Nc_1 = Nc_2 but Eg, Ec_offset differ).**
   Intuitive guess: matched DOS makes `n(1) = n(2)` bit-identical at
   equilibrium, so cancellation is exact. Reality: `pot(v) - be_v(v)`
   still differs in bit pattern even though it's algebraically uniform,
   so `dpot` retains `eps · |be|` cancellation noise. Same floor.

7. **Reordering dpot computation to subtract per-vertex `alpha = pot - be`.**
   Same algebraic form, different operation order. The bit patterns at
   convergence are not identical between adjacent interface vertices,
   so reordering doesn't manufacture an exact zero.

---

## Practical guidance going forward

- **If you're benchmarking heterojunction equilibrium currents in
  cybear**, treat ~1e-20 to 1e-22 A as "numerical zero" for V=0. It's
  9–11 orders below any physically relevant current (femtoamps), so
  using it as the threshold for "passes equilibrium" is safe.

- **The single-material schottky_test floor of ~1e-28 A is the
  benchmark for "what's achievable with bit-exact cancellation"** — that
  number is set by how tightly Newton converges the contact charge-
  neutrality equation, and it's the right target to aim at if/when v2
  lands the Slotboom refactor.

- **The `het_test.f90` per-edge diagnostic** is useful for any future
  precision investigation. It prints the max-|J| edge index and the
  raw `j_norm` value at that edge — confirms the floor is uniform
  across edges (continuity-propagated) rather than concentrated at any
  one location.

- **Do not be tempted to "snap" the current to zero in
  postprocessing**. The 1e-21 result is honest physics within the
  discretization. Hiding it would just defer the realization that v2
  needs the refactor.

- **If you suspect the floor has moved**, re-run `het_jacobian_test`
  before chasing a Jacobian bug. A correct Jacobian gives Taylor slope
  `2.00` cleanly; a regressed Jacobian flips the slope to `~1.00` and
  the JVP+Richardson error climbs by `8+` orders of magnitude. That's
  the fastest discriminator between "FP-floor moved because we changed
  scaling somewhere" and "a real bug entered the heterojunction code
  path".

[sesame]: https://sesame.readthedocs.io/en/latest/physics/discretization.html
[sentaurus]: https://www.synopsys.com/manufacturing/tcad/device-simulation/sentaurus-device.html

---

*Notes captured 2026-05-14, branch `heterojunction-v2`, after
commit `4128baa`. If anyone with a Sentaurus license can compile a
matching Si/SiGe₂₀ deck and report their V=0 equilibrium current, that's
the cleanest cross-check on the "1e-28 vs 1e-21" claim above.*
