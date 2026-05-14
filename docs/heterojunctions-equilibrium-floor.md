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

Closing the gap requires switching to **iref as a primary Newton
variable** (Slotboom-form, as in Sentaurus / Sesame). That's a multi-day
refactor across ~6 files and all the equation-system glue.

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
  O(10)-magnitude products that are *algebraically equal at equilibrium*
  but **cannot be bit-identical** because they're computed from different
  operands. The IEEE-FP residual is `≈ eps × |Bern·n| ≈ 1e-14` in
  normalized units.
- Newton enforces `∇·J = 0` everywhere, so this 1e-14 propagates to
  every edge in the chain. Contact current ≈ 1e-14 × denormalization
  factor ≈ 1e-20 to 1e-21 A.

I verified this empirically via the per-edge diagnostic in `het_test.f90`:
the max-|J| edge sits *inside* the Si bulk (around index 65 or 265 on
the 301-vertex grid), with `j_norm ≈ 5×10⁻¹⁷` — exactly the residual the
interface edge contributes, propagated via continuity.

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
and rewrite the system. Concretely the refactor touches at least:

| File | What changes |
|---|---|
| `src/device.f90` | Newton variable wiring; replace `ndens`/`pdens` with `n_iref`/`p_iref` in `sys_full` |
| `src/device_params.f90` | Output dataset includes iref; perhaps drop derived-dens caching |
| `src/dd.f90` | BC: `iref(contact) = V_applied` (currently `dens(contact)` is the Dirichlet BC) |
| `src/imref.f90` | The relation flips direction — dens becomes derived from `(pot, iref)`; current calc-imref equation removed or repurposed |
| `src/current_density.f90` | SG kernel rewritten in Sesame form using `iref` inputs; Jacobian re-derived |
| `src/continuity.f90` | Equation residual expressed in `iref` not `dens` |
| `src/contact.f90` | Ohmic contact BC rewritten in terms of iref |
| `src/schottky.f90` | Schottky BC rewritten (thermionic emission in terms of iref) |
| `src/approx.f90` | Initial guess in (pot, iref) space |
| `src/charge_density.f90` | Computes ρ from `dens = edos·exp(eta(pot, iref))` |
| `test/schottky_test.f90`, `test/het_test.f90` | Output variable changes |
| `lib/fortran-basic/.../equation/*` | Possibly the system-builder glue if it assumes specific variable structure |

Estimated effort: 3–5 focused days of work + validation against the
existing pn-diode / MOSFET / nanowire / Schottky / material regression
suite (all of which currently depend on `(pot, dens)` being primary).
Worth doing for v2 — heterojunctions, perovskite VFETs, and 2D-material
stacks are all going to demand the same equilibrium-current property.

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

[sesame]: https://sesame.readthedocs.io/en/latest/physics/discretization.html
[sentaurus]: https://www.synopsys.com/manufacturing/tcad/device-simulation/sentaurus-device.html

---

*Notes captured 2026-05-14, branch `heterojunction-v2`, after
commit `4128baa`. If anyone with a Sentaurus license can compile a
matching Si/SiGe₂₀ deck and report their V=0 equilibrium current, that's
the cleanest cross-check on the "1e-28 vs 1e-21" claim above.*
