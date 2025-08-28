# Schottky Contact Implementation: Technical Notes

## Fundamental Question: Is Direct FV Integration Equivalent to Robin BC?

**Answer: Yes, with critical requirements for correctness.**

## 1. Mathematical Equivalence

The direct finite volume integration approach **is** the Robin boundary condition:

- Robin BC mathematical form: `α n + β ∂n/∂n̂ = γ`
- FV boundary flux formulation: `j·n̂ = v(n - n₀)`
- FV integration: `∫_∂Ω j·dA = Σ A_face j·n̂`

Substituting the flux law directly into the boundary face integral **is the canonical FV implementation of Robin BC**.

## 2. Physical Correctness Requirements

### Necessary Conditions:
1. **Flux Definition**: j must be the total drift-diffusion particle flux, not diffusion alone
2. **Equilibrium Density**: n₀ must incorporate complete Schottky physics:
   ```
   n_inj(ψ) = N_c exp((ψ - V_c - Φ_B,eff)/V_T)
   Φ_B,eff = Φ_B - ΔΦ(E_n)  [with image force lowering]
   ```
3. **Discretization Consistency**: Boundary treatment must match interior scheme (Scharfetter-Gummel)

### Common Errors:
- Using diffusion-only flux: `j = -D ∇n` (incorrect under bias)
- Inconsistent unit systems (particle flux vs current density)
- Cell-center values instead of face interpolation
- Missing field dependence in barrier height

## 3. Critical Analysis of Standard Formulation

The standard formulation has inconsistencies that prevent absolute correctness:

### Unit System Inconsistency
**Problem**: Mixed notation between:
- Particle flux: `j_n` [m⁻²s⁻¹]
- Current density: `J_n = q·j_n` [A/m²]

**Solution**: Choose one consistently:
```
Particle flux form: ∂n/∂t + ∇·j_n = G - R
Current form: ∂n/∂t + (1/q)∇·J_n = G - R
```

### Scharfetter-Gummel Compatibility
**Problem**: Simple diagonal addition `A_face·s_n` assumes zero field approximation.

**Correct Implementation**:

For Slotboom variables `u_n = n exp(-ψ/V_T)`:
```
Boundary flux: j_n·n̂ = s_n[exp(ψ_f/V_T)·u_n - n_inj(ψ_f)]
Matrix diagonal: A_f·s_n·exp(ψ_f/V_T)
RHS contribution: -A_f·s_n·n_inj(ψ_f)
```

For standard variables with SG:
- Requires Bernoulli function weighting
- Diagonal contains both emission velocity and diffusion terms

## 4. Comprehensive Schottky Model Components

### Basic Thermionic Emission
```
s_n = (A*T²)/(q·N_c)  [Richardson formulation]
n_inj = N_c·exp(-(Φ_B - qV)/kT)
```

### Essential Physics Extensions

#### Image Force Barrier Lowering (IFBL)
```
ΔΦ = sqrt(q|E_n|/(4πε_s))  [SI units]
Φ_B,eff = Φ_B - ΔΦ
```
Update per Newton iteration using interface normal field.

#### Thermionic-Field Emission (TFE)
For N_D > 10¹⁸ cm⁻³:
```
E₀₀ = (q·ℏ/2)·sqrt(N_D/(m*·ε_s))
Implement Padovani-Stratton or WKB transmission
```

#### Interface States
SRH-like recombination:
```
R_interface = σ_n·v_th·N_st·(n·p - n_i²)/(n + p + 2n_i)
```

## 5. Implementation Strategies

### Option A: Direct FV with Explicit BC
```
For contact vertex k:
Ω_k ∂n_k/∂t + Σ(interior faces) J·dA + Σ(boundary faces) s_n(n_k - n₀)·A = Ω_k(G - R)
```

### Option B: Ghost Cell Elimination
1. Define ghost cell value n_ghost
2. Apply SG flux between cell and ghost
3. Enforce: `J_SG = q·s_n(n_boundary - n_inj)`
4. Eliminate n_ghost algebraically

Both yield identical matrix contributions when properly implemented.

## 6. Validation Criteria

### Numerical Requirements:
- Positive density preservation
- Monotonicity under high fields
- Newton convergence for Φ_B up to 1.0 eV

### Physical Benchmarks:
- Richardson equation fit in forward bias
- Saturation current: `J_s = A*T²·exp(-Φ_B/kT)`
- Ideality factor: n ≈ 1.0-1.1 (pure TE)
- Built-in potential extraction

### Limiting Cases:
- `s → ∞`: Recovers Dirichlet BC
- `Φ_B → 0`: Approaches ohmic contact
- Zero bias: Equilibrium with no net current

## 7. Conclusion

**The direct FV integration IS the correct Robin BC implementation.**

However, absolute correctness requires:
1. **Consistent flux formulation** (drift + diffusion)
2. **Proper discretization** (SG-compatible)
3. **Face value evaluation** (not cell centers)
4. **Complete barrier physics** (IFBL, TFE for advanced models)
5. **Consistent units** throughout

The approach is standard in production semiconductor simulators when these requirements are satisfied.