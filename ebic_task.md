## Task: Implement STEM-EBIC Generation Profile

Add electron-beam generation source for STEM-EBIC simulation. 200 keV beam through thin Si lamella (100-500 nm). Beam forms cone due to semi-convergence angle α = 10 mrad. At surface: σ₀ = 0.1 nm. At depth z: σ(z) = sqrt(σ₀² + (αz)²).

**Formula:** G(x,z) = G_tot × exp(-(x-x_B)²/(2σ(z)²)) / (sqrt(2π) × σ(z) × t), where G_tot = (I_beam/q) × (dE/dz × t) / E_ehp

**Parameters (Si):** I_beam=82.5pA, α=0.01rad, σ₀=0.1nm, dE/dz=580eV/μm, E_ehp=3.6eV

**Requirements:** Create module calculating G_array for given beam position. Integrate with mesh/solver. Verify conservation: sum(G×V) ≈ G_tot. Use CGS units.

**Verification:** t=500nm gives ~80 pairs/electron, G_tot~3e9/s, peak EBIC~0.5nA

Look at existing code first. Find where generation enters continuity equations.

When complete: <promise>EBIC_GEN_DONE</promise>
If stuck after 10 iterations: document blockers, output <promise>EBIC_GEN_BLOCKED</promise>
