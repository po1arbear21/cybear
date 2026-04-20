# Linearity Test Spec — EBIC Adjoint Prerequisite

## Objective

Verify that the forward map $F: G \mapsto I_\text{EBIC}$ can be treated as linear around the dark operating point $\mathbf{u}_0$, so that one adjoint solve yields a generation-independent collection efficiency $\eta(\mathbf{r})$ satisfying $I_\text{EBIC} = \int G(\mathbf{r})\, \eta(\mathbf{r})\, d\mathbf{r}$ to within experimental accuracy.

## Setup

**Device.** The STEM-EBIC Si p–n lamella with the actual W top / Pt bottom contact stack used in the SISPAD draft — not a simplified test structure. Linearity can hold on one geometry and fail on another, and only the production geometry is decision-relevant.

**Bias.** Short-circuit (collection condition). If multiple bias points appear in the paper, repeat the scaling test at each.

**Operating point.** Explicitly precompute and store the dark state $\mathbf{u}_0 = (\psi_0, n_0, p_0)$ from a converged Newton solve with $G \equiv 0$. All linearizations and perturbation metrics reference this state.

**Beam model.** Use the existing cybear STEM generation rate $G(\mathbf{r}; \mathbf{r}_\text{beam}, I_\text{beam})$. Confirm that $I_\text{beam}$ enters as a clean multiplicative prefactor on $G$ — this is assumed by Test 1 and should be checked in the code once.

## Tests

### Test 1 — Homogeneity (scaling)

Fix one beam position $\mathbf{r}_A$ near the junction (in the depletion region). Solve the forward problem at $I_\text{beam} \in \{10^{-2}, 10^{-1}, 1, 10, 10^2\} \times I_\text{exp}$, where $I_\text{exp}$ is your experimental beam current (~10–100 pA). Plot $I_\text{EBIC}$ vs $I_\text{beam}$ on log–log.

*Pass:* slope = 1.000 ± 0.001 over at least the two decades bracketing $I_\text{exp}$, zero intercept. Deviation at the high-$I_\text{beam}$ end locates where high injection breaks linearity and is useful information to report regardless.

### Test 2 — Superposition

Pick two beam positions:
- $\mathbf{r}_A$: quasi-neutral n-side, roughly one diffusion length from the junction
- $\mathbf{r}_B$: quasi-neutral p-side, similarly placed

Solve three configurations at $I_\text{beam} = I_\text{exp}$: $G_A$ alone, $G_B$ alone, $G_A + G_B$ simultaneously. Compute the residual
$$r = \frac{|I(G_A + G_B) - I(G_A) - I(G_B)|}{|I(G_A) + I(G_B)|}.$$

*Pass:* $r < 10^{-2}$ at $I_\text{exp}$. Repeat at $10\, I_\text{exp}$ and $100\, I_\text{exp}$ to find where superposition breaks — this is the cleanest direct test of linearity and the one I'd lead with when presenting to Prof. Jungemann.

### Test 3 — Small-signal condition

For each run in Test 1, compute over the full domain:
- $\max |\Delta n| / n_0$, $\max |\Delta p| / p_0$
- $\max |\Delta \psi| / V_T$
- quasi-Fermi splitting $\max (E_{Fn} - E_{Fp}) / k_B T$ in regions where it was ≤ 0 in the dark state

*Pass:* all ratios $< 0.1$ at $I_\text{exp}$. Map where the largest perturbations sit — if they concentrate at the Schottky contacts, that's directly relevant to the sign-reversal hypothesis and worth flagging separately.

### Test 4 — Jacobian-consistency check (recommended, cheap)

Compute $\delta \mathbf{u} = \mathbf{u}(G) - \mathbf{u}_0$ from a forward solve, then predict
$$\delta I_\text{EBIC}^\text{lin} = \left(\frac{\partial I_\text{EBIC}}{\partial \mathbf{u}}\right)^{\!\top}_{\!\mathbf{u}_0} \delta \mathbf{u}$$
and compare with the true $\delta I_\text{EBIC} = I(G) - I(0)$.

*Pass:* agreement to $< 10^{-2}$ relative error at $I_\text{exp}$. This is a Taylor-remainder test at the observable level and reuses the same Jacobian infrastructure you already built for the solver-level TRT. Failure here with Tests 1–3 passing would indicate the response functional is being evaluated inconsistently with the linearization point.

## Outputs

1. Scaling plot (log–log, with fitted slope and R²) for Test 1.
2. Table of $(I_A, I_B, I_{A+B}, r)$ at each tested $I_\text{beam}$ for Test 2.
3. Domain maps of $\Delta n / n_0$, $\Delta \psi / V_T$ at $I_\text{exp}$ and at the current where linearity first breaks.
4. One-paragraph conclusion: is the experimental operating point inside the linear regime, and by what margin?

## Pitfalls to watch

- **Linearization point.** If the adjoint Jacobian is assembled from the *final* Newton iterate of a beam-on solve rather than from $\mathbf{u}_0$, the whole analysis is self-inconsistent. Worth an explicit check in code.
- **Solver tolerance floor.** Newton residual tolerance must be tightened below the linearity deviation you're trying to measure, otherwise Test 2's residual is dominated by solver noise. Record residuals alongside results.
- **Contact BCs.** If IFBL $\beta$ is evaluated using beam-perturbed fields/densities at the Schottky contacts, you've got an additional nonlinearity localized at exactly the place that matters for your sign hypothesis. Check whether $\beta$ is held at dark-state values in the linearization.
- **Generation prefactor.** Confirm once that doubling `I_beam` in the input exactly doubles the assembled RHS — it should be, but a unit-conversion quirk could fake nonlinearity.
