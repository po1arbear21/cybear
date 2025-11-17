# Schottky Contact V3: Uniform Tsu-Esaki Model

**Author**: Chenyang Yu
**Date**: 2025-11-17
**Version**: 3.0 (Clean Slate Implementation)

---

## Core Physics: Uniform Tsu-Esaki Model

### Single Integral Formulation

The total current density at a Schottky contact is computed using a single integral:

```
J_contact = -(4πq·m*·m₀·kT/h³) ∫₀^∞ T(E) · N(E) dE
```

where:
- **T(E)**: WKB transmission probability
- **N(E)**: Fermi-Dirac occupancy difference
- **Negative sign**: Electron current convention (electrons leaving → negative current)

### WKB Transmission Function

For a triangular barrier under electric field F:

```
T(E) = {
  exp[-(4/3)√(2m*) · (φ_b - E)^(3/2) / |F|]  for E < φ_b
  1.0                                          for E ≥ φ_b
}
```

### Occupancy Difference

Using logarithmic form for numerical stability:

```
N(E) = ln[1 + exp(-E)] - ln[1 + exp(-E - φ_k)]
```

where φ_k is the local electrostatic potential at the contact.

### Integral Decomposition (Optional)

For analysis and optimization, the integral can be split at the barrier height:

```
J_total = J_tunnel + J_thermionic

J_tunnel = -(m*/(2π²)) ∫₀^φ_b T(E)·N(E) dE     [under-barrier]
J_thermionic = -(m*/(2π²)) ∫_{φ_b}^∞ N(E) dE    [over-barrier, T=1]
```

**Note**: This split is purely for computational efficiency, not a physical split model.

---

## Implementation Specification

### Module Structure

```fortran
module schottky_tsuesaki_m
  ! Single module for all Schottky transport

  use device_params_m
  use potential_m
  use quad_m       ! Quadrature integration
  use math_m       ! PI constant

  implicit none

  type :: schottky_current
    real :: J         ! Total current density
    real :: J_tunnel  ! Tunneling component (0 to φ_b)
    real :: J_therm   ! Thermionic component (φ_b to ∞)
  end type

contains

  function calculate_schottky_current(phi_b, F, m_star, phi_k) result(J)
    ! Main entry point - returns total current
  end function

  function wkb_transmission(E, phi_b, F, m_star) result(T)
    ! Pure WKB transmission probability
  end function

  function occupancy_difference(E, phi_k) result(N)
    ! Fermi-Dirac occupancy difference
  end function

  subroutine integrand(E, params, f, dfdE, dfdp)
    ! Integrand for quad module interface
  end subroutine

end module
```

### Boundary Condition Interface

```fortran
! In continuity.f90
subroutine apply_schottky_bc(...)
  use schottky_tsuesaki_m

  ! Get current density
  J = calculate_schottky_current(phi_b, F, m_star, phi_k)

  ! Apply as Neumann BC (current density boundary condition)
  ! No Robin BC, no split injection/recombination terms
  residual = ... + A_contact * J

end subroutine
```

### Parameter Structure

```fortran
type :: schottky_params
  real :: phi_b      ! Barrier height (eV)
  real :: m_tunnel   ! Tunneling effective mass (m*/m0)
  logical :: ifbl    ! Image force barrier lowering flag

  ! Integration control
  real :: rtol = 1e-6  ! Relative tolerance
  integer :: max_levels = 12  ! Maximum recursion depth
end type
```

---

## Key Simplifications from V2

### Removed Components

1. **NO separate thermionic emission model** - Everything through Tsu-Esaki
2. **NO Robin boundary conditions** - Direct current density BC
3. **NO n0B injection density** - Not needed with current BC
4. **NO surface recombination velocity** - Implicit in integral
5. **NO Richardson constant** - Absorbed in integral prefactor
6. **NO split TE/tunneling tracking** - Single unified current

### Clean Architecture

1. **Single physics module** - `schottky_tsuesaki_m`
2. **Single calculation function** - `calculate_schottky_current`
3. **Direct BC application** - No intermediate variables
4. **Minimal dependencies** - Only essential modules

---

## Numerical Implementation

### Integration Strategy

```fortran
! For tunneling region [0, φ_b]
! Use tanh-sinh quadrature (good for exponential decay)
call quad_tanhsinh(integrand_tunnel, 0.0, phi_b, ...)

! For thermionic region [φ_b, ∞]
! Use transformed variable u = exp(-E) → [exp(-φ_b), 1]
! Then E = -ln(u), dE = -du/u
call quad_gauss(integrand_therm_transformed, exp(-phi_b), 1.0, ...)
```

### Caching Strategy

```fortran
type :: current_cache
  real :: last_phi_b, last_F, last_phi_k
  real :: cached_J
  logical :: valid = .false.
end type

! Check cache before recalculating
if (cache_hit(phi_b, F, phi_k)) then
  J = cache%cached_J
else
  J = calculate_fresh(...)
  update_cache(...)
end if
```

---

## Validation Tests

### Test 1: Zero Field Limit
```
F → 0: J → J_thermionic_only
Should match Richardson equation
```

### Test 2: High Field Limit
```
F → ∞: J → J_field_emission
Should match Fowler-Nordheim
```

### Test 3: Energy Conservation
```
∫₀^∞ T(E)·N(E) dE = ∫₀^φ_b T(E)·N(E) dE + ∫_{φ_b}^∞ N(E) dE
```

### Test 4: Sign Convention
```
Electrons leaving semiconductor → J < 0
Verify with simple forward bias case
```

---

## Configuration Example

```ini
[contact]
name = "source"
type = "schottky"
phi_b = 0.7          : eV
m_tunnel = 0.4       : Effective mass ratio

[numerical]
rtol = 1e-6          : Integration tolerance
cache = true         : Enable current caching
```

---

## Physics Advantages

1. **Unified framework** - No artificial regime transitions
2. **Rigorous** - Captures all transport mechanisms naturally
3. **Self-consistent** - Fermi-Dirac statistics throughout
4. **Numerically stable** - Logarithmic occupancy formulation

## Implementation Advantages

1. **Simple** - Single calculation path
2. **Clean** - No legacy code or compatibility layers
3. **Fast** - Efficient quadrature with caching
4. **Maintainable** - Clear physics-to-code mapping

---

## Migration Path

1. **Phase 1**: Implement new `schottky_tsuesaki_m` module
2. **Phase 2**: Replace BC application in `continuity.f90`
3. **Phase 3**: Remove old `schottky.f90` and related code
4. **Phase 4**: Clean up device parameters (remove unused fields)

---

## References

1. Tsu & Esaki, "Tunneling in a finite superlattice," Appl. Phys. Lett., vol. 22, pp. 562-564, 1973.
2. ATLAS User Manual, "Schottky Barrier Models," Silvaco Inc.
3. Sze & Ng, "Physics of Semiconductor Devices," 3rd ed., Ch. 3, 2007.