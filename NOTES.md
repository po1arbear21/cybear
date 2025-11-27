# **NO FUCKING PLACEHOLDER ALLOWED**

# Schottky Contact V4: Pure Tsu-Esaki Implementation

**Author**: Chenyang Yu
**Date**: 2025-11-24
**Version**: 4.0 (Zero-Tolerance Implementation)

---
EVERY parameter comes from the ini file or the simulator, no hard coded constant!

## Core Principle: Direct Current Boundary Condition

The Schottky v4 implementation uses a **SINGLE** unified Tsu-Esaki integral that naturally handles both thermionic emission (over-barrier) and tunneling (under-barrier) transport. The current density is directly applied as a boundary condition in the continuity equation.

**ABSOLUTELY NO**:
- Placeholder functions or stub implementations
- "TODO" or "FIXME" comments
- Empty function bodies waiting for "future implementation"

---

## Physics Model: Pure Tsu-Esaki Integral

### The One and Only Current Formula

```
J_n = -(4πq·m*·kT/h³) ∫₀^∞ T(E) · ln[1 + exp((E_Fn - E)/kT)] dE

J_p = +(4πq·m*·kT/h³) ∫₀^∞ T(E) · ln[1 + exp((E - E_Fp)/kT)] dE
```

where:
- **T(E)**: WKB transmission coefficient (smooth transition at barrier height)
- **E_Fn, E_Fp**: Quasi-Fermi levels at the contact
- **Single integral**: Naturally captures all transport regimes

### WKB Transmission Coefficient

```
T(E) = exp(-2K)

where K = (2√(2m*)/ℏ) ∫ √(V(x) - E) dx

For triangular barrier approximation:
K = (4/3)√(2m*)/ℏF · (φ_b - E)^(3/2)  for E < φ_b
K = 0                                    for E ≥ φ_b
```
  subroutine wkb_transmission(E, phi_b, F, m_star, T, dT_dE, dT_dF)
    !! WKB transmission probability through triangular barrier
    real, intent(in) :: E, phi_b, F, m_star
    real, intent(out) :: T, dT_dE, dT_dF

    real :: F_smooth, gamma, coeff

    if (E >= phi_b) then
      T = 1.0
      dT_dE = 0.0
      dT_dF = 0.0
      return
    end if

    F_smooth = sqrt(F**2 + 1e-12**2)
    coeff = (4.0/3.0) * sqrt(2.0 * m_star)
    gamma = coeff * (phi_b - E)**1.5 / F_smooth

    T = exp(-gamma)
    dT_dE = T * coeff * 1.5 * sqrt(phi_b - E) / F_smooth
    dT_dF = T * coeff * (phi_b - E)**1.5 * F / F_smooth**3

  end subroutine

**Key**: T(E) = 1 for E ≥ φ_b (overbarrier) and T(E) < 1 for E < φ_b (tunneling)

### NO Barrier Lowering (Initially)

Version 4.0 implements the **clean physics first**:
- No image force barrier lowering (IFBL) in initial implementation
- Will be added as a separate, well-documented enhancement later
- φ_b remains constant for now

---

## Implementation Structure

### Main Blocks for `schottky.m`

```fortran
module schottky_m
  ! COMPLETE implementation - NO placeholders

  use constants_m
  use quadrature_m
  use device_m

  implicit none
  private

  public :: schottky_current_n
  public :: schottky_current_p

contains

  function schottky_current_n(phi_b, F, m_star, E_Fn, T) result(J_n)
    ! FULL implementation of electron current
    real(dp), intent(in) :: phi_b  ! Barrier height
    real(dp), intent(in) :: F      ! Electric field
    real(dp), intent(in) :: m_star ! Effective mass
    real(dp), intent(in) :: E_Fn   ! Electron quasi-Fermi level
    real(dp), intent(in) :: T      ! Temperature
    real(dp) :: J_n

    ! Actual quadrature integration
    J_n = integrate_tsu_esaki(phi_b, F, m_star, E_Fn, T, .true.)
  end function

  function schottky_current_p(phi_b, F, m_star, E_Fp, T) result(J_p)
    ! FULL implementation of hole current
    real(dp), intent(in) :: phi_b
    real(dp), intent(in) :: F
    real(dp), intent(in) :: m_star
    real(dp), intent(in) :: E_Fp   ! Hole quasi-Fermi level
    real(dp), intent(in) :: T
    real(dp) :: J_p

    ! Actual quadrature integration
    J_p = integrate_tsu_esaki(phi_b, F, m_star, E_Fp, T, .false.)
  end function

  function integrate_tsu_esaki(phi_b, F, m_star, E_F, T, is_electron) result(J)
    ! The REAL integration - NO shortcuts
    ! Adaptive quadrature over [0, ∞)
    ! Split at phi_b for numerical efficiency only
  end function

  function wkb_transmission(E, phi_b, F, m_star) result(T)
    ! COMPLETE WKB calculation
    real(dp) :: T
    if (E >= phi_b) then
      T = 1.0_dp
    else
      T = exp(-4.0_dp/3.0_dp * sqrt(2.0_dp*m_star*m_e) / (hbar * abs(F)) &
          * (phi_b - E)**1.5_dp)
    end if
  end function

end module
```

### Main Blocks for `continuity.m`

```fortran
module continuity_m
  use schottky_m

contains

  subroutine apply_boundary_conditions(...)
    ! Direct current BC application

    select case (contact_type)
    case (SCHOTTKY_CONTACT)
      ! Get the current density
      J_n = schottky_current_n(phi_b, F, m_n, E_Fn, T)
      J_p = schottky_current_p(phi_b, F, m_p, E_Fp, T)

      ! Apply as Neumann BC (current density)
      ! NO Robin BC, NO density calculations
      rhs_n = rhs_n - area * J_n / q
      rhs_p = rhs_p - area * J_p / q

    case (OHMIC_CONTACT)
      ! Standard Dirichlet BC

    end select

  end subroutine

  ! NO functions like:
  ! - calculate_n0B()  ← FORBIDDEN
  ! - calculate_p0B()  ← FORBIDDEN
  ! - robin_coefficient() ← FORBIDDEN
  ! - surface_recombination() ← FORBIDDEN

end module
```



## Numerical Implementation Details

### Integration Strategy

1. **Adaptive Quadrature**: Use Gauss-Kronrod or similar adaptive scheme
check quad.f90
ref : (REMOVE THE UNNECCESSARY INPUTS)
  subroutine tsu_esaki_current(phi_b, efield, m_tn, eta_semi_m, ci, E_g, J_tn, dJ_tn_dF, dJ_tn_deta)
    !! Calculate Tsu-Esaki tunneling current density through Schottky barrier
    !! Implements: J_tn = ±(m*/(2π²)) ∫₀^φ_eff T_wkb(E) N(E) dE
    !!
    !! Sign convention:
    !!   Electrons (ci=CR_ELEC): Negative prefactor (injection current)
    !!   Holes (ci=CR_HOLE):     Positive prefactor (extraction current)
    !!
    !! Barrier heights:
    !!   Electrons: phi_eff = phi_b (electron barrier to CB)
    !!   Holes:     phi_eff = E_g - phi_b (hole barrier to VB)
    !!
    !! This is the tunneling component only (under-barrier transport)
    !! Total current: J_total = J_TE + J_tn (ATLAS split model)
    !!
    !! All inputs/outputs in normalized units (kT/q normalization)

    real, intent(in)  :: phi_b      !! Schottky barrier height for electrons (with IFBL if applied)
    real, intent(in)  :: efield     !! Electric field magnitude (normalized)
    real, intent(in)  :: m_tn       !! Tunneling effective mass ratio (m*/m0)
    real, intent(in)  :: eta_semi_m !! Semiconductor quasi-Fermi relative to metal Fermi (constant)
    integer, intent(in) :: ci       !! Carrier index (CR_ELEC or CR_HOLE)
    real, intent(in)  :: E_g        !! Band gap (normalized)
    real, intent(out) :: J_tn       !! Tunneling current density (normalized)
    real, intent(out) :: dJ_tn_dF   !! Derivative dJ_tn/dF for Jacobian
    real, intent(out) :: dJ_tn_deta !! Derivative dJ_tn/deta_semi_m for density coupling

    real :: p(6)                   ! Parameter array for integrand
    real :: I                      ! Integration result
    real :: dIda, dIdb            ! Derivatives wrt bounds (not used)
    real :: dIdp(6)               ! Derivatives wrt parameters
    real :: err                   ! Integration error estimate
    integer :: ncalls             ! Number of function evaluations
    real :: prefactor             ! Normalization prefactor: ±m*/(2π²)
    real :: phi_eff               ! Effective barrier for this carrier

    ! Calculate effective barrier based on carrier type
    if (ci == CR_ELEC) then
      phi_eff = phi_b              ! Electrons: phi_bn (CB barrier)
    else  ! CR_HOLE
      phi_eff = E_g - phi_b        ! Holes: phi_bp = E_g - phi_bn (VB barrier)
    end if

    ! Check for degenerate cases
    if (phi_eff <= 0.0) then
      ! No barrier - no tunneling calculation needed
      J_tn = 0.0
      dJ_tn_dF = 0.0
      dJ_tn_deta = 0.0
      return
    end if

    if (abs(efield) < 1e-10) then
      ! Zero field - no tunneling possible
      J_tn = 0.0
      dJ_tn_dF = 0.0
      dJ_tn_deta = 0.0
      return
    end if

    ! Set up parameter array for integrand
    ! p(1) = phi_eff    : effective barrier (carrier-dependent)
    ! p(2) = efield     : electric field
    ! p(3) = m_tn       : tunneling mass
    ! p(4) = eta_semi_m : semiconductor quasi-Fermi relative to metal Fermi
    ! p(5) = ci         : carrier index
    ! p(6) = E_g        : band gap
    p(1) = phi_eff
    p(2) = efield
    p(3) = m_tn
    p(4) = eta_semi_m
    p(5) = real(ci)
    p(6) = E_g

    ! Perform integration (under-barrier only)
    ! Using adaptive quadrature with relative tolerance 1e-6
    if (ci == CR_ELEC) then
       call quad(tsu_esaki_integrand, 0.0, phi_eff, p, I, dIda, dIdb, dIdp, &
                rtol=1e-6, err=err, max_levels=12, ncalls=ncalls)
    else  ! CR_HOLE
       call quad(tsu_esaki_integrand, -phi_eff, 0.0, p, I, dIda, dIdb, dIdp, &
                rtol=1e-6, err=err, max_levels=12, ncalls=ncalls)
    end if

    ! call quad(tsu_esaki_integrand, 0.0, phi_eff, p, I, dIda, dIdb, dIdp, &
    !             rtol=1e-6, err=err, max_levels=12, ncalls=ncalls)

    ! Apply normalization prefactor with carrier-dependent sign
    ! Electrons: negative (injection)
    ! Holes: positive (extraction)
    if (ci == CR_ELEC) then
      prefactor = m_tn / (2.0 * PI**2)  ! Negative for electron injection
    else  ! CR_HOLE
      prefactor = m_tn / (2.0 * PI**2)  ! Positive for hole extraction
    end if

    ! Calculate tunneling current and derivatives
    J_tn = prefactor * I
    dJ_tn_dF = prefactor * dIdp(2)   ! dI/d(efield)
    dJ_tn_deta = prefactor * dIdp(4) ! dI/d(eta_semi_m)

  end subroutine


2. **Energy Grid**: (optional)
   - Dense near E = φ_b (transition region)
   - Logarithmic spacing for E >> φ_b (to phib + a few kT )
4. **No Caching Initially**: Get it working correctly first

### Field Extraction at Contact

```fortran
! Get field from potential gradient

REF (think carefully could be wrong):
  function get_schottky_contact_normal_dir(par, ict) result(normal_dir)
    !! Determine the normal direction of a Schottky contact for barrier lowering
    !! Analyzes vertex coordinates to find the constant dimension
    !! Returns: 1 for x-normal, 2 for y-normal, 3 for z-normal
    !! Returns: 0 if contact is not axis-aligned

    type(device_params), intent(in) :: par
    integer,             intent(in) :: ict     ! Contact index
    integer                         :: normal_dir

    integer :: i, dir, idx(par%g%idx_dim)
    integer :: n_vertices
    real    :: coord_min(3), coord_max(3), coord_range(3)
    real    :: vertex_pos(3)
    real    :: tolerance, min_range, char_length

    ! Initialize
    normal_dir = 0
    coord_min = 1e30
    coord_max = -1e30
    tolerance = 1e-6  ! Relative tolerance for planarity

    ! Get number of transport vertices for this contact
    n_vertices = par%transport_vct(ict)%n

    if (n_vertices == 0) then
      print *, "WARNING: Schottky contact ", ict, " has no transport vertices"
      return
    end if

    ! Analyze coordinate ranges
    do i = 1, n_vertices
      idx = par%transport_vct(ict)%get_idx(i)
      call par%g%get_vertex(idx, vertex_pos(1:par%g%dim))

      do dir = 1, par%g%dim
        coord_min(dir) = min(coord_min(dir), vertex_pos(dir))
        coord_max(dir) = max(coord_max(dir), vertex_pos(dir))
      end do
    end do

    ! Calculate range in each direction
    do dir = 1, par%g%dim
      coord_range(dir) = coord_max(dir) - coord_min(dir)
    end do

    ! Find direction with minimum range (normal to contact)
    min_range = huge(1.0)
    do dir = 1, par%g%dim
      if (coord_range(dir) < min_range) then
        min_range = coord_range(dir)
        normal_dir = dir
      end if
    end do

    ! Validate planarity
    if (par%g%dim > 1) then
      char_length = maxval(coord_range(1:par%g%dim))
      if (min_range > tolerance * char_length .and. char_length > 0.0) then
        print *, "WARNING: Schottky contact ", ict, " not perfectly axis-aligned for barrier lowering"
        print *, "  Range in normal direction: ", denorm(min_range, "nm"), " nm"
      end if
    end if

  end function


and then:
E_field = this%efield(normal_dir)%get(idx)



### Quasi-Fermi Level Calculation (check poisson)

```fortran
! From carrier density and potential
E_Fn = E_c - phi + kT * log(n/N_c)
E_Fp = E_v - phi - kT * log(p/N_v)

! These MUST be actual calculations, not zeros
```

---


## Configuration

```toml
[schottky_test]
schottky_diode.ini
```

---

## Implementation Checklist

- [ ] Complete `schottky_current()` with REAL integration
- [ ] Implement WKB transmission (no stubs)
- [ ] Implement field extraction at contacts
---

## Future Enhancements (After v4.0 Works)

1. **Image Force Barrier Lowering**: Add Schottky barrier lowering
2. **Anisotropic masses**: Different m* for tunneling vs transport
3. **Non-parabolic bands**: Energy-dependent effective mass
4. **Performance optimization**: Caching, lookup tables

---

## Schottky Jacobian Implementation Plan

### Problem
The lagged approach for Schottky current has slow convergence because the Jacobian doesn't include `dJ_sch/dn`. The Newton-Raphson solver needs this derivative for quadratic convergence.

### FVM Residual at Schottky Boundary
```
f = Σ(J_dd·A_ij) - J_sch·A_ct = 0
```

### Required Jacobian
```
df/dn = -dJ_sch/dn · A_ct
```

### Chain Rule for dJ_sch/dn
The Schottky current depends on density through the quasi-Fermi level:
- `eta = ln(n/Nc)` → `deta/dn = 1/n`
- For electrons: `eta_m = eta + phi_b` → `d(eta_m)/dn = 1/n`
- For holes: `eta_m = -eta - phi_b` → `d(eta_m)/dn = -1/n`

Therefore:
```
dJ_sch/dn = dJ_sch/d(eta_m) · d(eta_m)/dn
         = dJ_sch/d(eta_m) · (±1/n)
```

The quad integration already computes `dIdp(4) = dI/d(eta_m)`, so we can get `dJ_sch/d(eta_m)` directly.

### Implementation Steps

#### Step 1: Modify schottky_m

Add subroutine that returns both J and dJ/dn:

```fortran
subroutine schottky_current_with_deriv(par, ict, ci, E_normal, n_dens, J, dJ_dn)
  ! Same computation as schottky_current, but also returns derivative

  ! After quad integration:
  ! dJ/d(eta_m) = prefactor * dIdp(4) (with carrier sign)

  ! Chain rule:
  if (ci == CR_ELEC) then
    dJ_dn = prefactor * dIdp(4) / n_dens
  else
    dJ_dn = prefactor * dIdp(4) / n_dens  ! sign handled in prefactor
  end if
end subroutine
```

**Export:** Add `schottky_current_with_deriv` to public interface.

#### Step 2: Modify continuity_m

**2a. Add Jacobian pointer**
```fortran
type(jacobian), pointer :: jaco_schottky => null()
  !! Jacobian for Schottky current wrt density (non-const)
```

**2b. In init - Create Jacobian for Schottky contacts**
```fortran
! Only if Schottky contacts exist
if (any(par%contacts(1:par%nct)%type == CT_SCHOTTKY)) then
  this%jaco_schottky => this%init_jaco_f(idens, &
    & st = [this%st_em%get_ptr(), (st_schottky_ct(ict), ict = 1, par%nct)], &
    & const = .false., dtime = .false.)
end if
```

Where `st_schottky_ct(ict)` is `st_dir` for Schottky contacts, `st_em` otherwise.

**2c. In eval - Compute and set Jacobian entries**
```fortran
if (associated(this%jaco_schottky)) then
  call this%jaco_schottky%reset()

  ! Loop over Schottky contacts
  do ict = ...
    if (type == CT_SCHOTTKY) then
      ! Get J and dJ/dn
      call schottky_current_with_deriv(par, ict, ci, E_normal, n_dens, J_sch, dJ_dn)

      ! Set Jacobian entry: df/dn = -dJ_sch/dn * A_ct
      call this%jaco_schottky%set(idx, idx, -dJ_dn * A_ct)

      ! Add to residual (same as before)
      tmp(j) = tmp(j) - J_sch * A_ct
    end if
  end do

  call this%jaco_schottky%set_matr(const = .false., nonconst = .true.)
end if
```

### Files to Modify

1. **src/schottky.f90**
   - Add `schottky_current_with_deriv` subroutine
   - Export in public interface

2. **src/continuity.f90**
   - Add `jaco_schottky` pointer to type
   - Create Jacobian in init (with st_dir for Schottky, st_em otherwise)
   - Compute dJ/dn and set entries in eval
   - Call `set_matr` to materialize

### Key Considerations

1. **Sign convention:** Ensure consistent signs for electrons vs holes
2. **Numerical stability:** Handle small n (avoid division by zero)
3. **Stencil structure:** jaco_schottky uses diagonal stencil (st_dir) only for Schottky contacts
4. **const flag:** Must be `const = .false.` since entries change each iteration

---

## Convergence Failure Investigation & Fix Plan (2025-11-27)

### Problem Summary

**Observed Behavior:**
- Reverse bias: ✓ Converges perfectly
- Forward bias ≤0.191V: ✓ Converges
- Forward bias ≥0.192V: ✗ Newton solver fails (max iterations exceeded)
- **Jacobian verification results:**
  - At V=0.192V: "Jacobian may UNDERESTIMATE derivatives" (min ||F|| at t=1.15)
  - At V=1.0V: "Jacobian may OVERESTIMATE derivatives" (min ||F|| at t=0.35)
- Loosening rtol to 1e-6 allows convergence with physically reasonable results

### Root Cause Analysis

#### Finding 1: Integration Tolerance is Too Loose

**Location:** `src/schottky.f90:100`
```fortran
call quad(tsu_esaki_integrand, E_min, E_max, params, integral, &
          dIda, dIdb, dIdp, rtol=1.0e-6, err=err, ncalls=ncalls)
```

**Problem:**
- The quadrature routine computes BOTH the integral AND parameter derivatives simultaneously
- Parameter derivative `dIdp(4) = dI/d(eta_m)` is used directly in Jacobian computation
- At high forward bias (V>0.19V), the integrand develops sharp features:
  - Fermi-Dirac distributions become step-like
  - WKB transmission varies rapidly near E=phi_b
  - Integrand spans 10+ orders of magnitude

**Impact:**
- Integration error in dIdp(4) → error in dJ/ddens → incorrect Jacobian
- Accumulated errors prevent Newton solver from finding quadratic convergence region
- Error magnitude is bias-dependent, explaining regime-dependent Jacobian behavior

**Evidence:**
- Jacobian underestimates at 0.192V: derivatives too small → integration missed sharp features
- Loosening Newton rtol works: solver tolerates inaccurate Jacobian and converges slowly

#### Finding 2: WKB Discontinuity at E=phi_b

**Location:** `src/schottky.f90:224-229`
```fortran
if (E >= phi_b) then
  T = 1.0
  dT_dE = 0.0    ! ← Jump discontinuity
  dT_dphi = 0.0
  dT_dF = 0.0
  return
end if
```

**Problem:**
- Hard discontinuity in transmission derivatives at E=phi_b
- For phi_b=0.7eV (27kT), at V=0.19V:
  - Integration window: [0, phi_b+20] = [0, 47kT]
  - Window straddles the discontinuity
  - Adaptive quadrature struggles with discontinuous derivatives

**Impact:**
- Reduces integration accuracy near the discontinuity
- Parameter derivatives (dT_dphi, dT_dF) jump discontinuously
- Exacerbates errors in dIdp computation

#### Finding 3: The 0.192V Threshold is NOT Magic

**Physical interpretation:**
- At V=0.192V (~7.5kT): Current J ∝ exp(V/kT) ≈ 1800× equilibrium
- Jacobian dJ/ddens grows by same exponential factor
- Integration error that was negligible at low bias becomes significant
- Newton solver enters difficult regime where:
  - Residuals are ~1e-18 (very small)
  - Relative tolerance requirement is 1e-9
  - Any Jacobian inaccuracy prevents convergence

### Solution Plan

#### Phase 1: Tighten Integration Tolerance (PRIMARY FIX)

**Priority:** CRITICAL - Most likely root cause

**Implementation:**

Change `src/schottky.f90:100`:
```fortran
! OLD:
call quad(tsu_esaki_integrand, E_min, E_max, params, integral, &
          dIda, dIdb, dIdp, rtol=1.0e-6, err=err, ncalls=ncalls)

! NEW:
call quad(tsu_esaki_integrand, E_min, E_max, params, integral, &
          dIda, dIdb, dIdp, rtol=1.0e-8, err=err, ncalls=ncalls)
```

**Rationale:**
- 1e-8 matches Newton solver rtol=1e-9 requirement
- Ensures parameter derivatives accurate enough for Jacobian
- Minimal code change, easy to test and revert

**Testing:**
1. Run V_SCHOTTKY = 0.192V → should converge in <15 iterations
2. Run Jacobian verification → min ||F|| should occur near t=1.0
3. Check convergence across range: 0.18V, 0.19V, 0.195V, 0.2V, 0.5V, 1.0V
4. Monitor ncalls increase (expect 2-3× computational cost)

**Expected outcome:**
- ✓ Convergence at 0.192V and above
- ✓ Jacobian verification passes (min at t ≈ 0.95-1.05)
- ✓ Integration error drops below 1e-9
- ✓ J-V curve matches exp(qV/kT) behavior

#### Phase 2: Smooth WKB Discontinuity (IF NEEDED)

**Priority:** MEDIUM - Secondary fix if Phase 1 insufficient

**Implementation:**

Replace hard discontinuity with tanh smoothing in `src/schottky.f90:213-243`:

```fortran
pure subroutine wkb_transmission(E, phi_b, F, m_star, T, dT_dE, dT_dphi, dT_dF)
  real, intent(in)  :: E, phi_b, F, m_star
  real, intent(out) :: T, dT_dE, dT_dphi, dT_dF

  real :: gamma, dE, F_safe, coeff
  real :: smooth, dsmooth_dE, exp_gamma
  real, parameter :: eps = 1e-10
  real, parameter :: delta = 0.05  ! Smoothing width (kT units)

  dE = phi_b - E
  F_safe = sqrt(F**2 + eps**2)
  coeff = (4.0/3.0) * sqrt(2.0 * m_star)

  if (dE < -5.0*delta) then
    ! Far above barrier: pure thermionic (E >> phi_b)
    T = 1.0
    dT_dE = 0.0
    dT_dphi = 0.0
    dT_dF = 0.0
  else if (dE > 5.0*delta) then
    ! Far below barrier: pure tunneling (E << phi_b)
    gamma = coeff * dE**1.5 / F_safe
    exp_gamma = exp(-gamma)

    T = exp_gamma
    dT_dE   = exp_gamma * coeff * 1.5 * sqrt(dE) / F_safe
    dT_dphi = -dT_dE
    dT_dF   = exp_gamma * coeff * dE**1.5 * F / (F_safe**3)
  else
    ! Transition region: smooth interpolation
    ! smooth(dE) = 1 / (1 + exp(-dE/delta))
    smooth = 1.0 / (1.0 + exp(-dE/delta))
    dsmooth_dE = smooth * (1.0 - smooth) / delta

    gamma = coeff * max(dE, 0.0)**1.5 / F_safe
    exp_gamma = exp(-gamma)

    ! T interpolates: tunneling → thermionic
    T = 1.0 - smooth * (1.0 - exp_gamma)

    ! Derivatives with chain rule
    if (dE > 0.0) then
      dT_dE = -dsmooth_dE * (1.0 - exp_gamma) &
              + smooth * exp_gamma * coeff * 1.5 * sqrt(dE) / F_safe
      dT_dphi = -dT_dE
      dT_dF = smooth * exp_gamma * coeff * dE**1.5 * F / (F_safe**3)
    else
      dT_dE = -dsmooth_dE * (1.0 - exp_gamma)
      dT_dphi = -dT_dE
      dT_dF = 0.0
    end if
  end if
end subroutine wkb_transmission
```

**Rationale:**
- Eliminates derivative discontinuities
- delta=0.05kT gives narrow transition (~0.0013eV at 300K)
- Physically motivated: real barriers have smooth transitions
- Improves numerical integration robustness

**Testing:**
1. Compare J-V curves: new vs old implementation
2. Check difference is <1% for V<0.6V (should be negligible)
3. Verify convergence across full bias range
4. Jacobian verification should improve

**Expected outcome:**
- Robust convergence independent of bias point
- Smoother Jacobian behavior
- Slight change in J-V near barrier top (E≈phi_b region)

#### Phase 3: Add Diagnostic Output (INFRASTRUCTURE)

**Priority:** HIGH - Needed for debugging

**Implementation:**

Uncomment and enhance debug output in `src/schottky.f90`:

```fortran
! After line 100:
if (err > 1.0e-9) then
  print "(A,I2,A,ES12.4,A,I6)", &
    "  [SCHOTTKY WARNING] Contact ", ict, " err=", err, " ncalls=", ncalls
end if

! Uncomment lines 94-96, 120-122 for detailed diagnostics
```

Add to `src/continuity.f90` (uncomment lines 360-368):
```fortran
print "(A,I3,A,I3)", "DEBUG_SCHOTTKY: contact=", ict, " vertex=", i
print "(A,2ES14.6)", "  dens, DOS = ", dens, this%par%smc%edos(this%ci)
print "(A,2ES14.6)", "  E_normal, phi_b = ", E_normal, this%par%contacts(ict)%phi_b
print "(A,2ES14.6)", "  J_sch, A_ct = ", J_sch, A_ct
print "(A,ES14.6)", "  dJ_ddens = ", dJ_ddens
```

**Testing:**
- Run with diagnostics at V=0.192V
- Check for:
  - Integration errors >1e-9
  - Abnormal dJ_ddens values
  - Vertices with problematic Jacobian entries

#### Phase 4: Validation Tests

**Test suite to run after each fix:**

1. **Convergence test:**
   ```ini
   V_SCHOTTKY = 0.18, 0.185, 0.19, 0.191, 0.192, 0.195, 0.2, 0.5, 1.0 : V
   ```
   - All should converge in <20 iterations

2. **Jacobian verification:**
   - Run built-in test at each voltage
   - Check: min ||F|| occurs at t ∈ [0.95, 1.05]

3. **Physics validation:**
   - Plot log(J) vs V → should be linear with slope ≈ q/kT
   - Check current density magnitude matches expected ~A/cm² range

4. **Reverse bias check:**
   - V_SCHOTTKY = -0.5, -1.0V → ensure still works

5. **Computational cost:**
   - Compare total ncalls and wall time
   - Acceptable if <5× increase

### Critical Files

1. **src/schottky.f90**
   - Line 100: Integration tolerance (Phase 1)
   - Lines 213-243: WKB transmission (Phase 2)
   - Lines 94-96, 120-122: Diagnostics (Phase 3)

2. **src/continuity.f90**
   - Lines 360-368: Debug output (Phase 3)

3. **run_schottky.ini**
   - Add voltage sweep for systematic testing

### Implementation Order

**STEP 1:** Enable diagnostics (Phase 3) - 5 minutes
- Uncomment debug lines
- Run at V=0.192V
- Confirm integration error and Jacobian behavior

**STEP 2:** Apply Phase 1 fix - 2 minutes
- Change rtol to 1e-8
- Test convergence at 0.192V
- Run Jacobian verification
- **Decision point:** If converges → SUCCESS, cleanup diagnostics
- If not → proceed to Step 3

**STEP 3:** Apply Phase 2 fix (if needed) - 30 minutes
- Implement WKB smoothing
- Validate against low-bias reference
- Test full bias range
- Compare J-V curves

**STEP 4:** Full validation - 15 minutes
- Run complete test suite
- Document results
- Remove diagnostic output

### Success Criteria

**Must achieve:**
- ✓ Convergence at V=0.192V and all higher forward bias
- ✓ Jacobian verification passes (min at t≈1.0)
- ✓ J-V curve follows exp(qV/kT) behavior

**Should achieve:**
- ✓ Convergence in <20 Newton iterations
- ✓ No warnings about integration errors
- ✓ Computational cost <3× original

**Nice to have:**
- ✓ Robust convergence for any reasonable bias
- ✓ Clean diagnostic output (no anomalies)
