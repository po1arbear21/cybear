# Schottky Contact Implementation Specification

## Overview
Complete specification for Thermionic Emission (TE) and Tsu-Esaki tunneling implementation in the Cybear drift-diffusion simulator. This document enables duplication/optimization in another branch.

---

## NOTES ON IMPLEMENTATION

### 1. **TEST ARTIFACT**: Hardcoded barrier lowering (schottky.f90:379)
```fortran
phi_b = phi_b - norm(0.2, "eV")  ! INTENTIONAL TEST - remove for production
```
This 0.2 eV barrier reduction is **intentional for testing**. Remove when deploying to production.

### 2. Debug print statements to remove/comment for production:
- `schottky.f90:241-245` - Active debug prints in `schottky_velocity`
- `contact.f90:190-193` - Debug prints in `contact_set_phims_schottky`
- `region.f90:175,183-184` - Debug prints during config parsing

### 3. Hole sign convention (schottky.f90:387-391) - **VERIFIED CORRECT**
```fortran
! Both carriers use same formula - this IS correct:
if (ci == CR_ELEC) then
  eta_m = eta + phi_b
else
  eta_m = eta + phi_b  ! Correct! phi_b is already adjusted for holes at line 363
end if
```
**Why this is correct:**
- For electrons: φ_b = electron barrier height (from config)
- For holes: φ_b = E_g - electron_barrier (computed at line 363)
- η from `get_idist` is consistently defined for both carriers
- η_m = η + φ_b correctly positions quasi-Fermi relative to metal Fermi level

---

## FILE STRUCTURE

| File | Purpose |
|------|---------|
| `src/schottky.f90` | Core Schottky physics (TE velocity, n0B, tunneling) |
| `src/contact.f90` | Contact data structure + phims calculation |
| `src/continuity.f90` | Robin BC application in continuity equation |
| `src/region.f90` | Config parsing for Schottky parameters |
| `src/device_params.f90` | Contact initialization + surface area calc |

---

## CORE DATA STRUCTURES

### Contact Type (contact.f90:24-50)
```fortran
type contact
  character(:), allocatable :: name
  integer                   :: type          ! CT_OHMIC=1, CT_GATE=2, CT_SCHOTTKY=3
  real                      :: phims         ! Metal-semiconductor workfunction diff
  real                      :: phi_b         ! Barrier height (normalized to kT)
  real                      :: A_richardson_n = 112.0  ! A/cm²/K² for electrons
  real                      :: A_richardson_p = 32.0   ! A/cm²/K² for holes
  logical                   :: ifbl = .false.          ! Image force barrier lowering
  logical                   :: tunneling = .false.     ! Enable Tsu-Esaki tunneling
  real                      :: m_tunnel_n             ! Tunneling mass ratio (m*/m0) electrons
  real                      :: m_tunnel_p             ! Tunneling mass ratio (m*/m0) holes
contains
  procedure :: set_phims_schottky => contact_set_phims_schottky
end type
```

### Continuity Schottky Support (continuity.f90:42-52)
```fortran
! Additional members for Schottky:
type(vselector), allocatable :: efield(:)      ! E-field components
integer, allocatable :: schottky_normal(:)     ! Normal direction per contact
real, allocatable :: v_surf(:)                 ! Thermionic velocity per contact
real, allocatable :: accum_TE(:), accum_TN(:)  ! Accumulated currents (debug)
```

---

## CORE PHYSICS FUNCTIONS

### 1. Thermionic Velocity (schottky.f90:212-247)
```fortran
function schottky_velocity(par, ci, ict) result(v_th)
  !! v_th = A* T² / N_c  (Richardson velocity for Robin BC)

  ! Get Richardson constant based on carrier type
  if (ci == CR_ELEC) then
    A_star = par%contacts(ict)%A_richardson_n
  else
    A_star = par%contacts(ict)%A_richardson_p
  end if

  ! Compute velocity (normalized units)
  v_th = A_star * norm(par%T, "K")**2 / par%smc%edos(ci)
end function
```

### 2. Equilibrium Injection Density with IFBL (schottky.f90:252-316)
```fortran
subroutine schottky_n0b(par, ci, ict, E_normal, n0B, dn0B_dE)
  !! n0B = N_c exp(-φ_b_eff / kT)
  !! With IFBL: φ_b_eff = φ_b - sqrt(|E| / (4π))

  ! Get barrier height (normalized to kT)
  if (ci == CR_ELEC) then
    phi_b = par%contacts(ict)%phi_b
  else
    phi_b = par%smc%band_gap - par%contacts(ict)%phi_b
  end if
  phi_b_eff = phi_b

  ! Apply IFBL if enabled
  if (par%contacts(ict)%ifbl) then
    delta_phi = sqrt(abs(E_normal) / (4.0 * PI))
    phi_b_eff = phi_b - delta_phi

    ! Derivative for Jacobian: d(delta_phi)/dE = delta_phi/(2E)
    if (present(dn0B_dE)) then
      ddelta_dE = delta_phi / (2.0 * abs(E_normal))
    end if
  end if

  ! Equilibrium density
  n0B = par%smc%edos(ci) * exp(-phi_b_eff)

  ! IFBL derivative for Newton coupling
  if (present(dn0B_dE)) then
    dn0B_dE = n0B * ddelta_dE
  end if
end subroutine
```

### 3. Tunneling Current via Tsu-Esaki Integral (schottky.f90:321-419)
```fortran
subroutine schottky_tunneling(par, ci, ict, E_normal, dens, J_tn, dJ_tn_ddens)
  !! J_tn = (m*/2π²) ∫₀^{φ_b} T(E) · supply(E) dE

  ! Skip if tunneling disabled
  if (.not. par%contacts(ict)%tunneling) then
    J_tn = 0.0; dJ_tn_ddens = 0.0; return
  end if

  ! Get barrier and tunneling mass
  if (ci == CR_ELEC) then
    phi_b = par%contacts(ict)%phi_b
    m_tunnel = par%contacts(ict)%m_tunnel_n
  else
    phi_b = par%smc%band_gap - par%contacts(ict)%phi_b
    m_tunnel = par%contacts(ict)%m_tunnel_p
  end if

  ! Prefactor: m* / (2π²)
  prefactor = m_tunnel / (2.0 * PI**2)

  ! Apply IFBL to barrier
  if (par%contacts(ict)%ifbl) then
    phi_b = phi_b - sqrt(abs(E_normal) / (4.0 * PI))
  end if

  ! Compute η (quasi-Fermi level relative to band edge)
  call par%smc%get_idist(dens / par%smc%edos(ci), eta, detadF)

  ! η relative to metal Fermi level
  eta_m = eta + phi_b

  ! Integration parameters: [phi_b, |E|, m*, eta_m]
  params = [phi_b, E_normal, m_tunnel, eta_m]

  ! Integrate from 0 to phi_b (sub-barrier tunneling only)
  call quad(tsu_esaki_integrand, 0.0, phi_b, params, integral, &
            dIda, dIdb, dIdp, rtol=1.0e-9, err=err, ncalls=ncalls)

  ! Current density
  J_tn = -prefactor * integral

  ! Jacobian derivative via chain rule
  dJ_tn_ddens = -prefactor * dIdp(4) * detadF / par%smc%edos(ci)
end subroutine
```

### 4. Tsu-Esaki Integrand (schottky.f90:424-463)
```fortran
subroutine tsu_esaki_integrand(E, p, f, dfdE, dfdp)
  !! f(E) = T(E) * [ln(1+exp(η_m-E)) - ln(1+exp(-E))]

  phi_b  = p(1); efield = p(2); m_star = p(3); eta_m = p(4)

  ! WKB transmission coefficient
  call wkb_transmission(E, phi_b, efield, m_star, T, dT_dE, dT_dphi, dT_dF)

  ! Fermi-Dirac distributions
  f_metal = 1.0 / (1.0 + exp(E))
  f_semi  = 1.0 / (1.0 + exp(E - eta_m))
  dN_dE   = f_metal - f_semi

  ! Supply function (numerically stable form)
  N_supply = log1p(expm1(eta_m) * f_metal)

  ! Integrand and derivatives
  f = T * N_supply
  dfdE = dT_dE * N_supply + T * dN_dE
  dfdp(1) = dT_dphi * N_supply
  dfdp(2) = dT_dF * N_supply
  dfdp(3) = 0.0
  dfdp(4) = T * f_semi
end subroutine
```

### 5. WKB Transmission Coefficient (schottky.f90:468-501)
```fortran
pure subroutine wkb_transmission(E, phi_b, F, m_star, T, dT_dE, dT_dphi, dT_dF)
  !! T = exp(-γ) for E < phi_b, T = 1 for E >= phi_b
  !! γ = (4/3) * sqrt(2m*) * (phi_b - E)^(3/2) / F

  real, parameter :: eps = 1e-10

  dE = phi_b - E
  F_safe = sqrt(F**2 + eps**2)  ! Avoid division by zero
  coeff = (4.0/3.0) * sqrt(2.0 * m_star)

  if (dE < 0.0) then
    ! Above barrier: thermionic regime
    T = 1.0; dT_dE = 0.0; dT_dphi = 0.0; dT_dF = 0.0
  else
    ! Below barrier: tunneling regime
    gamma = coeff * dE**1.5 / F_safe
    exp_gamma = exp(-gamma)

    T = exp_gamma
    dT_dE   = exp_gamma * coeff * 1.5 * sqrt(dE) / F_safe
    dT_dphi = -dT_dE
    dT_dF   = exp_gamma * coeff * dE**1.5 * F / (F_safe**3)
  end if
end subroutine
```

---

## BOUNDARY CONDITION APPLICATION

### Robin BC in Continuity Equation (continuity.f90:340-441)

**Equation Form:**
```
F = div(J) - A_ct * [v_surf*(n - n0B) + J_tn] = 0
```

**Jacobian Entry:**
```
dF/d(dens) = v_surf*A_ct - A_ct*dJ_tn/ddens
```

**Key Code (continuity_eval):**
```fortran
! For each Schottky contact vertex:
if (par%contacts(ict)%type == CT_SCHOTTKY) then
  ! Get E-field and density at this vertex
  E_normal = efield_arr(j)
  dens = dens_arr(j)

  ! Get contact surface area and thermionic velocity
  A_ct = this%par%get_ct_surf(ict, idx)
  v_surf_ict = this%v_surf(ict)

  ! Compute n0B (with optional IFBL)
  call schottky_n0b(this%par, this%ci, ict, E_normal, n0B)

  ! Compute tunneling current and derivative
  call schottky_tunneling(this%par, this%ci, ict, E_normal, dens, J_tn, dJ_tn_ddens)

  ! Set Jacobian: dF/d(dens) = v_surf*A_ct - A_ct*dJ_tn/ddens
  call this%jaco_dens%set(idx, idx, v_surf_ict * A_ct - A_ct * dJ_tn_ddens)

  ! Add residual: F = div(J) + v_surf*A_ct*dens - v_surf*A_ct*n0B - A_ct*J_tn
  tmp(j) = tmp(j) + v_surf_ict * A_ct * dens - v_surf_ict * A_ct * n0B - A_ct * J_tn
end if

! Materialize non-constant Jacobian
call this%jaco_dens%set_matr(const = .false., nonconst = .true.)
```

---

## CONTACT SURFACE AREA CALCULATION

### get_ct_surf (device_params.f90:1362-1452)
```fortran
function get_ct_surf(this, ict, idx_vert) result(ct_surf)
  select case (this%g%dim)
  case (1)
    ct_surf = 1.0  ! Normalized cross-sectional area

  case (2)
    ! Sum half of edge lengths to neighboring contact vertices
    ct_surf = 0.0
    do idx_dir = 1, 2
      ! Get neighbor in this direction
      call this%g%get_neighb(..., idx_neighb, status)
      if (status .and. this%ict%get(idx_neighb) == ict) then
        edge_len = sqrt((p2(1)-p1(1))**2 + (p2(2)-p1(2))**2)
        ct_surf = ct_surf + 0.5 * edge_len
      end if
    end do

  case (3)
    ! Sum 1/4 of contacted face areas
    ct_surf = 0.0
    do idx_dir = 1, 3
      do m = 1, 2
        ! Check if face is fully contacted (all 4 vertices)
        if (all_contacted) then
          surf_area = 0.5 * norm(cross(vec1, vec2))
          ct_surf = ct_surf + 0.25 * surf_area
        end if
      end do
    end do
  end select
end function
```

---

## PHIMS CALCULATION FOR SCHOTTKY (contact.f90:169-197)

```fortran
subroutine contact_set_phims_schottky(this, smc)
  !! phims = -phi_b + ln(N_c/n_i) = -phi_b + 0.5*ln(N_c/N_v) + E_g/2

  ni_term = 0.5 * log(smc%edos(CR_ELEC) / smc%edos(CR_HOLE)) + 0.5 * smc%band_gap
  this%phims = -this%phi_b + ni_term
end subroutine
```

---

## CONFIG PARSING (region.f90:166-184)

```fortran
! Parse Schottky-specific parameters
call file%get(sid, "phi_b", this%phi_b, status = st)
call file%get(sid, "A_richardson_n", this%A_richardson_n, status = st)
if (.not. st) this%A_richardson_n = norm(112.0, "A/cm^2/K^2")  ! Default Si
call file%get(sid, "A_richardson_p", this%A_richardson_p, status = st)
if (.not. st) this%A_richardson_p = norm(32.0, "A/cm^2/K^2")   ! Default Si
call file%get(sid, "ifbl", this%ifbl, status = st)
if (.not. st) this%ifbl = .false.
call file%get(sid, "tunneling", this%tunneling, status = st)
if (.not. st) this%tunneling = .false.
call file%get(sid, "m_tunnel_n", this%m_tunnel_n, status = st)
if (.not. st) this%m_tunnel_n = 1.0
call file%get(sid, "m_tunnel_p", this%m_tunnel_p, status = st)
if (.not. st) this%m_tunnel_p = 1.0
```

---

## INITIALIZATION SEQUENCE

### 1. Region Parsing (region.f90)
Config file → region_contact structure

### 2. Contact Initialization (device_params.f90)
```fortran
! Transfer Schottky parameters from region to contact
contacts(ict)%phi_b = reg_ct(ri)%phi_b
contacts(ict)%A_richardson_n = reg_ct(ri)%A_richardson_n
contacts(ict)%A_richardson_p = reg_ct(ri)%A_richardson_p
contacts(ict)%ifbl = reg_ct(ri)%ifbl
contacts(ict)%tunneling = reg_ct(ri)%tunneling
contacts(ict)%m_tunnel_n = reg_ct(ri)%m_tunnel_n
contacts(ict)%m_tunnel_p = reg_ct(ri)%m_tunnel_p

! Calculate phims
call contacts(ict)%set_phims_schottky(smc)
```

### 3. Continuity Init (continuity.f90:118-143)
```fortran
if (has_schottky) then
  ! Store normal directions and thermionic velocity
  allocate(this%schottky_normal(par%nct))
  allocate(this%v_surf(par%nct))
  do ict = 1, par%nct
    if (par%contacts(ict)%type == CT_SCHOTTKY) then
      this%schottky_normal(ict) = get_normal_dir(par, ict)
      this%v_surf(ict) = schottky_velocity(par, ci, ict)
    end if
  end do

  ! Store E-field selectors
  allocate(this%efield(par%g%dim))
  do dir = 1, par%g%dim
    call this%efield(dir)%init(efield(dir), ...)
  end do
end if
```

---

## KEY NUMERICAL FEATURES

1. **Stable supply function**: Uses `log1p(expm1(eta_m) * f_metal)` to avoid cancellation
2. **Safe field computation**: `F_safe = sqrt(F² + eps²)` prevents division by zero
3. **Adaptive integration**: `quad()` with rtol=1e-9 for tunneling integral
4. **Non-constant Jacobian**: Updated each Newton iteration via `set_matr(const=.false., nonconst=.true.)`
5. **Full Newton coupling**: Derivatives computed for all parameters

---

## EXAMPLE CONFIG

```ini
[contact]
  name = "SCHOTTKY"
  type = "schottky"
  x = 0.0, 0.0 : nm
  y = 0.0, 1000 : nm
  phi_b = 0.7 : eV
  A_richardson_n = 112 : A/cm²/K²
  A_richardson_p = 32 : A/cm²/K²
  ifbl = true
  tunneling = true
  m_tunnel_n = 0.26
  m_tunnel_p = 0.39
```

---

## DEPENDENCIES

- `quad_m` - Adaptive quadrature with parameter derivatives
- `normalization_m` - Unit conversion (`norm`, `denorm`)
- `math_m` - `expm1`, `log1p`, `PI`
- `semiconductor_m` - `CR_ELEC`, `CR_HOLE`, `CR_CHARGE`, `get_dist`, `get_idist`

---

## SUMMARY: Required Changes for Another Branch

### Production Cleanup (optional - current branch works for testing):
1. **Line 379** in schottky.f90: Remove `phi_b = phi_b - norm(0.2, "eV")` (test artifact)
2. **Comment/remove debug prints** in schottky.f90, contact.f90, region.f90

### Verified - No Changes Needed:
- ✓ **Hole sign convention** in schottky_tunneling eta_m calculation is correct

### Files to Copy/Adapt:
1. `src/schottky.f90` - Core Schottky physics module
2. `src/contact.f90` - Contact data structure + phims calculation
3. `src/continuity.f90` - Robin BC application (Schottky-specific parts)
4. `src/region.f90` - Config parsing for Schottky parameters
5. `src/device_params.f90` - `get_ct_surf` function + contact initialization

### Dependencies Required:
- `quad_m` module with parameter derivative support
- `normalization_m` for `norm`, `denorm`
- `math_m` for `expm1`, `log1p`, `PI`
- `semiconductor_m` for `CR_ELEC`, `CR_HOLE`, `CR_CHARGE`, `get_dist`, `get_idist`
