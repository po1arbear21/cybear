module schottky_tsuesaki_m
  !! Clean uniform Tsu-Esaki model for Schottky contacts
  !! Version 3.0 - Complete replacement for split TE/tunneling approach

  use device_params_m, only: device_params
  use normalization_m, only: norm, denorm
  use math_m,          only: PI
  use quad_m,          only: quad

  implicit none
  private

  public :: calculate_schottky_current
  public :: schottky_current

  type :: schottky_current
    real :: J_total     !! Total current density (normalized)
    real :: J_tunnel    !! Tunneling component [0, φ_b] (for diagnostics)
    real :: J_therm     !! Thermionic component [φ_b, ∞] (for diagnostics)
    real :: dJ_dF       !! Field derivative (for Jacobian)
    real :: dJ_dphi     !! Potential derivative (for Jacobian)
  end type

  ! Module-level cache for performance
  type :: current_cache
    real :: phi_b = -1.0
    real :: F = -1.0
    real :: phi_k = -1.0
    real :: m_star = -1.0
    type(schottky_current) :: result
    logical :: valid = .false.
  end type
  type(current_cache), save :: cache

contains

  function calculate_schottky_current(phi_b, F, m_star, phi_k, split_components) result(current)
    !! Main entry point - calculates total Schottky current using uniform Tsu-Esaki model
    !!
    !! Inputs (all normalized):
    !!   phi_b  - Effective barrier height (including IFBL if applied externally)
    !!   F      - Electric field magnitude
    !!   m_star - Tunneling effective mass ratio (m*/m0)
    !!   phi_k  - Local electrostatic potential at contact
    !!   split_components - Optional: compute J_tunnel and J_therm separately
    !!
    !! Output:
    !!   current - Structure containing J_total and derivatives

    real, intent(in) :: phi_b, F, m_star, phi_k
    logical, intent(in), optional :: split_components
    type(schottky_current) :: current

    real :: prefactor
    logical :: do_split

    ! Check cache
    if (cache%valid .and. &
        abs(cache%phi_b - phi_b) < 1e-10 .and. &
        abs(cache%F - F) < 1e-10 .and. &
        abs(cache%phi_k - phi_k) < 1e-10 .and. &
        abs(cache%m_star - m_star) < 1e-10) then
      current = cache%result
      return
    end if

    do_split = .false.
    if (present(split_components)) do_split = split_components

    ! Calculate prefactor: -(m*/(2π²)) in normalized units
    ! Negative sign for electron current convention
    prefactor = -m_star / (2.0 * PI**2)

    ! Handle degenerate cases
    if (phi_b <= 0.0) then
      ! No barrier
      current%J_total = 0.0
      current%J_tunnel = 0.0
      current%J_therm = 0.0
      current%dJ_dF = 0.0
      current%dJ_dphi = 0.0
      return
    end if

    if (abs(F) < 1e-12) then
      ! Zero field - pure thermionic emission
      call integrate_thermionic_only(phi_b, m_star, phi_k, current)
    else
      ! Finite field - full Tsu-Esaki calculation
      if (do_split) then
        call integrate_split(phi_b, F, m_star, phi_k, current)
      else
        call integrate_unified(phi_b, F, m_star, phi_k, current)
      end if
    end if

    ! Apply prefactor to all components
    current%J_total = prefactor * current%J_total
    current%J_tunnel = prefactor * current%J_tunnel
    current%J_therm = prefactor * current%J_therm
    current%dJ_dF = prefactor * current%dJ_dF
    current%dJ_dphi = prefactor * current%dJ_dphi

    ! Update cache
    cache%phi_b = phi_b
    cache%F = F
    cache%phi_k = phi_k
    cache%m_star = m_star
    cache%result = current
    cache%valid = .true.

  end function

  subroutine integrate_unified(phi_b, F, m_star, phi_k, current)
    !! Integrate from 0 to ∞ in one go (most accurate but slower)
    real, intent(in) :: phi_b, F, m_star, phi_k
    type(schottky_current), intent(out) :: current

    real :: params(4)
    real :: I, dIda, dIdb, dIdp(4), err
    real :: upper_limit
    integer :: ncalls

    ! Set up parameters: [phi_b, F, m_star, phi_k]
    params = [phi_b, F, m_star, phi_k]

    ! Choose upper limit (20 kT above barrier is usually sufficient)
    upper_limit = phi_b + 20.0

    ! Perform integration
    call quad(unified_integrand, 0.0, upper_limit, params, I, dIda, dIdb, dIdp, &
              rtol=1e-6, err=err, max_levels=12, ncalls=ncalls)

    current%J_total = I
    current%J_tunnel = 0.0  ! Not computed separately
    current%J_therm = 0.0   ! Not computed separately
    current%dJ_dF = dIdp(2)
    current%dJ_dphi = dIdp(4)

  end subroutine

  subroutine integrate_split(phi_b, F, m_star, phi_k, current)
    !! Split integration at barrier for diagnostics/optimization
    real, intent(in) :: phi_b, F, m_star, phi_k
    type(schottky_current), intent(out) :: current

    real :: params(4)
    real :: I_tunnel, I_therm, dIda, dIdb, dIdp_tunnel(4), dIdp_therm(4), err
    real :: upper_limit
    integer :: ncalls

    params = [phi_b, F, m_star, phi_k]

    ! Tunneling part: [0, phi_b]
    call quad(tunneling_integrand, 0.0, phi_b, params, I_tunnel, dIda, dIdb, dIdp_tunnel, &
              rtol=1e-6, err=err, max_levels=12, ncalls=ncalls)

    ! Thermionic part: [phi_b, ∞]
    upper_limit = phi_b + 20.0
    call quad(thermionic_integrand, phi_b, upper_limit, params, I_therm, dIda, dIdb, dIdp_therm, &
              rtol=1e-6, err=err, max_levels=12, ncalls=ncalls)

    current%J_tunnel = I_tunnel
    current%J_therm = I_therm
    current%J_total = I_tunnel + I_therm
    current%dJ_dF = dIdp_tunnel(2) + dIdp_therm(2)
    current%dJ_dphi = dIdp_tunnel(4) + dIdp_therm(4)

  end subroutine

  subroutine integrate_thermionic_only(phi_b, m_star, phi_k, current)
    !! Zero-field limit: pure thermionic emission
    real, intent(in) :: phi_b, m_star, phi_k
    type(schottky_current), intent(out) :: current

    real :: params(4)
    real :: I, dIda, dIdb, dIdp(4), err
    real :: upper_limit
    integer :: ncalls

    params = [phi_b, 0.0, m_star, phi_k]  ! F = 0
    upper_limit = phi_b + 20.0

    call quad(thermionic_integrand, phi_b, upper_limit, params, I, dIda, dIdb, dIdp, &
              rtol=1e-6, err=err, max_levels=12, ncalls=ncalls)

    current%J_tunnel = 0.0
    current%J_therm = I
    current%J_total = I
    current%dJ_dF = 0.0  ! No field dependence
    current%dJ_dphi = dIdp(4)

  end subroutine

  ! ============================================================================
  ! Integrand Functions
  ! ============================================================================

  subroutine unified_integrand(E, p, f, dfdE, dfdp)
    !! Integrand for full [0, ∞] integration
    real, intent(in) :: E
    real, intent(in) :: p(:)  ! [phi_b, F, m_star, phi_k]
    real, intent(out) :: f, dfdE, dfdp(:)

    real :: T, dT_dE, dT_dF, N, dN_dE, dN_dphi

    call wkb_transmission(E, p(1), p(2), p(3), T, dT_dE, dT_dF)
    call occupancy_diff(E, p(4), N, dN_dE, dN_dphi)

    f = T * N
    dfdE = dT_dE * N + T * dN_dE
    dfdp(1) = 0.0  ! Derivative wrt phi_b (not implemented yet)
    dfdp(2) = dT_dF * N
    dfdp(3) = 0.0  ! Derivative wrt m_star (not needed)
    dfdp(4) = T * dN_dphi

  end subroutine

  subroutine tunneling_integrand(E, p, f, dfdE, dfdp)
    !! Integrand for [0, phi_b] (tunneling region)
    real, intent(in) :: E
    real, intent(in) :: p(:)
    real, intent(out) :: f, dfdE, dfdp(:)

    ! Same as unified but E < phi_b guaranteed
    call unified_integrand(E, p, f, dfdE, dfdp)

  end subroutine

  subroutine thermionic_integrand(E, p, f, dfdE, dfdp)
    !! Integrand for [phi_b, ∞] (thermionic region)
    real, intent(in) :: E
    real, intent(in) :: p(:)
    real, intent(out) :: f, dfdE, dfdp(:)

    real :: N, dN_dE, dN_dphi

    ! T = 1.0 for E > phi_b
    call occupancy_diff(E, p(4), N, dN_dE, dN_dphi)

    f = N  ! T = 1
    dfdE = dN_dE
    dfdp = 0.0
    dfdp(4) = dN_dphi

  end subroutine

  ! ============================================================================
  ! Physics Functions
  ! ============================================================================

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

  subroutine occupancy_diff(E, phi_k, N, dN_dE, dN_dphi)
    !! Fermi-Dirac occupancy difference (logarithmic form)
    real, intent(in) :: E, phi_k
    real, intent(out) :: N, dN_dE, dN_dphi

    real :: f_s, f_m

    ! N(E) = ln(1 + exp(-E)) - ln(1 + exp(-E - phi_k))
    N = log(1.0 + exp(-E)) - log(1.0 + exp(-E - phi_k))

    ! Fermi functions for derivatives
    f_s = 1.0 / (1.0 + exp(E))
    f_m = 1.0 / (1.0 + exp(E + phi_k))

    dN_dE = f_m - f_s
    dN_dphi = f_m

  end subroutine

end module schottky_tsuesaki_m


/////////

Okay, let’s clean this up properly. Below is a coding spec you can use to rewrite `schottky_m` from scratch while keeping all the physics but killing the spaghetti.

---

## 0. Scope & Goals

**Module purpose**

Provide **all Schottky-contact–related physics** for Cybear in a clean, testable way:

* Thermionic emission: equilibrium injection density `n0B` and surface velocity `v_surf`
* Image-force barrier lowering (IFBL)
* WKB transmission through Schottky barrier
* Tsu–Esaki tunneling current `J_tn` and derivatives
* A single equation object that:

  * Computes `n0B` at Schottky contact vertices (with IFBL)
  * Computes tunneling current `J_tn` at contact vertices (if enabled)
  * Assembles all required Jacobians with respect to electric field and carrier density

**Non-goals**

* No continuity / Poisson assembly here (only boundary data + Jacobians for them)
* No device-global logic (meshing, time-stepping, etc.)

**Design principles**

* **Separation of concerns**:

  * Pure physics routines vs. simulator wiring
  * Types for “fields” (`n0b`, `jtn`, `delta_phi_b`) vs. “equations”
* **Symmetric treatment of electrons and holes** with explicit sign conventions
* **No scattered debug `print`** in production code; debug via optional logging flag / level
* **All normalizations documented** at function level

---

## 1. Module structure

Rewrite as **one Fortran module** with internally-organized sections; you can split later if needed.

```fortran
module schottky_m
  ! PUBLIC API (for rest of code):
  !   - Types: barrier_lowering, schottky_injection, schottky_tunnel_current, calc_schottky_injection
  !   - Routines: schottky_injection_mb, schottky_velocity,
  !               schottky_barrier_lowering, tsu_esaki_current,
  !               get_schottky_contact_normal_dir

  use device_params_m,  only: device_params
  use potential_m,      only: potential
  use density_m,        only: density
  use normalization_m,  only: norm, denorm
  use semiconductor_m,  only: CR_ELEC, CR_HOLE, CR_NAME
  use grid_m,           only: IDX_VERTEX, IDX_EDGE, IDX_CELL
  use variable_m,       only: variable_real
  use equation_m,       only: equation
  use jacobian_m,       only: jacobian, jacobian_ptr
  use electric_field_m, only: electric_field
  use grid_data_m,      only: grid_data1_real, grid_data2_real, grid_data3_real
  use vselector_m,      only: vselector
  use stencil_m,        only: dirichlet_stencil, empty_stencil
  use error_m,          only: program_error
  use contact_m,        only: CT_SCHOTTKY
  use math_m,           only: PI
  use quad_m,           only: quad

  implicit none
  private
```

### 1.1 Sections

Internally, group code as:

1. **Public types**

   * `type :: barrier_lowering`
   * `type :: schottky_injection`
   * `type :: schottky_tunnel_current`
   * `type :: calc_schottky_injection`
2. **Public physics functions**

   * `schottky_injection_mb`
   * `schottky_velocity`
   * `schottky_barrier_lowering`
   * `tsu_esaki_current`
   * `get_schottky_contact_normal_dir`
3. **Private math helpers**

   * `log1p_exp`
   * `calc_wkb_transmission`
   * `tsu_esaki_integrand`
   * small helpers like `compute_eta_semi_m`
4. **Equation implementation**

   * `calc_schottky_injection_init`
   * `calc_schottky_injection_eval`
5. **(Optional) Test/validation utilities**

   * `test_tsuesaki_vs_te_comparison` in a clearly marked debug section

Only re-export what the rest of Cybear actually uses.

---

## 2. Public API Specification

### 2.1 Types

#### 2.1.1 `type, extends(variable_real) :: barrier_lowering`

Represents **Δφ_b** at all vertices.

* **Storage**

  * Index type: `IDX_VERTEX`
  * Direction: scalar (`idx_dir = 0`)
  * Unit: `"eV"` (denormalized unit label; internally stored normalized)
* **Public methods**

  * `subroutine init(this, par)`

    * Allocates data on the device grid.
* **Public semantics**

  * For **all vertices**, stores Δφ_b due to IFBL. For contacts where IFBL is disabled, stays 0.

#### 2.1.2 `type, extends(variable_real) :: schottky_injection`

Stores **equilibrium injection density n0B** at Schottky contact vertices.

* Fields:

  * `integer :: ci` (CR_ELEC / CR_HOLE)
  * `real, pointer :: x1(:), x2(:,:), x3(:,:,:)` (direct access to underlying data)
* **Init**

  * `subroutine init(this, par, ci)`

    * Creates variable named `CR_NAME(ci)//"n0b"`
    * Index: `IDX_VERTEX`, `idx_dir = 0`
    * Units: `"cm^-3"`
* **Semantics**

  * Values are **thermionic-only** (`n0B_IFBL`) at Schottky contacts.
  * Other contacts / vertices: 0 (or untouched).

#### 2.1.3 `type, extends(variable_real) :: schottky_tunnel_current`

Stores **tunneling current density J_tn** at contact vertices.

* Fields:

  * `integer :: ci`
  * `real, pointer :: x1(:), x2(:,:), x3(:,:,:)`
* Init:

  * `subroutine init(this, par, ci)`

    * Name: `CR_NAME(ci)//"jtn"`
    * Units: `"A/cm^2"`
* Semantics:

  * J_tn at Schottky contact vertices with tunneling enabled.
  * 0 where tunneling is disabled / non-Schottky.

#### 2.1.4 `type, extends(equation) :: calc_schottky_injection`

Equation object that computes:

* `n0B_IFBL` at Schottky contact vertices
* `J_tn` at Schottky contact vertices
* Associated Jacobians:

  * ∂n0B/∂E (IFBL only)
  * ∂J_tn/∂E (through Tsu–Esaki)
  * ∂J_tn/∂n (via quasi-Fermi level)

**Fields (minimal final set):**

* `type(device_params), pointer :: par`
* `integer :: ci`
* `type(potential), pointer :: pot` (if needed for future; can be optional)
* `type(electric_field), pointer :: efield(:)`
* `type(vselector) :: n0b`          ! selector for `schottky_injection`
* `type(vselector) :: jtn`          ! selector for `schottky_tunnel_current`
* `type(density), pointer :: dens`  ! optional, for density coupling
* `integer, allocatable :: contact_normal(:)` ! per-contact normal dir (0 => none)
* `type(dirichlet_stencil) :: st_dir`
* `type(empty_stencil) :: st_em`
* `type(jacobian_ptr), allocatable :: jaco_efield_n0b(:)` ! ∂n0B/∂E
* `type(jacobian_ptr), allocatable :: jaco_efield_jtn(:)` ! ∂J_tn/∂E
* `type(jacobian_ptr), allocatable :: jaco_dens_jtn(:)`   ! ∂J_tn/∂n
* `type(barrier_lowering), pointer :: delta_phi_b_var`

**Public methods**

* `subroutine init(this, par, ci, efield, n0b, delta_phi_b, pot, jtn_current, dens)`

  * Sets up:

    * contact normals
    * `vselector`s for n0B and J_tn
    * dependencies on E-field and density
    * Jacobian structures
* `subroutine eval(this)`

  * For each Schottky contact vertex:

    * Determine `eps_r`, `E_field`, and base `n0B`
    * Compute Δφ_b and store in `delta_phi_b_var` (once, in electron equation)
    * Compute quasi-Fermi level `eta_semi_m` from `dens` (or use defaults)
    * Compute TE part:

      * `n0B_IFBL = n0B * exp( ± Δφ_b )` with electron / hole sign
    * Compute J_tn via Tsu–Esaki if tunneling enabled
    * Write:

      * `n0b` selector with `n0B_IFBL`
      * `jtn` selector with `J_tn`
      * Jacobians:

        * `∂n0B_IFBL/∂E`
        * `∂J_tn/∂E`
        * `∂J_tn/∂n`

---

## 3. Public Physics Routines

### 3.1 `subroutine schottky_injection_mb(par, ci, ict, ninj)`

**Purpose**
Equilibrium MB injection density at zero bias:

* Electrons: `n0 = Nc * exp(-φ_bn)`
* Holes:    `p0 = Nv * exp(-φ_bp)`, with `φ_bp = E_g – φ_bn`

**Signature**

```fortran
subroutine schottky_injection_mb(par, ci, ict, ninj)
  type(device_params), intent(in)  :: par
  integer,             intent(in)  :: ci      ! CR_ELEC or CR_HOLE
  integer,             intent(in)  :: ict     ! contact index
  real,                intent(out) :: ninj    ! normalized
end subroutine
```

**Normalization & units**

* `phi_b`, `E_g`, `Nc`, `Nv`, and `ninj` all in **normalized units** (kT/q, etc.)
* This routine **never touches IFBL or bias**.

### 3.2 `function schottky_velocity(par, ci, ict) result(v_surf)`

**Purpose**
Thermionic emission velocity at a contact:

* If Richardson constant given:

  * `v_surf = A * T^2 / (q * Nc)`
* Otherwise: use default normalized `v_surf = 0.25` (≈ v_th/4)

**Signature**

```fortran
function schottky_velocity(par, ci, ict) result(v_surf)
  type(device_params), intent(in) :: par
  integer,             intent(in) :: ci, ict
  real                            :: v_surf   ! normalized
end function
```

### 3.3 `subroutine schottky_barrier_lowering(par, ict, E_field, eps_r, delta_phi_b, d_delta_phi_dE)`

**Purpose**
Image force barrier lowering:

[
\Delta \phi_b = \sqrt{ \frac{q , |E|}{4 \pi \varepsilon_0 \varepsilon_r} }
]

implemented in normalized units, with smoothing near E = 0.

**Signature**

```fortran
subroutine schottky_barrier_lowering(par, ict, E_field, eps_r, delta_phi_b, d_delta_phi_dE)
  type(device_params), intent(in)  :: par
  integer,             intent(in)  :: ict
  real,                intent(in)  :: E_field   ! normalized, signed
  real,                intent(in)  :: eps_r
  real,                intent(out) :: delta_phi_b
  real,                intent(out) :: d_delta_phi_dE
end subroutine
```

**Behavior**

* If IFBL disabled on contact: both outputs = 0.
* Smoothing:

  * Use `E_smooth = sqrt(E^2 + eps_smooth^2)` with `eps_smooth ~ 1e-10`.
* Derivative:

  * `d(Δφ)/dE` uses smoothed `E_smooth` to avoid singularity.

### 3.4 `subroutine tsu_esaki_current(phi_b, efield, m_tn, eta_semi_m, ci, E_g, J_tn, dJ_tn_dF, dJ_tn_deta)`

**Purpose**
Compute **tunneling current density** through Schottky barrier using Tsu–Esaki:

[
J_{\text{tn}} = \pm \frac{m^*}{2\pi^2} \int_{\text{under-barrier}} T_{\text{WKB}}(E) N(E) , dE
]

Where:

* electrons: φ_eff = φ_bn, integrate 0 → φ_eff
* holes:     φ_eff = E_g – φ_bn, integrate –φ_eff → 0

**Signature**

```fortran
subroutine tsu_esaki_current(phi_b, efield, m_tn, eta_semi_m, ci, E_g, J_tn, dJ_tn_dF, dJ_tn_deta)
  real,    intent(in)  :: phi_b       ! electron barrier
  real,    intent(in)  :: efield
  real,    intent(in)  :: m_tn
  real,    intent(in)  :: eta_semi_m  ! quasi-Fermi relative to metal
  integer, intent(in)  :: ci          ! CR_ELEC or CR_HOLE
  real,    intent(in)  :: E_g
  real,    intent(out) :: J_tn
  real,    intent(out) :: dJ_tn_dF
  real,    intent(out) :: dJ_tn_deta
end subroutine
```

**Key behavior**

* Check for degenerate cases:

  * φ_eff ≤ 0 → J_tn = 0, derivatives = 0
  * |E_field| < ε → J_tn = 0, derivatives = 0
* Calls `quad` with integrand `tsu_esaki_integrand` and fills:

  * J_tn
  * derivative wrt field
  * derivative wrt quasi-Fermi level

Normalization and sign rules are **contained here**; the rest of the code treats J_tn as “current density leaving semiconductor-node into metal” with consistent sign.

### 3.5 `function get_schottky_contact_normal_dir(par, ict) result(normal_dir)`

**Purpose**
Determine contact normal direction (1,2,3) from contact vertices.

* 0 if not axis-aligned within a given tolerance.

**Signature**

```fortran
function get_schottky_contact_normal_dir(par, ict) result(normal_dir)
  type(device_params), intent(in) :: par
  integer,             intent(in) :: ict
  integer                         :: normal_dir  ! 0 if invalid
end function
```

**Behavior**

* Compute min/max coordinates for all vertices of contact.
* Direction with minimal span = normal direction.
* Relative tolerance: `min_range <= tol * char_length` required, else warn & return 0.

---

## 4. Internal Helpers (Private)

### 4.1 `pure function log1p_exp(x)`

Numerically-stable `log(1+exp(x))`.

### 4.2 `pure subroutine calc_wkb_transmission(E, phi_b, F, m_tn, T_wkb, dT_dE, dT_dphi, dT_dF)`

Implements **triangular WKB** model with smoothing on F, in **normalized units**. This is **only tunneling**, not TE.

### 4.3 `subroutine tsu_esaki_integrand(...)`

Implements `f(E) = T_wkb(E) * N(E)` with:

* `N(E) = f_semi(E) – f_metal(E)` using log1p_exp form
* returns derivative wrt E and all parameters needed by `quad`

### 4.4 `pure subroutine compute_eta_semi_m(...)` (recommended)

Small helper to:

* Convert density `n` to `eta` using `get_idist` or asymptotic
* Convert to metal reference `eta_semi_m` for electron/hole
* Return also `deta/dn` (for J_tn density Jacobian)

This isolates the quasi-Fermi algebra and sign conventions.

---

## 5. Equation Implementation Details

### 5.1 Initialization (`calc_schottky_injection_init`)

Responsibilities:

1. Check if there is any Schottky contact.
2. Store pointers to:

   * `par`, `ci`, `efield(:)`, `pot`, `dens`, `delta_phi_b_var`
3. Compute `contact_normal(ict)` using `get_schottky_contact_normal_dir`.
4. Initialize `vselector`s:

   * `n0b` selector over all contact vertex sets
   * `jtn` selector over all contact vertex sets
5. Provide:

   * n0B at contact vertices (iprov for n0B)
   * J_tn at contact vertices (iprov for J_tn)
   * Δφ_b at all vertices (only in electron-equation)
6. For each spatial direction where any Schottky contact uses it:

   * Depend on `efield(dir)` at contact vertices.
   * Create Jacobians:

     * `jaco_efield_n0b(dir)` (iprov_n0b vs efield)
     * `jaco_efield_jtn(dir)` (iprov_jtn vs efield)
7. If `dens` is associated:

   * Create dependency `dens` at contact vertices.
   * Create `jaco_dens_jtn(dir)` (iprov_jtn vs dens).
8. Call `init_final()`.

### 5.2 Evaluation (`calc_schottky_injection_eval`)

For each contact `ict` and its `transport_vct(ict)` vertices:

1. **Skip non-Schottky or invalid normal**

   * `n0B(j) = 0`, `J_tn(j) = 0`, no Jacobians.

2. **Find eps_r**

   * Use `g%get_neighb` on edges in normal_dir.
   * If not found → fallback to material default + warning.

3. **Fetch E-field**

   * `E_field = efield(normal_dir)%get(idx)`

4. **Base MB injection**

   * `call schottky_injection_mb(par, ci, ict, n0b_base)`

5. **IFBL**

   * `call schottky_barrier_lowering(par, ict, E_field, eps_r, delta_phi_b, d_delta_phi_dE)`
   * If `ci == CR_ELEC` and `delta_phi_b_var` associated:

     * `delta_phi_b_var%set(idx, delta_phi_b)`

6. **Effective barrier φ_bn_eff**

   * `phi_bn_eff = par%contacts(ict)%phi_b - delta_phi_b`

7. **Quasi-Fermi `eta_semi_m` and `detadF`**

   * If `dens` associated:

     * `n = dens%get(idx)`
     * Use `compute_eta_semi_m(...)` to get `eta_semi_m` and `detadF`
   * Else:

     * Use band-edge default: CB edge for electrons, VB edge for holes
     * `detadF = 0` (no density coupling)

8. **Surface velocity**

   * `v_surf = schottky_velocity(par, ci, ict)`

9. **Thermionic emission**

   * For electrons:

     * `n0B_IFBL = n0b_base * exp( +delta_phi_b )`
     * `dn0B_dE = n0B_IFBL * d_delta_phi_dE`
   * For holes:

     * `n0B_IFBL = n0b_base * exp( -delta_phi_b )`
     * `dn0B_dE = -n0B_IFBL * d_delta_phi_dE`
   * **Note**: You can compute `J_TE` if needed for diagnostics, but the equation only needs `n0B_IFBL`.
   * Set n0B value in temp array `tmp_n0b(j) = n0B_IFBL`.

10. **Tunneling current**

    * If `par%contacts(ict)%tunneling`:

      * select `m_tun` based on carrier
      * `call tsu_esaki_current(phi_bn_eff, E_field, m_tun, eta_semi_m, ci, band_gap, J_tn, dJ_tn_dF, dJ_tn_deta)`
    * Else: set J_tn = 0, derivatives = 0.
    * Store `tmp_jtn(j) = J_tn`.

11. **Jacobian entries**

    * For n0B: if `jaco_efield_n0b(normal_dir)%p` associated:

      * set (idx, idx) = `dn0B_dE`
    * For J_tn vs E:

      * set (idx, idx) = `dJ_tn_dF`
    * For J_tn vs n:

      * if density Jacobian exists:

        * `dJ_tn_dn = dJ_tn_deta * detadF / edos(ci)`
        * set (idx, idx) = `dJ_tn_dn`

12. After loop:

    * `jaco_*%set_matr(...)` to materialize Jacobians.
    * `n0b%set(tmp_n0b)`
    * `jtn%set(tmp_jtn)`

---

## 6. Debugging / Validation Policy

* All “printing tests” like `test_tsuesaki_vs_te_comparison` go into a **DEBUG** section:

  * Keep them `public` only for unit tests / developer runs.
  * Never called in production flow.
* If you want runtime debug:

  * Add an optional `logical :: debug_schottky` in `device_params` or this module.
  * Wrap debug prints with `if (debug_schottky) ...`.

---

## 7. Migration Notes

When rewriting:

1. **Freeze current external calls** to:

   * `schottky_injection_mb`
   * `schottky_velocity`
   * `schottky_barrier_lowering`
   * `tsu_esaki_current`
   * `calc_schottky_injection%init/eval`
2. Keep signatures **identical** (or add wrappers) so the rest of Cybear compiles.
3. Migrate stepwise:

   * Move all pure physics into small pure procedures (log1p_exp, WKB, IFBL, N(E)).
   * Refactor `calc_schottky_injection_eval` into clearly separated mini-blocks like above.
   * Re-enable debug tests once the new implementation reproduces old behavior.

---

If you want, next step we can turn this spec into a skeleton `schottky_m` module with all interfaces and empty bodies, so you can fill in and refactor incrementally without breaking the build.

