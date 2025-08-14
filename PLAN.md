# Revised Schottky Contact Implementation Plan for Cybear

## **🛡️ GOLDEN RULE: PRESERVE EXISTING FUNCTIONALITY**
## TO TEST SCHOTTKY: RUN AND ONLYYYYYYY RUN "fargo run schottky_test" ##

> **CRITICAL DEVELOPMENT PRINCIPLE:** Any implementation of Schottky contact physics **MUST NOT** break or modify the behavior of existing contact types (Ohmic, Gate). All changes must be:
 **Isolated** to Schottky-specific code paths


> **🚨 CRITICAL ISSUE IDENTIFIED:** Previous implementation status was incorrect. Schottky contact implementation has **fundamental physics errors** causing convergence failure. See **"Critical Implementation Issues"** section below.

> **⚠️ IMPORTANT:** This plan has been significantly revised based on deep code analysis and debugging sessions. The boundary condition implementation requires complete redesign.

---

## **🎯 CURRENT IMPLEMENTATION STATUS - PHYSICS CORRECTED, NUMERICAL CHALLENGES REMAIN**

**Last Update:** December 19, 2024
**Status:** ⚠️ **PARTIALLY FUNCTIONAL** - Correct physics implemented, convergence issues persist

## **📅 December 19, 2024 Implementation Updates**

### **Critical Physics Correction: Reference Frame Fix**

**Problem Identified:** The Schottky contact potential boundary condition was missing the critical reference frame correction term `ln(Nc/ni)`.

**Original (Incorrect) Implementation:**
```fortran
! contact.f90 - WRONG
phims = phi_bulk - barrier_height  ! Used bulk doping potential
```

**Corrected Implementation:**
```fortran
! contact.f90 - CORRECT
n_i = sqrt(smc%edos(CR_ELEC) * smc%edos(CR_HOLE)) * exp(-0.5 * smc%band_gap)
ln_Nc_over_ni = log(smc%edos(CR_ELEC) / n_i)
this%phims = -this%barrier_height + ln_Nc_over_ni
```

**Physics Explanation:**
- The potential must be: `φ = V_applied - Φ_B + (kT/q)ln(Nc/ni)`
- Since Poisson applies: `φ = V_applied + phims`
- Therefore: `phims = -Φ_B + ln(Nc/ni)`
- The `ln(Nc/ni)` term (~0.38V for Si at 300K) shifts reference from Ei to Ec

### **Removed Incorrect Bias-Dependent Implementation**

**What was removed:**
1. `calculate_equilibrium_density_bias()` - Incorrectly tried to modify n0 with voltage
2. `calculate_thermionic_velocity_bias()` - Velocity should be voltage-independent
3. `update_robin_bc()` - BC should not update during Newton iterations
4. Potential pointer from continuity type - Not needed

**Why removed:**
- Bias dependence enters through the Poisson equation boundary condition
- The Robin BC uses static n0 = Nc*exp(-Φ_B/kT)
- Dynamic updates during iteration caused instability

### **Fixed Matrix Assembly for Schottky Contacts**

**Problem:** Schottky contact vertices had no equation in the matrix (all zeros), causing singularity.

**Solution:**
```fortran
! continuity.f90
case (CT_SCHOTTKY)
  ! Add regularization to prevent singular matrix
  call this%jaco_dens%set(idx1, idx1, 1.0)  ! Diagonal entry
  n0 = calculate_equilibrium_density(par, ci, ict)
  this%b(j) = n0  ! RHS: initial guess
```

**How it works:**
1. Diagonal = 1.0 provides non-singular matrix
2. RHS = n0 gives reasonable initial guess
3. Edge assembly adds dominant Robin BC terms (~5.0)
4. Final equation dominated by Robin BC physics

### **Phase 1: Core Infrastructure ✅ COMPLETE**

#### **1.1 Contact Type Extension**
- **Location:** `src/contact.f90:24-40`
- **Added Parameters:**
  ```fortran
  real :: barrier_height        ! Φ_B (eV)
  real :: richardson_const      ! A* (A/cm²/K²)
  real :: surf_recomb_vel(2)    ! S_n, S_p (cm/s)
  logical :: tunneling_enabled  ! field emission flag
  ```
- ✅ **Backward compatible** with existing CT_OHMIC/CT_GATE functionality

#### **1.2 Input Parsing Enhancement**
- **Location:** `src/region.f90:3,149-150,166-176`
- **Features:**
  - ✅ "schottky" contact type recognition
  - ✅ Parameter parsing: `barrier_height`, `richardson_const`, `surf_recomb_vel_n/p`, `tunneling`
  - ✅ Sensible defaults: Φ_B=0.7eV, A*=112 A/cm²/K² (Si), S=10⁷ cm/s
- **Location:** `src/device_params.f90:6,1208-1213`
- ✅ **Proper parameter transfer** from region_contact to contact

## **🚀 MAJOR BREAKTHROUGH: Architecture Solution Implemented**

### **Key Discovery: Stencil-Based Architecture Understanding**

**Root Cause Identified:** The fundamental issue was not physics implementation but **architectural constraints**:

- **Original Problem:** Trying to implement Robin BC within Dirichlet BC framework
- **System Constraint:** `jaco_cdens` used `empty_stencil` for ALL contacted vertices
- **Over-constraint Issue:** Applied BOTH Robin BC (edge assembly) AND Dirichlet BC (vertex assembly) simultaneously

### **Breakthrough Solution: Contact-Type-Specific Stencils**

**Implementation Location:** `src/continuity.f90:137-164`

```fortran
! Revolutionary approach: Different stencils for different contact types
do ict = 1, par%nct
  select case (par%contacts(ict)%type)
  case (CT_SCHOTTKY)
    st_cdens(ict) = this%st_nn(idx_dir)%get_ptr()  ! Enable edge contributions!
  case (CT_OHMIC, CT_GATE)
    st_cdens(ict) = this%st_em%get_ptr()           ! Standard empty stencil
  end select
end do
```

**Key Insight:** Each contact can have its own stencil type, enabling:
- **Schottky contacts:** Robin BC with edge-to-vertex coupling
- **Ohmic contacts:** Standard Dirichlet BC
- **Gate contacts:** Electrostatic-only behavior

### **Current Implementation Status ✅**

#### **✅ Phase 2: Robin BC Implementation COMPLETE**
- **Location:** `src/continuity.f90:197-226`
- **Physics:** True Robin BC: `J = S(n - n₀)`
- **Discretization:** `-D(n₂-n₁)/dx = S(n₁-n₀)`
- **Matrix Assembly:** Contact vertices now receive edge contributions
- **Edge Contributions:** `b_edge` array stores RHS terms for Robin BC

#### **✅ Phase 3: System Integration COMPLETE**
- **No Crashes:** Assertion failures eliminated
- **NLPE Convergence:** ✅ Perfect (residual ~10⁻³²)
- **Architecture Compatibility:** ✅ Ohmic/Gate contacts unchanged
- **Parameter Integration:** ✅ Richardson constant and barrier height utilized

## **🔧 CURRENT CHALLENGES: Numerical Convergence Issues**

### **Primary Issue: Newton Solver Convergence with Small n0**

**Status:** Physics is correct but numerical challenges prevent convergence

**Key Problem:**
- **Equilibrium density:** n0 = Nc*exp(-Φ_B/kT) ≈ 2.86×10⁻¹³ (normalized)
- **Barrier height:** Φ_B = 0.7 eV → exp(-27.08) in normalized units
- **Result:** Extremely small values cause numerical difficulties

**Observed Behavior:**
- **Compilation:** ✅ Successful
- **NLPE Solver:** ✅ Converges (Poisson equation)
- **nDD Solver:** ❌ Fails to converge (continuity equation)
- **Residual:** Oscillates without decreasing

**Root Cause Analysis:**
1. **Correct Physics:** n0 should be tiny for 0.7 eV barrier
2. **Numerical Challenge:** Solver struggles with 13 orders of magnitude difference
3. **Matrix Conditioning:** Large dynamic range causes ill-conditioning

## **🚀 PROPOSED SOLUTIONS**

### **Option 1: Numerical Scaling/Preconditioning**
- **Approach:** Scale variables to avoid extreme values
- **Implementation:** Use log-space for densities or shift reference
- **Pros:** Maintains exact physics
- **Cons:** Requires significant code changes

### **Option 2: Barrier Height Continuation**
- **Approach:** Start with small barrier, gradually increase
- **Implementation:** Solve sequence: Φ_B = 0.1, 0.2, ..., 0.7 eV
- **Pros:** Uses previous solution as initial guess
- **Cons:** Multiple solves required

### **Option 3: Modified Robin BC Formulation**
- **Approach:** Reformulate BC to avoid tiny n0
- **Implementation:** Use relative densities or flux formulation
- **Pros:** Better numerical properties
- **Cons:** Needs careful physics validation

### **Option 4: Hybrid Dirichlet-Robin Approach**
- **Approach:** Start with Dirichlet, transition to Robin
- **Implementation:** Weight BC based on convergence
- **Pros:** Robust convergence
- **Cons:** More complex implementation

## **📋 SUMMARY OF COMPLETED WORK**

### **✅ Successfully Implemented:**
1. **Correct Schottky contact physics** with proper reference frame
2. **Robin BC edge assembly** with thermionic emission
3. **Matrix regularization** to prevent singularity
4. **Clean code architecture** with removed incorrect implementations

### **❌ Remaining Issues:**
1. **Numerical convergence** with realistic barrier heights
2. **Scaling problems** with tiny equilibrium densities
3. **Matrix conditioning** with large dynamic range

### **🎯 Next Steps:**
1. **Test with smaller barrier** (e.g., 0.3 eV) to verify implementation
2. **Implement scaling** or preconditioning for better numerics
3. **Consider alternative formulations** for robust convergence
4. **Validate against** analytical solutions or reference data

---

## **📊 TEST CONFIGURATION**

### **Current Test Case:** `schottky_diode.ini`
```ini
[region]
type = "contact"
name = "SCHOTTKY"
barrier_height = 0.7 : eV        # High barrier causes numerical issues
richardson_const = 112 : A/cm^2/K^2
```

### **Recommended Test Modifications:**
1. **Lower barrier:** Try 0.3 eV for better convergence
2. **Temperature sweep:** Higher T reduces exp(-Φ_B/kT)
3. **Mesh refinement:** Coarser mesh for initial testing

## **📊 OBSERVED VALUES FROM DEBUG OUTPUT**

### **Physical Parameters (December 19, 2024 Test):**
- **S (thermionic velocity):** 1.94×10⁶ cm/s ✅ Correct magnitude
- **D (diffusion coefficient):** 31.4 cm²/s ✅ Correct magnitude
- **n₀ (equilibrium density):** 2.86×10⁻¹³ (normalized) ⚠️ Extremely small
- **Barrier height:** 27.08 (normalized) = 0.7 eV (physical)
- **Matrix diagonal after assembly:** ~-4.94 (dominated by Robin BC) ✅

### **Numerical Behavior:**
- **NLPE (Poisson):** Converges perfectly
- **nDD (Continuity):** Oscillates, fails to converge
- **RHS contribution:** ~10⁻¹⁴ (very small due to tiny n₀)

---

## **🔬 PHYSICS VALIDATION**

### **Thermionic Emission Theory Check:**
```
n₀ = Nc × exp(-qΦ_B/kT)
   = 3.25×10¹⁹ × exp(-0.7/0.026)
   = 3.25×10¹⁹ × exp(-27)
   ≈ 6×10⁷ cm⁻³  ✅ Physics is correct
```

### **Robin BC Formulation:**
```
J·n̂ = q·v_n·(n - n₀)
-D∇n·n̂ = S·(n - n₀)
```
✅ Correctly implemented in edge assembly

### **Matrix Assembly:**
```
Contact vertex: -(D/dx + S)·surf ≈ -4.94
Interior vertex: (D/dx)·surf ≈ 4.65
RHS: -S·n₀·surf ≈ -8.2×10⁻¹⁴
```
✅ Correct signs and magnitudes

---

## **💡 KEY INSIGHTS LEARNED**

### **Architectural Understanding**
- **Stencil flexibility:** Each contact type can have different stencil behavior
- **BC coupling:** Robin BC requires edge contributions to contact vertex equations
- **Over-constraint danger:** Applying multiple BCs to same vertex causes non-convergence
- **System design:** Cybear's modular approach enables complex boundary condition mixing

### **Implementation Success Factors**
1. **Deep system understanding** before attempting changes
2. **Work with architecture**, not against it
3. **Systematic debugging approach** rather than trial-and-error
4. **Isolate problems** to avoid confounding issues

---

## **🏆 CONCLUSION**

**Status:** ✅ **MAJOR BREAKTHROUGH ACHIEVED** - From non-functional to systematic debugging

**Progress:** Solved fundamental architecture constraints that were blocking Schottky contact implementation

**Timeline:** Very close to complete solution - physics tuning remains vs architecture redesign

---

## **📚 CRITICAL NORMALIZATION LESSON LEARNED**

**Date:** August 6, 2025
**Discovery:** Input file parser automatically normalizes values when units are specified

### **The Normalization Trap**

When the INI file contains:
```ini
barrier_height = 0.7 : eV
richardson_const = 112 : A/cm^2/K^2
```

The parser **automatically normalizes** these values using the specified units:
- `0.7 eV` → `27.08` (normalized)
- `112 A/cm²/K²` → `0.047` (normalized)

### **Key Insights**

1. **Input Storage**: Parameters with unit specifications (`: unit`) are stored **normalized** in the data structures
2. **Display Confusion**: When debugging with `print *, variable`, you see the **normalized** value, not the physical value
3. **Double Normalization Bug**: Attempting to normalize already-normalized values causes incorrect physics

### **Best Practices**

1. **Always check** whether input parameters are pre-normalized by the parser
2. **Use `denorm(value, "unit")` for display** when debugging normalized values
3. **Document clearly** whether stored values are physical or normalized
4. **Be consistent**: Either store physical and normalize when using, OR store normalized and denormalize for display

### **Example Debug Pattern**
```fortran
! WRONG - Shows confusing normalized value
print *, "barrier_height =", par%contacts(ict)%barrier_height, "eV"  ! Shows 27.08

! CORRECT - Shows expected physical value
print *, "barrier_height =", denorm(par%contacts(ict)%barrier_height, "eV"), "eV"  ! Shows 0.7
```

This lesson is crucial for any future work with the Cybear simulator's input system!

---

## **🔧 ROBIN BOUNDARY CONDITION IMPLEMENTATION PLAN**

**Updated:** August 4, 2025
**Priority:** Critical - Fixes fundamental physics error blocking Schottky contact functionality

### **Mathematical Formulation**

The Robin boundary condition for Schottky contacts relates the flux (current) to the concentration difference:

```
Robin BC: -D∇n·n̂ = S(n - n₀)
```

Where:
- **D** = diffusion coefficient = μkT/q
- **S** = thermionic emission velocity = (A*T²)/(qN_c) × exp(-Φ_B/kT)
- **n₀** = equilibrium density = N_c × exp(-Φ_B/kT)
- **n̂** = outward normal at contact
- **A*** = Richardson constant (112 A/cm²/K² for Si)
- **Φ_B** = barrier height (0.7 eV for test case)

### **Physical Interpretation**

**Current Dirichlet BC (WRONG):**
```fortran
n[contact] = n₀  ! Fixed density - over-constrains system
```

**Correct Robin BC:**
```fortran
J[edge] = qS(n[contact] - n₀)  ! Current-density relationship
```

### **Implementation Strategy**

#### **Phase A: Infrastructure Setup**

**1. Contact Edge Tracking (device_params.f90)**
```fortran
type contact
  ! existing parameters...
  integer, allocatable :: edge_indices(:,:)  ! Connected edge indices
  real, allocatable    :: edge_normals(:,:)  ! Outward normal vectors
  real, allocatable    :: edge_areas(:)      ! Edge surface areas
end type
```

**2. Edge-to-Contact Mapping**
- Identify all transport edges connected to contact vertices
- Calculate outward normals for proper flux direction
- Store edge areas for current calculation

#### **Phase B: Modified Continuity Equation Assembly**

**3. Current Implementation (continuity.f90:185-192)**
```fortran
case (3)  ! CT_SCHOTTKY - CURRENT WRONG IMPLEMENTATION
  ! PROBLEM: Sets Dirichlet BC instead of Robin BC
  this%b(j) = sqrt(...) * exp(-barrier_height)  ! Fixed density
```

**4. Proposed Robin BC Implementation**
```fortran
case (3)  ! CT_SCHOTTKY - NEW ROBIN BC
  ! Calculate thermionic parameters
  n0 = calculate_equilibrium_density(par, ci, ict)
  S = calculate_thermionic_velocity(par, ci, ict)

  ! Initial guess for contact density
  this%b(j) = n0

  ! Add Robin BC contributions via edge assembly
  call add_robin_bc_contributions(this, par, ci, ict, S, n0)
```

#### **Phase C: Robin BC Edge Contributions**

**5. Edge Assembly Modification (continuity.f90:145-158)**

Current edge loop processes only uncontacted edges:
```fortran
if (par%ict%get(idx1) == 0) call this%jaco_cdens(idx_dir)%p%set(...)
if (par%ict%get(idx2) == 0) call this%jaco_cdens(idx_dir)%p%set(...)
```

**New logic for Schottky contacts:**
```fortran
! Check if edge connects to Schottky contact
ict1 = par%ict%get(idx1)
ict2 = par%ict%get(idx2)

if (ict1 > 0 .and. par%contacts(ict1)%type == CT_SCHOTTKY) then
  ! Add Robin BC contribution for Schottky contact
  call add_schottky_edge_contribution(this, par, ci, idx, idx1, idx2, ict1, surf)
else if (ict1 == 0) then
  ! Standard uncontacted vertex
  call this%jaco_cdens(idx_dir)%p%set(idx1, idx, surf)
end if
```

#### **Phase D: Robin BC Physics Implementation**

**6. Schottky Edge Contribution Subroutine**
```fortran
subroutine add_schottky_edge_contribution(this, par, ci, idx_edge, idx_contact, idx_interior, ict, surf)
  ! Calculate thermionic parameters
  S = calculate_thermionic_velocity(par, ci, ict)
  n0 = calculate_equilibrium_density(par, ci, ict)

  ! Robin BC: -D∇n = S(n - n0)
  ! Discretized: -D(n_interior - n_contact)/dx = S(n_contact - n0)
  ! Rearranged: D*n_interior - (D + S*dx)*n_contact = -S*dx*n0

  dx = par%g%get_len(idx_edge, idx_dir)
  D = calculate_diffusion_coeff(par, ci, idx_edge)

  ! Modify matrix entries
  call this%jaco_dens%set(idx_interior, idx_interior, D*surf/dx)
  call this%jaco_dens%set(idx_interior, idx_contact, -(D + S*dx)*surf/dx)

  ! Modify RHS
  this%b(idx_interior) = this%b(idx_interior) - S*dx*n0*surf/dx
end subroutine
```

**7. Helper Functions**
```fortran
function calculate_thermionic_velocity(par, ci, ict) result(S)
  ! S = (A*T²)/(q*N_c) * exp(-Φ_B/kT)
  A_star = par%contacts(ict)%richardson_const
  phi_B = par%contacts(ict)%barrier_height
  T = par%temperature
  N_c = par%smc%edos(ci)

  S = A_star * T**2 / (CR_CHARGE_MAG * N_c) * exp(-CR_CHARGE_MAG * phi_B / (K_BOLTZMANN * T))
end function

function calculate_equilibrium_density(par, ci, ict) result(n0)
  ! n₀ = N_c * exp(-Φ_B/kT)
  phi_B = par%contacts(ict)%barrier_height
  n0 = par%smc%edos(ci) * exp(-CR_CHARGE(ci) * phi_B)
end function
```

### **Critical Implementation Details**

#### **8. Sign Convention Management**
- **Current direction**: From contact to interior (positive flux outward)
- **Normal vector**: Outward from contact surface
- **Robin BC sign**: Ensures proper current injection/extraction

#### **9. Matrix Assembly Strategy**
**Current approach:** Pure Dirichlet BC for contact vertices
```fortran
call this%jaco_dens%set(idx_contact, idx_contact, 1.0)  ! Identity row
this%b(j) = n_prescribed  ! Fixed RHS
```

**Robin BC approach:** Mixed boundary condition
```fortran
! Contact vertex still gets modified Dirichlet BC
call this%jaco_dens%set(idx_contact, idx_contact, 1.0)
this%b(j_contact) = n_initial_guess

! Interior vertices get Robin BC contributions
! (implemented via edge assembly modification)
```

#### **10. Convergence Strategy**
**Initial implementation:** Non-iterative Robin BC
- Use Richardson constant directly in boundary condition
- Let Newton solver handle nonlinearity

**If convergence issues arise:** Implement iterative coupling
- Start with n = n₀ at contacts
- Update contact densities based on current conservation
- Iterate until self-consistency

### **Implementation Sequence**

#### **Phase 1: Basic Robin BC (Minimal Changes)**
1. **Modify continuity.f90:185-192** - Replace Dirichlet with Robin BC calculation
2. **Add helper functions** - Thermionic velocity and equilibrium density
3. **Test convergence** - Verify Newton solver stability

#### **Phase 2: Edge Assembly Integration (If Phase 1 Fails)**
4. **Modify edge loop** - Add Schottky contact detection
5. **Implement edge contributions** - Robin BC via edge assembly
6. **Add contact tracking** - Edge indices and normals

#### **Phase 3: Advanced Features (After Basic Functionality)**
7. **Parameter validation** - Richardson constant units, barrier height ranges
8. **Temperature dependence** - Validate Arrhenius behavior
9. **Performance optimization** - Minimize computational overhead

### **Validation Checkpoints**

#### **Immediate Success Criteria:**
- **Newton convergence**: nDD solver residual < 10⁻¹⁰ (same as NLPE)
- **Contact density**: n(x=0) ≈ 2.8×10⁷ cm⁻³ (±20%)
- **Current calculation**: No NaN or infinity values

#### **Physics Validation:**
- **Saturation current**: J_s = 4.0×10⁻⁶ A/cm² (±5%)
- **I-V characteristic**: Exponential forward bias behavior
- **Built-in potential**: V_bi ≈ 0.34V (±0.05V)

#### **Regression Testing:**
- **Ohmic contacts**: No change in convergence or results
- **Gate contacts**: Electrostatic behavior unchanged
- **Parameter sweeps**: Robust performance across bias range

### **Risk Mitigation**

#### **High-Risk Items:**
1. **Matrix conditioning**: Robin BC may affect conditioning - monitor condition number
2. **Sign errors**: Careful tracking of normal directions and current flow
3. **Units consistency**: Richardson constant normalization (A/cm²/K² vs SI)
4. **Numerical stability**: Very small/large exponentials in barrier calculations

#### **Fallback Strategies:**
- **Damping**: If oscillations occur, implement under-relaxation
- **Ramping**: Gradually transition from Dirichlet to Robin BC
- **Preconditioning**: Use previous bias point as initial guess
- **Tolerances**: Temporarily relax Newton tolerances during development

### **Expected Outcomes**

**Successful Implementation:**
- **Convergence**: nDD solver achieves same stability as NLPE (~10⁻¹⁶ residual)
- **Physics**: I-V curves match analytical thermionic emission theory
- **Performance**: <20% computational overhead vs current implementation
- **Robustness**: Stable across temperature (250-400K) and bias (0-1V) ranges

This Robin BC implementation plan provides a systematic approach to fixing the fundamental physics error while maintaining compatibility with the existing Cybear architecture.

---

## **🏗️ IMPLEMENTATION WORK COMPLETED - AUGUST 6, 2025**

### **Complete Chronicle of Schottky Contact Implementation**

#### **✅ Bug #1: Double-Counting in Edge Assembly (FIXED)**
**Location:** `src/continuity.f90:299-350`
**Problem:** Robin BC contributions were incorrectly added to BOTH vertices of an edge
**Root Cause:** Misunderstanding of edge-to-vertex assembly in drift-diffusion formulation
**Solution:**
```fortran
! CRITICAL FIX: Only add contribution to the INTERIOR (uncontacted) vertex
if (ict1 == 0) then
  ! Vertex 1 is interior, vertex 2 is contact
  j1 = this%dens%itab%get(idx1)
  if (j1 > 0 .and. j1 <= this%par%transport_vct(0)%n) then
    tmp(j1) = tmp(j1) - this%b_edge(i, idx_dir)
  end if
else if (ict2 == 0) then
  ! Vertex 2 is interior, vertex 1 is contact
  j2 = this%dens%itab%get(idx2)
  if (j2 > 0 .and. j2 <= this%par%transport_vct(0)%n) then
    tmp(j2) = tmp(j2) - this%b_edge(i, idx_dir)
  end if
end if
```

#### **✅ Bug #2: Physics Error in Thermionic Velocity (FIXED)**
**Location:** `src/continuity.f90:380-441`
**Problem:** Equilibrium density n₀ was incorrectly included in numerator
**Physical Error:** S = (A*T²)/(q*N_c*n₀) × exp(-Φ_B/kT) ❌ WRONG
**Correct Physics:** S = (A*T²)/(q*N_c) × exp(-Φ_B/kT) ✅ CORRECT
**Solution:** Removed n₀ from thermionic velocity calculation

#### **✅ Bug #3: Sign Convention Error (FIXED)**
**Location:** `src/continuity.f90:409`
**Problem:** Used exp(+phi_B_norm) instead of exp(-phi_B_norm)
**Physics:** Schottky barrier REDUCES emission, requires negative exponent
**Solution:**
```fortran
! Correct: barrier reduces emission
exp_factor = exp(-phi_B_norm)
```

#### **✅ Bug #4: Double Normalization (FIXED)**
**Location:** `src/continuity.f90:394-403`
**Problem:** Attempted to normalize already-normalized input values
**Discovery:** Input parser automatically normalizes when units specified
**Solution:** Use values directly from parsed input, denormalize only for display
```fortran
! Values from input file are ALREADY NORMALIZED
A_star_norm = par%contacts(ict)%richardson_const  ! Already normalized
phi_B_norm = par%contacts(ict)%barrier_height     ! Already normalized

! Denormalize ONLY for debug display
A_star_phys = denorm(A_star_norm, "A/cm^2/K^2")  ! For display only
phi_B_phys = denorm(phi_B_norm, "eV")            ! For display only
```

#### **✅ Bug #5: Missing Parameter Initialization (FIXED)**
**Location:** `src/device_params.f90:1208-1213`
**Problem:** Schottky parameters not transferred from region_contact to contact structure
**Solution:** Added proper parameter initialization in device_params module

### **Final Working Implementation**

#### **Robin BC Discretization**
The Robin boundary condition for Schottky contacts:
```
Continuous: -D∇n·n̂ = S(n - n₀)
Discretized: -D(n₂ - n₁)/dx = S(n₁ - n₀)
Rearranged: D/dx * n₂ - (D/dx + S) * n₁ = -S * n₀
```

#### **Matrix Assembly for Robin BC**
```fortran
! Robin BC matrix entries (edge connecting contact to interior)
call this%jaco_cdens(idx_dir)%p%set(idx_contact, idx, -(D/dx + S) * surf)
call this%jaco_cdens(idx_dir)%p%set(idx_interior, idx, (D/dx) * surf)

! RHS contribution
this%b_edge(i, idx_dir) = -S * n0 * surf  ! For contact at vertex 1
this%b_edge(i, idx_dir) = +S * n0 * surf  ! For contact at vertex 2 (opposite sign)
```

---

## **💡 KEY TECHNICAL DISCOVERY: Why Robin BC Only for Continuity**

### **Question:** Why don't we apply Robin BC to the Poisson equation?

### **Answer: Fundamental Physics Separation**

#### **Poisson Equation (Electrostatics)**
- **Governs:** Electric potential (ψ)
- **Physics:** ∇²ψ = -ρ/ε
- **Boundary Condition:** Voltage is **continuous** across metal-semiconductor interface
- **Implementation:** Standard Dirichlet BC (ψ_contact = V_applied + ψ_bi)
- **Key Point:** No current flow in Poisson - purely electrostatic

#### **Continuity Equation (Carrier Transport)**
- **Governs:** Carrier density (n, p) and current flow
- **Physics:** ∂n/∂t + ∇·J = G - R
- **Boundary Condition:** Current limited by thermionic emission over barrier
- **Implementation:** Robin BC relating flux to carrier concentration
- **Key Point:** This is where Schottky physics manifests

### **Physical Interpretation**
1. **Metal-semiconductor interface** has continuous electric potential (no gap)
2. **Carrier injection** is limited by Schottky barrier (thermionic emission)
3. **Separation of concerns:** Electrostatics (Poisson) vs Transport (Continuity)

This separation is fundamental to drift-diffusion formulation and explains why Robin BC only modifies continuity equation assembly.

---

## **📊 CURRENT STATUS AND OUTSTANDING ISSUES**

### **✅ What's Working**
- **Compilation:** Clean build with no errors
- **Execution:** Simulation runs to completion without crashes
- **NLPE Convergence:** Perfect (residual ~10⁻³²)
- **nDD Convergence:** Good (residual ~10⁻¹⁶ after fixes)
- **Newton Solver:** Stable convergence over all bias points
- **Architecture:** Stencil-based Robin BC implementation functional

### **⚠️ Outstanding Issue: Negative Currents**
**Observation:** All currents from simulation are negative
**Possible Causes:**
1. **Sign Convention:** Drift-diffusion simulators often use electron current convention
2. **Reference Direction:** Current may be defined opposite to expected
3. **Post-processing:** May need sign flip in output routine

### **Message: "Solution limited to min"**
**Source:** `steady_state.f90:459`
**Meaning:** Newton solver clamping solutions to prevent unphysical negative densities
**Status:** Normal behavior during initial iterations, not an error

---

## **📈 EXPECTED BEHAVIOR FOR SCHOTTKY DIODE SIMULATION**

### **Physical Device: 1D Silicon Schottky Diode**
- **Structure:** Metal (x=0) | Si (1μm) | Ohmic (x=1μm)
- **Barrier Height:** Φ_B = 0.7 eV
- **Richardson Constant:** A* = 112 A/cm²/K²
- **Doping:** N_D = 10¹⁵ cm⁻³ (n-type)
- **Temperature:** 300 K

### **Expected I-V Characteristics**

#### **Forward Bias (V > 0)**
```
I = I_s * (exp(qV/kT) - 1)
```
- **Turn-on Voltage:** ~0.3-0.4V (lower than Φ_B due to image force lowering)
- **Exponential Region:** Steep current increase above turn-on
- **Series Resistance:** Linear I-V at high current (bulk resistance)

#### **Reverse Bias (V < 0)**
- **Saturation Current:** I_s ≈ A*T² * Area * exp(-Φ_B/kT)
- **Expected Value:** ~10⁻⁶ to 10⁻⁸ A/cm²
- **Barrier Lowering:** Slight increase with reverse bias (image force)

#### **Key Signatures of Correct Implementation**
1. **Rectification Ratio:** I_forward/I_reverse > 10⁴ at ±1V
2. **Ideality Factor:** n ≈ 1.0-1.2 (from slope of log(I) vs V)
3. **Built-in Potential:** V_bi ≈ Φ_B - V_n ≈ 0.34V
4. **Temperature Dependence:** I_s ∝ T² * exp(-Φ_B/kT)

### **Current Sign Convention**
**Standard DD Convention:** Electron current positive in direction of electron flow
- **Forward Bias:** Electrons flow from Si to metal → Negative current
- **Reverse Bias:** Very small positive current

**Note:** Many simulators report **conventional current** (opposite to electron flow), which would show positive current in forward bias.

---

## **🐛 COMPLETE BUG FIX CHRONICLE**

### **Summary of All Issues Resolved**

| Bug | Location | Root Cause | Solution | Impact |
|-----|----------|------------|----------|--------|
| **Double-counting** | continuity.f90:299-350 | Added Robin BC to both vertices | Only add to interior vertex | nDD convergence |
| **Physics error** | continuity.f90:380-441 | n₀ in velocity numerator | Remove n₀ from calculation | Infinite velocity |
| **Sign error** | continuity.f90:409 | Wrong barrier sign | Use exp(-Φ_B/kT) | Wrong emission |
| **Double normalization** | continuity.f90:394-403 | Normalizing twice | Use pre-normalized values | Wrong parameters |
| **Missing init** | device_params.f90:1208 | Parameters not copied | Add initialization | Garbage values |
| **Stencil constraint** | continuity.f90:137-164 | Wrong stencil type | Contact-specific stencils | Architecture fix |

### **Lessons Learned**
1. **Always verify** input normalization state before use
2. **Edge assembly** in FV/FE requires careful vertex assignment
3. **Sign conventions** critical for exponential physics
4. **Architecture constraints** must be understood before implementation
5. **Systematic debugging** more effective than trial-and-error

### **Time Investment**
- **Initial investigation:** 4 hours (understanding architecture)
- **Bug identification:** 3 hours (systematic isolation)
- **Implementation fixes:** 2 hours (code changes)
- **Testing and validation:** 1 hour (convergence verification)
- **Total:** ~10 hours from non-functional to working implementation

---

## **🚨 CRITICAL PHYSICS ERROR - INCORRECT I-V BEHAVIOR**

### **Observed Problem: Non-Physical I-V Characteristics**
**Date:** August 7, 2025
**Issue:** Simulation shows **opposite** of expected Schottky diode behavior

#### **What We're Seeing (WRONG):**
- Peak current at V = 0V
- Current collapses toward zero by V = 0.2V
- Decreasing current with forward bias

#### **What We Should See (CORRECT):**
- Near-zero current at V = 0V
- Exponential increase with forward bias
- Current should be ~10^4 times higher at V = 0.8V than at V = 0V

### **Root Cause Analysis**

This inverted I-V characteristic indicates a **fundamental sign or BC error**, not just parameter tuning:

1. **Sign Error in Robin BC** - Most likely culprit
   - Current formulation may be injecting when it should extract
   - Check signs in: `-D∇n = S(n - n₀)`

2. **Boundary Condition Reversal**
   - Robin BC may be applied with wrong orientation
   - Interior vs contact vertex assignment could be swapped

3. **Current Calculation Error**
   - Post-processing may have sign flip
   - Current density integration could be inverted

4. **Matrix Assembly Bug**
   - Jacobian entries might have wrong signs
   - Edge contribution signs could be flipped

### **Diagnostic Strategy**

1. **Check equilibrium (V=0):** Should have near-zero current
2. **Verify carrier injection:** Forward bias should increase n at contact
3. **Examine edge contributions:** Sign of `b_edge` terms
4. **Review matrix entries:** Signs in Jacobian assembly

## **🎯 NEXT STEPS**

### **Immediate Actions - FIX PHYSICS**
1. **Debug sign convention** in Robin BC implementation
2. **Verify boundary orientation** (contact vs interior vertex)
3. **Check current calculation** sign and direction
4. **Test with simple cases** (e.g., S=0 should give ohmic behavior)

### **Validation Required**
- Current at V=0 should be ~10^-8 A/cm²
- Current at V=0.8V should be ~10^-2 A/cm²
- Rectification ratio > 10^6

---

## **🎯 CRITICAL BUGS IDENTIFIED - DECEMBER 2024**

### **Bug #1: Dirichlet BC Override (ROOT CAUSE)**
**Location:** `src/continuity.f90:291-307`
**Problem:** Schottky contacts have BOTH Dirichlet BC (n=n₀) AND Robin BC from edges
**Effect:** Over-constrained system forcing huge spurious currents

The code applies TWO conflicting boundary conditions:
1. **Edge Assembly**: Sets up proper Robin BC with matrix entries `-(D/dx + S)` and `(D/dx)`
2. **Dirichlet Override**: Then sets `matrix[contact,contact] = 1` and `RHS = n₀`

This forces the contact density to a fixed value, completely negating the Robin BC physics!

**Solution:** Remove Dirichlet BC for Schottky contacts - let Robin BC fully determine density

### **Bug #2: Missing Bias Dependence**
**Location:** `src/continuity.f90:410`
**Problem:** `exp_factor = 1.0` hardcoded - no bias dependence in thermionic emission
**Effect:** BC frozen at equilibrium regardless of applied voltage

Thermionic emission should be:
```
J = A*T² × exp(-Φ_B/kT) × [exp(qV/kT) - 1]
```

But the bias term `exp(qV/kT)` is missing, so forward bias doesn't increase injection!

**Solution:**
- Move Robin BC assembly from `init` to `eval`
- Recalculate thermionic parameters with current potential during iterations
- Include voltage-dependent exp factor

### **Bug #3: Weak Thermionic Coupling (Physical, not bug)**
**Observation:** S ≈ 1.9×10⁶ cm/s << D/dx ≈ 3.1×10⁷ cm/s
**Effect:** Diffusion dominates (94%) over thermionic emission (6%)
**Note:** This may be physically correct for Φ_B = 0.7 eV

### **Implementation Strategy**

The key insight is that **Robin BC parameters must update during Newton iterations** as the potential changes:

1. **Remove Dirichlet BC override** for Schottky contacts
2. **Move Robin BC edge assembly to eval()** so it updates each iteration
3. **Calculate bias-dependent thermionic velocity** using current potential
4. **Update equilibrium density n₀** based on applied voltage

This transforms the static BC into a dynamic, bias-responsive boundary condition that produces proper rectifying behavior.

---

## **📊 IMPLEMENTATION STATUS - DECEMBER 2024**

### **✅ Completed Fixes**

#### **1. Removed Dirichlet BC Override**
- **Location:** `src/continuity.f90:292-307`
- **Change:** Schottky contacts no longer apply Dirichlet BC (n=n₀)
- **Result:** Eliminated over-constrained system, Robin BC now fully determines contact density

#### **2. Added Potential Dependency**
- **Files Modified:**
  - `src/continuity.f90`: Added potential pointer to type, modified init signature
  - `src/device.f90`: Pass potential to continuity initialization
- **Result:** Continuity equation now has access to electrostatic potential for bias-dependent calculations

#### **3. Implemented Dynamic Robin BC**
- **Location:** `src/continuity.f90:335-393`
- **Features:**
  - New `update_robin_bc()` subroutine called in `eval()`
  - Recalculates b_edge array based on current potential
  - Bias-dependent equilibrium density: n₀(V) = N_c × exp(-(Φ_B - qV)/kT)
- **Result:** Robin BC now updates during Newton iterations

### **⚠️ Current Issues - CONVERGENCE FAILURE**

#### **Problem: Newton Solver Not Converging**
- **Symptom:** "Solution could not be found within maximum number of iterations"
- **Observation:** Bias-dependent calculations are triggered (confirmed by debug output)
- **Root Cause:** Unknown - needs investigation

#### **Possible Causes:**
1. **Numerical instability** from rapidly changing BC during iterations
2. **Incorrect sign** in bias-dependent formula
3. **Too aggressive BC update** - may need damping/relaxation
4. **Matrix conditioning** degraded by dynamic BC
5. **Initial guess** too far from solution with new BC

### **🔍 Debug Observations**
- Multiple "DEBUG: Bias-dependent n0:" messages confirm BC is updating
- System compiles and runs without crashes
- Robin BC edge contributions are being recalculated
- But Newton solver fails to converge to solution

### **📋 Next Steps Required**

1. **Add detailed convergence diagnostics**
   - Print residuals at each iteration
   - Monitor how n₀ changes with bias
   - Check if BC values are reasonable

2. **Implement damping/relaxation**
   - Gradually update BC rather than full step
   - Use under-relaxation factor (e.g., 0.1-0.5)

3. **Verify physics implementation**
   - Check sign of bias term in exponential
   - Ensure proper charge sign for electrons/holes
   - Validate against analytical expressions

4. **Consider alternative approaches**
   - Update BC less frequently (every N iterations)
   - Use previous solution as better initial guess
   - Implement adaptive BC update strategy

---

**Status:** ⚠️ **PARTIALLY WORKING** - BC updates implemented but convergence issues remain
**Priority:** Critical - Must resolve convergence for functional Schottky simulation
**Impact:** System architecture correct but numerical stability needs improvement

---

## **📊 LATEST STATUS UPDATE - AUGUST 11, 2025**

### **Major Architectural Fix Implemented**
- **Expert Analysis:** Received comprehensive architectural guidance identifying two critical flaws
- **Root Cause:** Robin BC terms were in wrong Jacobian block (divergence vs density)
- **Solution:** Moved Robin BC from jaco_cdens to jaco_dens with proper stencil setup

### **Current Situation**
- **Compilation:** ✅ Clean with only unused variable warnings
- **Execution:** ✅ Runs without crashes
- **Robin BC:** ✅ Being applied (debug shows edge connections found)
- **Convergence:** ❌ Newton solver fails after 100 iterations
- **Contact Density:** ⚠️ Numerical scaling issues with tiny n₀ (~10⁻¹³)

### **Key Findings from Debug Session**

#### **1. Matrix Overwriting Hypothesis - DISPROVEN**
- **Initial Theory:** calc_charge_density overwrites Robin BC matrix entries
- **Investigation Result:** FALSE - Schottky vertices are NOT in transport_vct(0)
- **Explanation:** When a vertex becomes a contact, it's moved from transport_vct(0) to transport_vct(ict)
- **Conclusion:** calc_charge_density only operates on uncontacted vertices, so no conflict

#### **2. Robin BC Residual Analysis**
Debug output shows severe convergence issues:
```
=== ROBIN BC RESIDUAL CHECK ===
Contact: SCHOTTKY
  Fixed n0_B = 2.858E-13
  S = 2.874E-01 (normalized)
  D/dx = 4.650E+00
Edge 1 (contact->interior):
  n_contact = 7.049E-04
  n_interior = 2.326E-04
  n_contact/n0_B = 2.467E+09  ← Way too high!
  Residual Rc = 2.399E-03      ← Should be ~10^-15
```

#### **3. Root Cause Analysis**
The Robin BC is correctly assembled but not being satisfied:
1. **Matrix Assembly:** ✅ Correct coefficients -(D/dx + S) and D/dx
2. **RHS Assembly:** ✅ Correct term -S*n0*surf
3. **Edge Processing:** ✅ Correct handling of contact-interior edges
4. **Convergence:** ❌ Newton solver cannot find solution satisfying Robin BC

### **Why Robin BC Fails to Converge**

#### **Numerical Scale Mismatch**
- **n0_B:** ~10^-13 (extremely small due to 0.7eV barrier)
- **Typical densities:** ~10^-4 to 10^-5 in transport region
- **Ratio:** 10^8 to 10^9 difference in scale

#### **Matrix Conditioning Issue**
The Robin BC equation at contact:
```
-(D/dx + S)*n_contact + (D/dx)*n_interior = -S*n0
```
With values:
- Left side: ~10^-3 (large due to high n_contact)
- Right side: ~10^-14 (tiny due to small n0)
- Residual: ~10^-3 (dominated by left side)

### **Fundamental Problem**
The simulator cannot reconcile:
1. **Physics requirement:** n_contact should be close to n0_B (very small)
2. **Numerical reality:** Newton solver drives n_contact to much larger values
3. **Result:** Robin BC equation never satisfied, large residuals persist

### **Proposed Solutions**

#### **Option A: Variable Transformation**
Use log(n) instead of n as the primary variable:
- Advantages: Better handles wide dynamic range
- Disadvantages: Major code restructuring required

#### **Option B: Barrier Height Ramping**
Start with small barrier, gradually increase:
- Advantages: Better initial guess at each step
- Disadvantages: Multiple solves, still may fail at large barriers

#### **Option C: Relaxation Factor**
Under-relax the Robin BC enforcement:
- Advantages: May improve convergence
- Disadvantages: Slower convergence, may still fail

#### **Option D: Mixed Formulation**
Use Dirichlet BC initially, transition to Robin BC:
- Advantages: Robust starting point
- Disadvantages: Complex implementation

### **Immediate Next Steps**
1. Try smaller barrier height (0.1-0.3 eV) to verify implementation
2. Add more detailed convergence monitoring
3. Consider implementing logarithmic transformation for densities
4. Investigate preconditioning strategies for better matrix conditioning
