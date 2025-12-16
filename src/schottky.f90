module schottky_m
  !! Schottky contact module - Unified Tsu-Esaki implementation
  !! Simple lagged approach: compute J_schottky from E-field and density

  use device_params_m,  only: device_params
  use semiconductor_m,  only: CR_ELEC, CR_HOLE
  use math_m,           only: PI, expm1, log1p
  use quad_m,           only: quad
  use normalization_m,  only: norm, denorm

  implicit none
  private

  ! public :: schottky_current
  public :: get_normal_dir
  public :: schottky_velocity
  public :: schottky_n0b
  public :: schottky_tunneling

contains

  ! !============================================================================
  ! ! Main Schottky Current Function
  ! !============================================================================
  ! function schottky_current(par, ict, ci, E_normal, dens, dJ_ddens) result(J)
  !   !! Calculate Schottky current density using Tsu-Esaki integral
  !   !!
  !   !! Inputs (all normalized):
  !   !!   par      - device parameters
  !   !!   ict      - contact index
  !   !!   ci       - carrier index (CR_ELEC/CR_HOLE)
  !   !!   E_normal - normal electric field component at contact
  !   !!   dens     - carrier density at contact (n for electrons, p for holes)
  !   !!
  !   !! Output:
  !   !!   J        - current density (normalized, A/cm²)
  !   !!             Sign convention: inward injection is positive
  !   !!   dJ_ddens - (optional) derivative dJ/d(dens) for Jacobian

  !   type(device_params), intent(in) :: par
  !   integer, intent(in) :: ict, ci
  !   real, intent(in) :: E_normal, dens
  !   real, optional, intent(out) :: dJ_ddens
  !   real :: J

  !   real :: phi_b, m_tunnel, eta, eta_m, detadF
  !   real :: params(4), integral, dIda, dIdb, dIdp(4), err
  !   real :: E_min, E_max, prefactor
  !   real :: dens_eq  ! Equilibrium density at barrier
  !   real :: integral1, integral2, err1, err2
  !   real :: dIda1, dIdb1, dIdp1(4), dIda2, dIdb2, dIdp2(4)
  !   integer :: ncalls, ncalls1, ncalls2

  !   ! print "(A,I3,A,I3)", "DEBUG: schottky_current called, ict=", ict, " ci=", ci

  !   ! Get barrier height (normalized to kT)
  !   if (ci == CR_ELEC) then
  !     phi_b = par%contacts(ict)%phi_b
  !     m_tunnel = par%contacts(ict)%m_tunnel_n
  !   else
  !     phi_b = par%smc%band_gap - par%contacts(ict)%phi_b
  !     m_tunnel = par%contacts(ict)%m_tunnel_p
  !   end if


  !   ! Prefactor for current: m* / (2π²) in normalized units
  !   ! prefactor = m_tunnel / (2.0 * PI**2)

  !   if (ci == CR_ELEC) then
  !     prefactor = m_tunnel / (2.0 * PI**2)

  !   else
  !     prefactor = m_tunnel / (2.0 * PI**2)
  !   end if


  !   ! Apply image force barrier lowering if enabled
  !   if (par%contacts(ict)%ifbl .and. abs(E_normal) > 1e-10) then
  !     phi_b = phi_b - sqrt(abs(E_normal) / (4.0 * PI))
  !   end if

  !   ! Equilibrium density at this barrier height
  !   dens_eq = par%smc%edos(ci) * exp(-phi_b)

  !   ! Compute η for current calculation

  !   call par%smc%get_idist(dens / par%smc%edos(ci), eta, detadF)


  !   ! Quasi-Fermi relative to metal Fermi level
  !   if (ci == CR_ELEC) then
  !     eta_m = eta + phi_b
  !   else
  !     eta_m = -eta - phi_b
  !   end if

  !   ! print "(A,4ES14.6)", "  eta = ", eta

  !   ! Integration parameters: [phi_b, |E|, m*, eta_m]
  !   params = [phi_b, E_normal, m_tunnel, eta_m]

  !   ! Integration bounds
  !   E_min = 0.0
  !   E_max = phi_b + 20.0  ! Several kT above barrier
  !   if (.not. par%contacts(ict)%tunneling) then
  !     E_min = phi_b  ! Thermionic emission only
  !   end if

  !   ! Debug: show integration parameters
  !   ! print "(A)", "DEBUG_SCHOTTKY_INTEGRAL:"
  !   ! print "(A,2ES14.6)", "  E_min, E_max = ", E_min, E_max
  !   ! print "(A,4ES14.6)", "  params [phi_b, |E|, m*, eta_m] = ", params

  !   ! Perform integration - split at barrier for better accuracy
  !   if (par%contacts(ict)%tunneling) then
  !     ! Split integration at barrier height
  !     ! Part 1: Below barrier (tunneling regime)
  !     call quad(tsu_esaki_integrand, 0.0, phi_b, params, integral1, &
  !               dIda1, dIdb1, dIdp1, rtol=1.0e-9, err=err1, ncalls=ncalls1)

  !     ! Part 2: Above barrier (thermionic regime)
  !     call quad(tsu_esaki_integrand, phi_b, E_max, params, integral2, &
  !               dIda2, dIdb2, dIdp2, rtol=1.0e-9, err=err2,  ncalls=ncalls2)
  !     ! Combine results
  !     integral = integral1 + integral2
  !     dIdp = dIdp1 + dIdp2
  !     err = sqrt(err1**2 + err2**2)
  !     ncalls = ncalls1 + ncalls2

  !     ! print "(A,2ES14.6,A,2I6)", "  integral = ", integral1, integral2, " ncalls = ", ncalls1, ncalls2
  !   else
  !     ! Single integration (small barrier or no tunneling)
  !     call quad(tsu_esaki_integrand, E_min, E_max, params, integral, &
  !               dIda, dIdb, dIdp, rtol=1.0e-6, err=err, ncalls=ncalls)

  !     ! print "(A,ES14.6,A,I6)", "  integral = ", integral, ", ncalls = ", ncalls
  !   end if

  !   ! ! Check integration error
  !   ! if (err > 1.0e-9) then
  !   !   print "(A,I2,A,ES12.4,A,I6)", &
  !   !     "  [SCHOTTKY WARNING] Contact ", ict, " err=", err, " ncalls=", ncalls
  !   ! end if

  !   ! Current density (prefactor already computed above) ******CR_CHARGE
  !   if (ci == CR_ELEC) then
  !     J = -prefactor * integral   ! Positive = electrons entering semiconductor
  !   else
  !     J = prefactor * integral  ! Positive = holes entering semiconductor
  !   end if

  !   ! Compute derivative dJ/d(dens) if requested (for Jacobian)
  !   ! Chain rule: dJ/d(dens) = dJ/d(eta_m) * d(eta_m)/d(eta) * d(eta)/d(dens)
  !   ! d(eta_m)/d(eta) = ±1, d(eta)/d(dens) = detadF / edos (from get_idist)
  !   ! For electrons: dJ/d(eta_m) = -prefactor * dIdp(4)
  !   ! For holes: dJ/d(eta_m) = +prefactor * dIdp(4)
  !   ! Both cases: dJ/d(dens) = -prefactor * dIdp(4) * detadF / edos
  !   if (present(dJ_ddens)) then
  !     dJ_ddens = -prefactor * dIdp(4) * detadF / par%smc%edos(ci)

  !     print "(A)", "DEBUG_JACOBIAN:"
  !     print "(A,3ES14.6)", "  prefactor, dIdp(4), dens = ", prefactor, dIdp(4), dens
  !   end if

  ! end function schottky_current

  !============================================================================
  ! Get Normal Direction for a Contact
  !============================================================================
  function get_normal_dir(par, ict) result(dir)
    !! Determine which direction is normal to a contact
    !! Returns: 1 for x-normal, 2 for y-normal, 3 for z-normal

    type(device_params), intent(in) :: par
    integer, intent(in) :: ict
    integer :: dir

    integer :: i, d, n
    integer, allocatable :: idx(:)
    real :: cmin(3), cmax(3), range(3), pos(3)

    dir = 1  ! Default
    cmin = 1e30
    cmax = -1e30

    allocate(idx(par%g%idx_dim))
    n = par%transport_vct(ict)%n
    if (n == 0) return

    ! Find coordinate ranges
    do i = 1, n
      idx = par%transport_vct(ict)%get_idx(i)
      call par%g%get_vertex(idx, pos(1:par%g%dim))
      do d = 1, par%g%dim
        cmin(d) = min(cmin(d), pos(d))
        cmax(d) = max(cmax(d), pos(d))
      end do
    end do

    ! Normal direction has minimum range
    range = cmax - cmin
    dir = 1
    do d = 2, par%g%dim
      if (range(d) < range(dir)) dir = d
    end do

  end function get_normal_dir

  !============================================================================
  ! Robin BC Helper: Thermionic Velocity (Richardson Velocity)
  !============================================================================
  function schottky_velocity(par, ci, ict) result(v_th)
    !! Calculate thermionic velocity for Robin BC
    !!
    !! v_th = A* T² / (q N_c)
    !!
    !! This is the "surface recombination velocity" for the thermionic
    !! emission component of Schottky BC.
    !!
    !! Returns: normalized velocity (cm/s in physical units)

    type(device_params), intent(in) :: par
    integer, intent(in) :: ci   ! carrier index (CR_ELEC/CR_HOLE)
    integer, intent(in) :: ict  ! contact index
    real :: v_th

    real :: A_star

    ! Get Richardson constant (A/cm²/K²)
    if (ci == CR_ELEC) then
      A_star = par%contacts(ict)%A_richardson_n
    else
      A_star = par%contacts(ict)%A_richardson_p
    end if


    ! Velocity = J_th / (q * N_c) = J_th / edos (in normalized units)
    ! First normalize J_th, then divide by edos
    v_th = A_star * norm(par%T, "K")**2 / par%smc%edos(ci)

    print "(A,I2,A,I2)", "DEBUG_SCHOTTKY_VELOCITY: ci=", ci, " ict=", ict
    print "(A,ES14.6,A)", "  A* = ", A_star, " A/cm²/K²"
    print "(A,ES14.6,A)", "  T = ", norm(par%T, "K"), " K"
    print "(A,ES14.6,A)", "  A*T² = ", A_star * norm(par%T, "K")**2, " A/cm²"
    print "(A,ES14.6,A)", "  v_th = ", denorm(v_th, "cm/s"), " cm/s"

  end function schottky_velocity

  !============================================================================
  ! Robin BC Helper: Equilibrium Injection Density
  !============================================================================
  subroutine schottky_n0b(par, ci, ict, E_normal, n0B, dn0B_dE)
    !! Calculate equilibrium injection density at Schottky barrier
    !!
    !! n0B = N_c exp(-φ_b_eff / kT)
    !!
    !! With optional image force barrier lowering (IFBL):
    !! φ_b_eff = φ_b - sqrt(q|E| / (4πε))
    !!
    !! Inputs:
    !!   par      - device parameters
    !!   ci       - carrier index
    !!   ict      - contact index
    !!   E_normal - normal electric field (normalized)
    !!
    !! Outputs:
    !!   n0B      - equilibrium density (normalized)
    !!   dn0B_dE  - (optional) derivative w.r.t. E-field for IFBL Jacobian

    type(device_params), intent(in) :: par
    integer, intent(in) :: ci, ict
    real, intent(in) :: E_normal
    real, intent(out) :: n0B
    real, optional, intent(out) :: dn0B_dE

    real :: phi_b, phi_b_eff, delta_phi, ddelta_dE

    ! Get barrier height (normalized to kT)
    if (ci == CR_ELEC) then
      phi_b = par%contacts(ict)%phi_b
    else
      phi_b = par%smc%band_gap - par%contacts(ict)%phi_b
    end if

    phi_b_eff = phi_b

    ! Apply image force barrier lowering if enabled
    if (par%contacts(ict)%ifbl) then
      delta_phi = sqrt(abs(E_normal) / (4.0 * PI))
      phi_b_eff = phi_b - delta_phi

      ! Derivative: d(delta_phi)/dE = 1/(2*sqrt(4πE)) = delta_phi/(2E)
      if (present(dn0B_dE)) then
        ddelta_dE = delta_phi / (2.0 * abs(E_normal))
      end if
    else
      delta_phi = 0.0
      if (present(dn0B_dE)) ddelta_dE = 0.0
    end if

    ! Equilibrium density at effective barrier
    n0B = par%smc%edos(ci) * exp(-phi_b_eff)

    ! Derivative w.r.t. E-field (for full Newton coupling)
    ! dn0B/dE = n0B * d(phi_b_eff)/dE = n0B * (-ddelta_dE)
    if (present(dn0B_dE)) then
      dn0B_dE = n0B * ddelta_dE
    end if

    ! print "(A,I2,A,I2)", "DEBUG_SCHOTTKY_N0B: ci=", ci, " ict=", ict
    ! print "(A,ES14.6)", "  phi_b = ", phi_b
    ! print "(A,ES14.6)", "  delta_phi (IFBL) = ", delta_phi
    ! print "(A,ES14.6)", "  phi_b_eff = ", phi_b_eff
    ! print "(A,ES14.6,A)", "  n0B = ", denorm(n0B, "cm^-3"), " cm^-3"

  end subroutine schottky_n0b

  !============================================================================
  ! Robin BC Helper: Tunneling Current (Below-Barrier Integral)
  !============================================================================
  subroutine schottky_tunneling(par, ci, ict, E_normal, dens, J_tn, dJ_tn_ddens)
    !! Calculate tunneling current from WKB integral below barrier
    !!
    !! J_tn = (m*/2π²) ∫₀^{φ_b} T(E) · supply(E) dE
    !!
    !! This is the sub-barrier tunneling component only.
    !! The above-barrier thermionic component is handled by Robin BC.
    !!
    !! Inputs:
    !!   par      - device parameters
    !!   ict      - contact index
    !!   ci       - carrier index
    !!   E_normal - normal electric field (normalized)
    !!   dens     - carrier density at contact (normalized)
    !!
    !! Outputs:
    !!   J_tn     - tunneling current density (normalized)
    !!   dJ_tn_ddens - derivative dJ_tn/d(dens) for Jacobian

    type(device_params), intent(in) :: par
    integer, intent(in) :: ict, ci
    real, intent(in) :: E_normal, dens
    real, intent(out) :: J_tn
    real, intent(out) :: dJ_tn_ddens

    real :: phi_b, m_tunnel, eta, eta_m, detadF
    real :: params(4), integral, dIda, dIdb, dIdp(4), err
    real :: prefactor
    integer :: ncalls

    ! Skip if tunneling disabled
    if (.not. par%contacts(ict)%tunneling) then
      J_tn = 0.0
      dJ_tn_ddens = 0.0
      return
    end if

    ! Get barrier height (normalized to kT)
    if (ci == CR_ELEC) then
      phi_b = par%contacts(ict)%phi_b
      m_tunnel = par%contacts(ict)%m_tunnel_n
    else
      phi_b = par%smc%band_gap - par%contacts(ict)%phi_b
      m_tunnel = par%contacts(ict)%m_tunnel_p
    end if

    ! Prefactor for current (same as in schottky_current)
    if (ci == CR_ELEC) then
      prefactor = m_tunnel / (2.0 * PI**2)
    else
      prefactor = m_tunnel / (2.0 * PI**2)
    end if

    ! Apply image force barrier lowering if enabled
    if (par%contacts(ict)%ifbl) then
      phi_b = phi_b - sqrt(abs(E_normal) / (4.0 * PI))
    end if

    phi_b = phi_b - norm(0.2, "eV")

    ! Compute η for current calculation
    call par%smc%get_idist(dens / par%smc%edos(ci), eta, detadF)

    ! print "(A,1ES14.6)", " eta = ", eta

    ! Quasi-Fermi relative to metal Fermi level
    if (ci == CR_ELEC) then
      eta_m = eta + phi_b
    else
      eta_m = eta + phi_b
    end if

    ! Integration parameters: [phi_b, |E|, m*, eta_m]
    params = [phi_b, E_normal, m_tunnel, eta_m]
    ! print "(A,4ES14.6)", "  params [phi_b, |E|, m*, eta_m] = ", params
    ! Integrate from 0 to phi_b (tunneling regime only)
    call quad(tsu_esaki_integrand, 0.0, phi_b, params, integral, &
              dIda, dIdb, dIdp, rtol=1.0e-9, err=err, ncalls=ncalls)

    ! print "(A,ES14.6)", "  err = ", err

    ! Current density
    if (ci == CR_ELEC) then
      J_tn = -prefactor * integral
    else
      J_tn = -prefactor * integral
    end if

    ! Derivative dJ_tn/d(dens)
    ! Same chain rule as in schottky_current
    dJ_tn_ddens = -prefactor * dIdp(4) * detadF / par%smc%edos(ci)

    ! print "(A,I2,A,I2)", "DEBUG_SCHOTTKY_TUNNELING: ci=", ci, " ict=", ict
    ! print "(A,ES14.6)", "  phi_b = ", phi_b
    ! print "(A,ES14.6)", "  integral = ", integral
    ! print "(A,ES14.6,A)", "  J_tn = ", denorm(J_tn, "A/cm^2"), " A/cm²"
    ! print "(A,ES14.6)", "  dJ_tn/ddens = ", dJ_tn_ddens

  end subroutine schottky_tunneling

  !============================================================================
  ! Tsu-Esaki Integrand
  !============================================================================
  subroutine tsu_esaki_integrand(E, p, f, dfdE, dfdp)
    !! f(E) = T(E) * [f_semi(E) - f_metal(E)]
    !! T = WKB transmission coefficient

    real, intent(in)  :: E, p(:)
    real, intent(out) :: f, dfdE, dfdp(:)

    real :: phi_b, efield, m_star, eta_m
    real :: T, dT_dE, dT_dphi, dT_dF
    real :: N_supply, dN_dE, f_metal, f_semi

    phi_b  = p(1)
    efield = p(2)
    m_star = p(3)
    eta_m  = p(4)

    ! WKB transmission
    call wkb_transmission(E, phi_b, efield, m_star, T, dT_dE, dT_dphi, dT_dF)

    ! Derivatives
    f_metal = 1.0 / (1.0 + exp(E))
    f_semi  = 1.0 / (1.0 + exp(E - eta_m))
    dN_dE   = f_metal - f_semi

    ! Supply function: ln(1+exp(eta_m-E)) - ln(1+exp(-E))
    ! Rewritten as log1p(expm1(eta_m) * f_metal) for numerical stability
    ! This avoids catastrophic cancellation and eliminates branching
    N_supply = log1p(expm1(eta_m) * f_metal)

    ! Integrand
    f = T * N_supply
    dfdE = dT_dE * N_supply + T * dN_dE

    ! Parameter derivatives
    dfdp(1) = dT_dphi * N_supply
    dfdp(2) = dT_dF * N_supply
    dfdp(3) = 0.0
    dfdp(4) = T * f_semi

  end subroutine tsu_esaki_integrand

  !============================================================================
  ! WKB Transmission (with smooth transition at barrier)
  !============================================================================
  pure subroutine wkb_transmission(E, phi_b, F, m_star, T, dT_dE, dT_dphi, dT_dF)
    !! Standard WKB transmission coefficient
    !! T = exp(-γ) for E < phi_b (tunneling)
    !! T = 1 for E >= phi_b (thermionic)
    !! γ = (4/3) * sqrt(2m*) * (phi_b - E)^(3/2) / F

    real, intent(in)  :: E, phi_b, F, m_star
    real, intent(out) :: T, dT_dE, dT_dphi, dT_dF

    real :: gamma, dE, F_safe, coeff, exp_gamma
    real, parameter :: eps = 1e-10

    dE = phi_b - E
    F_safe = sqrt(F**2 + eps**2)
    coeff = (4.0/3.0) * sqrt(2.0 * m_star)

    if (dE < 0.0) then
      ! Above the barrier: thermionic (E >= phi_b)
      T = 1.0
      dT_dE = 0.0
      dT_dphi = 0.0
      dT_dF = 0.0
    else
      ! Below barrier: tunneling (E < phi_b)
      gamma = coeff * dE**1.5 / F_safe
      exp_gamma = exp(-gamma)

      T = exp_gamma
      dT_dE   = exp_gamma * coeff * 1.5 * sqrt(dE) / F_safe
      dT_dphi = -dT_dE
      dT_dF   = exp_gamma * coeff * dE**1.5 * F / (F_safe**3)
    end if

  end subroutine wkb_transmission

  !============================================================================
  ! Stable log(1 + exp(x))
  !============================================================================
  pure function log1p_exp(x) result(y)
    real, intent(in) :: x
    real :: y

    if (x > 30.0) then
      y = x
    else if (x < -30.0) then
      y = exp(x)
    else
      y = log(1.0 + exp(x))
    end if
  end function log1p_exp

end module schottky_m
