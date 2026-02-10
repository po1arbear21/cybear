module schottky_m
  !! Schottky contact physics module
  !!
  !! Provides Robin BC helpers for thermionic emission and WKB tunneling
  !! at metal-semiconductor Schottky barriers.
  !!
  !! Robin BC formulation: J_schottky = v_surf * (n - n0B) + J_tn
  !!   - v_surf: thermionic velocity (Richardson equation)
  !!   - n0B: equilibrium injection density (with optional IFBL)
  !!   - J_tn: sub-barrier tunneling current (Tsu-Esaki integral)

  use device_params_m,  only: device_params
  use semiconductor_m,  only: CR_ELEC, CR_HOLE
  use math_m,           only: PI, expm1, log1p
  use quad_m,           only: quad
  use normalization_m,  only: norm, denorm

  implicit none
  private

  public :: get_normal_dir
  public :: schottky_velocity
  public :: schottky_n0b
  public :: schottky_tunneling

contains

  !============================================================================
  ! Get Normal Direction for a Contact
  !============================================================================
  function get_normal_dir(par, ict) result(dir)
    !! Determine which coordinate direction is normal to a contact surface.
    !! The normal direction is the one with minimum coordinate range across
    !! the contact vertices (i.e., the flat direction).

    type(device_params), intent(in) :: par
    integer, intent(in) :: ict
    integer :: dir

    integer :: i, d, n
    integer, allocatable :: idx(:)
    real :: cmin(3), cmax(3), range(3), pos(3)

    dir = 1
    cmin = 1e30
    cmax = -1e30

    allocate(idx(par%g%idx_dim))
    n = par%transport_vct(ict)%n
    if (n == 0) return

    do i = 1, n
      idx = par%transport_vct(ict)%get_idx(i)
      call par%g%get_vertex(idx, pos(1:par%g%dim))
      do d = 1, par%g%dim
        cmin(d) = min(cmin(d), pos(d))
        cmax(d) = max(cmax(d), pos(d))
      end do
    end do

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
    !! Calculate thermionic injection velocity for Robin BC.
    !!
    !! v_th = A* T^2 / (q N_c)
    !!
    !! This is the surface recombination velocity coefficient in the
    !! Robin BC: J_TE = v_surf * (n - n0B)

    type(device_params), intent(in) :: par
    integer, intent(in) :: ci   !! carrier index (CR_ELEC/CR_HOLE)
    integer, intent(in) :: ict  !! contact index
    real :: v_th

    real :: A_star

    if (ci == CR_ELEC) then
      A_star = par%contacts(ict)%A_richardson_n
    else
      A_star = par%contacts(ict)%A_richardson_p
    end if

    v_th = A_star * norm(par%T, "K")**2 / par%smc%edos(ci)
  end function schottky_velocity

  !============================================================================
  ! Robin BC Helper: Equilibrium Injection Density
  !============================================================================
  subroutine schottky_n0b(par, ci, ict, E_normal, n0B, dn0B_dE)
    !! Calculate equilibrium injection density at Schottky barrier.
    !!
    !! n0B = N_c exp(-phi_b_eff / kT)
    !!
    !! With optional image force barrier lowering (IFBL):
    !!   phi_b_eff = phi_b - sqrt(q|E| / (4 pi epsilon))

    type(device_params), intent(in) :: par
    integer, intent(in) :: ci, ict
    real, intent(in) :: E_normal
    real, intent(out) :: n0B
    real, optional, intent(out) :: dn0B_dE

    real :: phi_b, phi_b_eff, delta_phi, ddelta_dE

    if (ci == CR_ELEC) then
      phi_b = par%contacts(ict)%phi_b
    else
      phi_b = par%smc%band_gap - par%contacts(ict)%phi_b
    end if

    phi_b_eff = phi_b

    if (par%contacts(ict)%ifbl .and. abs(E_normal) > 1e-10) then
      delta_phi = sqrt(abs(E_normal) / (4.0 * PI))
      phi_b_eff = phi_b - delta_phi
      if (present(dn0B_dE)) ddelta_dE = delta_phi / (2.0 * abs(E_normal))
    else
      delta_phi = 0.0
      if (present(dn0B_dE)) ddelta_dE = 0.0
    end if

    n0B = par%smc%edos(ci) * exp(-phi_b_eff)

    if (present(dn0B_dE)) then
      dn0B_dE = n0B * ddelta_dE
    end if
  end subroutine schottky_n0b

  !============================================================================
  ! Robin BC Helper: Tunneling Current (Below-Barrier Integral)
  !============================================================================
  subroutine schottky_tunneling(par, ci, ict, E_normal, dens, J_tn, dJ_tn_ddens)
    !! Calculate sub-barrier tunneling current via Tsu-Esaki integral.
    !!
    !! J_tn = (m* / 2 pi^2) * integral_0^{phi_b} T(E) * supply(E) dE
    !!
    !! The above-barrier thermionic component is handled separately
    !! by the Robin BC velocity term.

    type(device_params), intent(in) :: par
    integer, intent(in) :: ict, ci
    real, intent(in) :: E_normal, dens
    real, intent(out) :: J_tn
    real, intent(out) :: dJ_tn_ddens

    real :: phi_b, m_tunnel, eta, eta_m, detadF
    real :: params(4), integral, dIda, dIdb, dIdp(4), err
    real :: prefactor
    integer :: ncalls

    if (.not. par%contacts(ict)%tunneling) then
      J_tn = 0.0
      dJ_tn_ddens = 0.0
      return
    end if

    ! barrier height and tunneling mass
    if (ci == CR_ELEC) then
      phi_b = par%contacts(ict)%phi_b
      m_tunnel = par%contacts(ict)%m_tunnel_n
    else
      phi_b = par%smc%band_gap - par%contacts(ict)%phi_b
      m_tunnel = par%contacts(ict)%m_tunnel_p
    end if

    prefactor = m_tunnel / (2.0 * PI**2)

    ! apply IFBL if enabled
    if (par%contacts(ict)%ifbl .and. abs(E_normal) > 1e-10) then
      phi_b = phi_b - sqrt(abs(E_normal) / (4.0 * PI))
    end if

    ! compute eta from density via inverse distribution function
    call par%smc%get_inv_dist(dens / par%smc%edos(ci), eta, detadF)

    ! quasi-Fermi level relative to metal Fermi level
    eta_m = eta + phi_b

    ! integration parameters: [phi_b, |E|, m*, eta_m]
    params = [phi_b, E_normal, m_tunnel, eta_m]

    ! integrate from 0 to phi_b (tunneling regime only)
    call quad(tsu_esaki_integrand, 0.0, phi_b, params, integral, &
              dIda, dIdb, dIdp, rtol=1.0e-9, err=err, ncalls=ncalls)

    ! current density (negative = injection into semiconductor)
    J_tn = -prefactor * integral

    ! derivative dJ_tn/d(dens) via chain rule:
    ! dJ/d(dens) = dJ/d(eta_m) * d(eta_m)/d(eta) * d(eta)/d(dens)
    dJ_tn_ddens = -prefactor * dIdp(4) * detadF / par%smc%edos(ci)
  end subroutine schottky_tunneling

  !============================================================================
  ! Tsu-Esaki Integrand
  !============================================================================
  subroutine tsu_esaki_integrand(E, p, f, dfdE, dfdp)
    !! Integrand: f(E) = T(E) * supply(E)
    !! where supply = ln(1+exp(eta_m-E)) - ln(1+exp(-E))
    !! and T = WKB transmission coefficient

    real, intent(in)  :: E, p(:)
    real, intent(out) :: f, dfdE, dfdp(:)

    real :: phi_b, efield, m_star, eta_m
    real :: T, dT_dE, dT_dphi, dT_dF
    real :: N_supply, dN_dE, f_metal, f_semi

    phi_b  = p(1)
    efield = p(2)
    m_star = p(3)
    eta_m  = p(4)

    call wkb_transmission(E, phi_b, efield, m_star, T, dT_dE, dT_dphi, dT_dF)

    f_metal = 1.0 / (1.0 + exp(E))
    f_semi  = 1.0 / (1.0 + exp(E - eta_m))
    dN_dE   = f_metal - f_semi

    ! supply function: rewritten for numerical stability
    ! log1p(expm1(eta_m) * f_metal) avoids catastrophic cancellation
    N_supply = log1p(expm1(eta_m) * f_metal)

    f = T * N_supply
    dfdE = dT_dE * N_supply + T * dN_dE

    dfdp(1) = dT_dphi * N_supply
    dfdp(2) = dT_dF * N_supply
    dfdp(3) = 0.0
    dfdp(4) = T * f_semi
  end subroutine tsu_esaki_integrand

  !============================================================================
  ! WKB Transmission Coefficient
  !============================================================================
  pure subroutine wkb_transmission(E, phi_b, F, m_star, T, dT_dE, dT_dphi, dT_dF)
    !! Standard WKB transmission coefficient for triangular barrier.
    !!
    !! T = exp(-gamma)  for E < phi_b (tunneling)
    !! T = 1            for E >= phi_b (thermionic, above barrier)
    !!
    !! gamma = (4/3) * sqrt(2 m*) * (phi_b - E)^(3/2) / F

    real, intent(in)  :: E, phi_b, F, m_star
    real, intent(out) :: T, dT_dE, dT_dphi, dT_dF

    real :: gamma, dE, F_safe, coeff, exp_gamma
    real, parameter :: eps = 1e-10

    dE = phi_b - E
    F_safe = sqrt(F**2 + eps**2)
    coeff = (4.0/3.0) * sqrt(2.0 * m_star)

    if (dE < 0.0) then
      T = 1.0
      dT_dE = 0.0
      dT_dphi = 0.0
      dT_dF = 0.0
    else
      gamma = coeff * dE**1.5 / F_safe
      exp_gamma = exp(-gamma)
      T = exp_gamma
      dT_dE   =  exp_gamma * coeff * 1.5 * sqrt(dE) / F_safe
      dT_dphi = -dT_dE
      dT_dF   =  exp_gamma * coeff * dE**1.5 * F / (F_safe**3)
    end if
  end subroutine wkb_transmission

end module schottky_m
