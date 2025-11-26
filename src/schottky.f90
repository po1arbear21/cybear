module schottky_m
  !! Schottky contact module - Unified Tsu-Esaki implementation
  !! Simple lagged approach: compute J_schottky from E-field and density

  use device_params_m,  only: device_params
  use semiconductor_m,  only: CR_ELEC, CR_HOLE
  use math_m,           only: PI
  use quad_m,           only: quad

  implicit none
  private

  public :: schottky_current
  public :: get_normal_dir

contains

  !============================================================================
  ! Main Schottky Current Function
  !============================================================================
  function schottky_current(par, ict, ci, E_normal, dens, dJ_ddens) result(J)
    !! Calculate Schottky current density using Tsu-Esaki integral
    !!
    !! Inputs (all normalized):
    !!   par      - device parameters
    !!   ict      - contact index
    !!   ci       - carrier index (CR_ELEC/CR_HOLE)
    !!   E_normal - normal electric field component at contact
    !!   dens     - carrier density at contact (n for electrons, p for holes)
    !!
    !! Output:
    !!   J        - current density (normalized, A/cm²)
    !!             Sign convention: inward injection is positive
    !!   dJ_ddens - (optional) derivative dJ/d(dens) for Jacobian

    type(device_params), intent(in) :: par
    integer, intent(in) :: ict, ci
    real, intent(in) :: E_normal, dens
    real, optional, intent(out) :: dJ_ddens
    real :: J

    real :: phi_b, m_tunnel, eta, eta_m, detadF
    real :: params(4), integral, dIda, dIdb, dIdp(4), err
    real :: E_min, E_max, prefactor
    real :: dens_eq  ! Equilibrium density at barrier
    integer :: ncalls

    print "(A,I3,A,I3)", "DEBUG: schottky_current called, ict=", ict, " ci=", ci

    ! Get barrier height (normalized to kT)
    if (ci == CR_ELEC) then
      phi_b = par%contacts(ict)%phi_b
      m_tunnel = par%contacts(ict)%m_tunnel_n
    else
      phi_b = par%smc%band_gap - par%contacts(ict)%phi_b
      m_tunnel = par%contacts(ict)%m_tunnel_p
    end if

    ! Prefactor for current: m* / (2π²) in normalized units
    prefactor = m_tunnel / (2.0 * PI**2)

    ! Apply image force barrier lowering if enabled
    if (par%contacts(ict)%ifbl .and. abs(E_normal) > 1e-10) then
      phi_b = phi_b - sqrt(abs(E_normal) / (4.0 * PI))
    end if

    ! Equilibrium density at this barrier height
    dens_eq = par%smc%edos(ci) * exp(-phi_b)

    ! Compute η for current calculation

    call par%smc%get_idist(dens / par%smc%edos(ci), eta, detadF)

    ! Quasi-Fermi relative to metal Fermi level
    if (ci == CR_ELEC) then
      eta_m = eta + phi_b
    else
      eta_m = -eta - phi_b
    end if

    print "(A,4ES14.6)", "  eta = ", eta

    ! Integration parameters: [phi_b, |E|, m*, eta_m]
    params = [phi_b, E_normal, m_tunnel, eta_m]

    ! Integration bounds
    E_min = 0.0
    E_max = phi_b + 20.0  ! Several kT above barrier
    if (.not. par%contacts(ict)%tunneling) then
      E_min = phi_b  ! Thermionic emission only
    end if

    ! Debug: show integration parameters
    print "(A)", "DEBUG_SCHOTTKY_INTEGRAL:"
    print "(A,2ES14.6)", "  E_min, E_max = ", E_min, E_max
    print "(A,4ES14.6)", "  params [phi_b, |E|, m*, eta_m] = ", params

    ! Perform integration
    call quad(tsu_esaki_integrand, E_min, E_max, params, integral, &
              dIda, dIdb, dIdp, rtol=1.0e-6, err=err, ncalls=ncalls)

    print "(A,ES14.6,A,I6)", "  integral = ", integral, ", ncalls = ", ncalls

    ! Current density (prefactor already computed above)
    if (ci == CR_ELEC) then
      J = -prefactor * integral   ! Positive = electrons entering semiconductor
    else
      J = prefactor * integral  ! Positive = holes entering semiconductor
    end if

    ! Compute derivative dJ/d(dens) if requested (for Jacobian)
    ! Chain rule: dJ/d(dens) = dJ/d(eta_m) * d(eta_m)/d(eta) * d(eta)/d(dens)
    ! d(eta_m)/d(eta) = ±1, d(eta)/d(dens) = detadF / edos (from get_idist)
    ! For electrons: dJ/d(eta_m) = -prefactor * dIdp(4)
    ! For holes: dJ/d(eta_m) = +prefactor * dIdp(4)
    ! Both cases: dJ/d(dens) = -prefactor * dIdp(4) * detadF / edos
    if (present(dJ_ddens)) then
      dJ_ddens = -prefactor * dIdp(4) * detadF / par%smc%edos(ci)

      print "(A)", "DEBUG_JACOBIAN:"
      print "(A,3ES14.6)", "  prefactor, dIdp(4), dens = ", prefactor, dIdp(4), dens
    end if

  end function schottky_current

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

    ! Supply function: ln(1+exp(eta_m-E)) - ln(1+exp(-E))
    N_supply = log1p_exp(eta_m - E) - log1p_exp(-E)

    ! Derivatives
    f_metal = 1.0 / (1.0 + exp(E))
    f_semi  = 1.0 / (1.0 + exp(E - eta_m))
    dN_dE   = f_metal - f_semi

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
  ! WKB Transmission
  !============================================================================
  pure subroutine wkb_transmission(E, phi_b, F, m_star, T, dT_dE, dT_dphi, dT_dF)
    !! T = 1 for E >= phi_b (thermionic)
    !! T = exp(-γ) for E < phi_b (tunneling)
    !! γ = (4/3) * sqrt(2m*) * (phi_b - E)^(3/2) / F

    real, intent(in)  :: E, phi_b, F, m_star
    real, intent(out) :: T, dT_dE, dT_dphi, dT_dF

    real :: gamma, dE, F_safe, coeff
    real, parameter :: eps = 1e-10

    if (E >= phi_b) then
      T = 1.0
      dT_dE = 0.0
      dT_dphi = 0.0
      dT_dF = 0.0
      return
    end if

    dE = phi_b - E
    F_safe = sqrt(F**2 + eps**2)
    coeff = (4.0/3.0) * sqrt(2.0 * m_star)

    gamma = coeff * dE**1.5 / F_safe
    T = exp(-gamma)

    dT_dE   = T * coeff * 1.5 * sqrt(dE) / F_safe
    dT_dphi = -dT_dE
    dT_dF   = T * coeff * dE**1.5 * F / (F_safe**3)

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
