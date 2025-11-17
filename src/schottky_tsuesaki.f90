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